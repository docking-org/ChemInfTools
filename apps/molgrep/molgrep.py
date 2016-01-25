#!/usr/bin/env python
from __future__ import print_function

import argparse
import contextlib
from cStringIO import StringIO
import functools
import itertools
import logging
import os
import sys
import tempfile

from werkzeug import secure_filename
from flask import (
    abort,
    Flask, 
    redirect,
    render_template,
    request,
    Response,
    stream_with_context,
    url_for,
)

from rdkit.Chem import (
    EditableMol,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmiles,
    RDKFingerprint, 
    ReplaceSidechains,
    SmilesWriter,
)
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.DataStructs import (
    BulkTanimotoSimilarity, 
    BulkDiceSimilarity,
    BulkTverskySimilarity,
)

DEBUG = True


DESCRIPTORS = {
    'path': RDKFingerprint,
    'ecfp4': lambda mol: GetMorganFingerprintAsBitVect(mol, radius=2),
}


COEFFICIENTS = {
    'tanimoto': lambda x, ys, *args: BulkTanimotoSimilarity(x, ys),
    'dice': lambda x, ys, *args: BulkDiceSimilarity(x, ys),
    'tversky': lambda x, ys, a, b, *args: BulkTverskySimilarity(x, ys, a, b), 
}


def mol_parse(it, parser=MolFromSmiles):
    for num, line in enumerate(it, start=1):
        try:
            tokens = str(line).split()
            tokens.append(str(num))
            smiles, cid = tokens[:2]
            if not smiles:
                raise ValueError("No acceptable input to parse")
            mol = parser(smiles)
            if mol is not None:
                mol.SetProp('_Name', cid)
                yield mol
            else:
                raise ValueError("Parsing failed to yield a result")
        except Exception as e:
            cid = cid or ''
            logging.warning("Failed to parse/load ${:d}: {cid}. Reason: {!r}".format(num, cid, e))


def get_matching_parts(mol, *matches):
    indices = reversed(range(mol.GetNumAtoms()))
    for match in matches:
        match = set(match)
        part = EditableMol(mol)
        for idx in indices:
            if idx not in match:
                part.RemoveAtom(idx)
        yield part.GetMol()


def run_smarts_filter(needles, haystack, invert=False, annotate=False):
    for needle in needles:
        if annotate:
            for hay in haystack:
                matches = needle.GetSubstructMatches(hay, useChirality=True)
                if matches:
                    matching_parts = get_matching_parts(needle, *matches)
                    matching_smiles = [MolToSmiles(part, isomericSmiles=True) for part in matching_parts]
                    annotation = ';'.join(matching_smiles)
                    needle.SetProp('match', annotation)
                    yield needle
                elif invert:
                    needle.SetProp('match', annotaion)
                    yield needle
        elif not any(needle.HasSubstructMatch(hay, useChirality=True) for hay in haystack) == invert:
            yield needle


def run_similarity_filter(needles, haystack, coefficient, descriptor, threshold, invert=False, annotate=False):
    for needle_mol, needle_fp in needles:
        similarities = coefficient(needle_fp, haystack)
        matches = (similarity >= threshold for similarity in similarities)
        if not any(matches) == invert:
            if annotate:
                annotation = max(similarities)
                needle_mol.SetProp('match', '{:0.2f}'.format(annotation))
            yield needle_mol


def stream_smiles_results(operation, query, source, annotate=False):
    buf = StringIO()
    filtered = operation(query, source, annotate=annotate)
    results = SmilesWriter(buf, isomericSmiles=True, includeHeader=False)
    if annotate:
        results.SetProps(['match'])
    for result in filtered:
        buf.truncate(0)
        results.write(result)
        results.flush()
        yield buf.getvalue()
        

def write_smiles_results(operation, needles, haystack, dest, annotate=False):
    filtered = operation(needles, haystack, annotate=annotate)
    results = SmilesWriter(dest, isomericSmiles=True, includeHeader=False)
    if annotate:
        results.SetProps(['match'])
    for result in filtered:
        results.write(result)
    return results.NumMols()


def query_loader(parser, wrap=list, mode='r', reader_cls=argparse.FileType):
    reader = reader_cls(mode)

    def loader(path):
        source = reader(path)
        mols = mol_parse(source, parser=parser)
        mols = wrap(mols)
        return mols

    return loader


def smiles_reader(smiles, **kwargs):
    kwargs.setdefault('sanitize', True)
    return MolFromSmiles(smiles, **kwargs)


def smarts_reader(smarts, **kwargs):
    kwargs.setdefault('mergeHs', True)
    return MolFromSmarts(smarts, **kwargs)


def fp_reader(smiles, descriptor):
    mol = smiles_reader(smiles)
    fp = descriptor(mol)
    return mol


app = Flask(__name__)
app.config['DEBUG'] = DEBUG
app.config['UPLOAD_FOLDER'] = '/tmp'


@contextlib.contextmanager
def temporary_upload(uploaded):
    name = secure_filename(uploaded.filename)
    with tempfile.NamedTemporaryFile(dir=app.config.get('UPLOAD_FOLDER', '/tmp')) as tmp:
        uploaded.save(tmp.name)
        uploaded.seek(0)
        yield tmp


def process_upload(smarts, uploaded):
    with temporary_upload(uploaded) as source:
        results = stream_smiles_results(run_smarts_filter, smarts, source)
        for result in results:
            yield result


@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html')


@app.route('/filter', methods=['GET', 'POST'], endpoint='filter')
def filter_():
    smarts = request.form.get('pattern')
    smiles = request.files.get('source')
    if smarts and smiles:
        name = os.path.basename(smiles.filename)
        headers = {
            'Content-Disposition': 'attachment; filename={}'.format(name),
        }
        results = process_upload(smarts, smiles)
        try:
            results = itertools.chain([next(results)], results)
            code = 200
        except StopIteration:
            results = ()
            code = 404
        return Response(results, 
                        status=code,
                        headers=headers,
                        mimetype='chemical/x-daylight-smiles')
    else:
        return redirect(url_for('index'))


def main(params):
    source = params.input
    dest = params.output

    if params.command == 'smarts':
        query = params.smarts
        operation = functools.partial(run_smarts_filter,
                                      invert=params.invert)

    elif params.command == 'similarity':
        coefficient = COEFFICIENTS[params.coefficient]
        descriptor = DESCRIPTORS[params.descriptor]
        logging.info("Caching fingerprints")
        query = [descriptor(mol) for mol in params.smiles]  # Enumerate all
        source = ((mol, descriptor(mol)) for mol in source)  # Stream with original mol
        operation = functools.partial(run_similarity_filter,
                                      coefficient=coefficient,
                                      descriptor=descriptor,
                                      invert=params.invert,
                                      threshold=params.threshold)

    matched = write_smiles_results(operation=operation,
                                   needles=source, 
                                   haystack=query,
                                   dest=dest,
                                   annotate=params.annotate)
    return matched == 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest='command')
    commands.add_parser('serve')

    cli_substruct = commands.add_parser('smarts')
    cli_substruct.add_argument('-v', '--invert', dest='invert', action='store_true', default=False,
                               help='Only return those not matching all queries instead of those matching any')
    cli_substruct.add_argument('-a', '--annotate', dest='annotate', action='store_true', default=False,
                               help='Annotate results with match data')

    cli_substruct_input = cli_substruct.add_mutually_exclusive_group(required=True)
    cli_substruct_input.add_argument('-s', '--smarts', nargs='?',
                                     type=lambda smarts: [smarts_reader(smarts)],
                                     help='SMARTS or SMILES pattern to match')
    cli_substruct_input.add_argument('-f', '--file', dest='smarts', nargs='?',
                                     type=query_loader(smarts_reader),
                                     help='File containing smarts patterns to match')
    cli_substruct.add_argument('input', nargs='?', 
                               default=query_loader(smiles_reader, wrap=iter)('-'),
                               type=query_loader(smiles_reader),
                               help='Source SMILES to search [default: stdin]')
    cli_substruct.add_argument('output', nargs='?', default='-', type=str, 
                               help='Destination to write matching SMILES to [default: stdout]')

    cli_similarity = commands.add_parser('similarity')
    cli_similarity.add_argument('-v', '--invert', dest='invert', action='store_true', default=False,
                                help='Only return those not matching all queries instead of those matching any')
    cli_similarity.add_argument('-a', '--annotate', dest='annotate', action='store_true', default=False,
                               help='Annotate results with match data')
    cli_similarity.add_argument('-t', '--threshold', type=float, 
                                help='Similarity threshold to match')
    cli_similarity.add_argument('-c', '--coefficient', nargs='?', choices=COEFFICIENTS.keys(), default='tanimoto',
                                help='Similarity measure to use [default: %(default)s]')
    cli_similarity.add_argument('-d', '--descriptor', nargs='?', choices=DESCRIPTORS.keys(), default='path',
                                help='Fingerprint (descriptor) to use [default: %(default)s]')
    cli_similarity_input = cli_similarity.add_mutually_exclusive_group(required=True)
    cli_similarity_input.add_argument('-s', '--smiles', nargs='?',
                                      type=lambda smiles: [smiles_reader(smiles)],
                                      help='SMILES to generate fingerprints from')
    cli_similarity_input.add_argument('-f', '--file', dest='smiles', nargs='?',
                                      type=query_loader(smiles_reader), 
                                      help='File containing smiles patterns to generate finterprints from')
    cli_similarity.add_argument('input', nargs='?', 
                                default=query_loader(smiles_reader, wrap=iter)('-'),
                                type=query_loader(smiles_reader),
                                help='Source SMILES to search [default: stdin]')
    cli_similarity.add_argument('output', nargs='?', default='-', type=str, 
                               help='Destination to write matching SMILES to [default: stdout]')
    params = parser.parse_args()
    if params.command == 'serve':
        sys.exit(app.run(host='0.0.0.0', port=8081, debug=DEBUG))
    else:
        sys.exit(main(params))

