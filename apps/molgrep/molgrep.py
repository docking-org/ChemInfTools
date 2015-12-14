#!/usr/bin/env python
from __future__ import print_function

import argparse
import contextlib
from cStringIO import StringIO
import functools
import itertools
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
    RDKFingerprint, 
    MolFromSmarts,
    SmilesWriter,
    MolFromSmiles,
)
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.DataStructs import (
    TanimotoSimilarity, 
    DiceSimilarity,
    TverskySimilarity,
)

DEBUG = True


DESCRIPTORS = {
    'path': RDKFingerprint,
    'ecfp4': lambda mol: GetMorganFingerprintAsBitVect(mol, radius=2),
}


COEFFICIENTS = {
    'tanimoto': lambda x, y, *args: TanomotoSimilarity(x, y),
    'dice': lambda x, y, *args: DiceSimilarity(x, y),
# TODO: Add parsing of alpha and beta
#    'tversky': lambda x, y, a, b, *args: TverskySimilarity(x, y, a, b), 
}


def parse_smiles(it):
    for line in it:
        smiles, cid = str(line).split()[:2]
        try:
            mol = MolFromSmiles(smiles)
            if mol is not None:
                mol.SetProp('_Name', cid)
                yield mol
        except Exception as e:
            logging.warning("Failed to parse/load {cid}. Reason: {!r}".format(cid, e))


def run_smarts_filter(smarts, smiles):
    needle = MolFromSmarts(smarts)
    haystack = parse_smiles(smiles)
    return (mol for mol in haystack if mol.HasSubstructMatch(needle))


def run_similarity_filter(needle, haystack, coefficient, descriptor, threshold):
    query = descriptor(MolFromSmiles(needle))
    mols = parse_smiles(haystack)
    for mol in mols:
        fp = descriptor(mol)
        if coefficient(query, fp) >= threshold:
            yield mol


def stream_smiles_results(operation, query, source):
    buf = StringIO()
    filtered = operation(query, source)
    results = SmilesWriter(buf, isomericSmiles=True, includeHeader=False)
    for result in filtered:
        buf.truncate(0)
        results.write(result)
        results.flush()
        yield buf.getvalue()
        

def write_smiles_results(operation, query, source, dest):
    filtered = operation(query, source)
    results = SmilesWriter(dest, isomericSmiles=True, includeHeader=False)
    for result in filtered:
        results.write(result)
    return results.NumMols()



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
    if params.command == 'smarts':
        operation = run_smarts_filter
    elif params.command == 'similarity':
        coefficient = COEFFICIENTS[params.coefficeint]
        descriptor = DESCRIPTORS[params.descriptor]
        operation = functools.partial(run_similarity_filter,
                                      coefficient=coefficient,
                                      descriptor=descriptor,
                                      threshold=params.threshold)
    matched = write_smiles_results(operation=operation,
                                   smarts=params.smarts, 
                                   source=params.input, 
                                   dest=params.output)
    return matched == 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest='command')
    commands.add_parser('serve')
    cli_substruct = commands.add_parser('smarts')
    cli_substruct.add_argument('smarts', help='SMARTS or SMILES pattern to match')
    cli_substruct.add_argument('input', nargs='?', default=sys.stdin, type=argparse.FileType('r'),
                               help='Source SMILES to search [default: stdin]')
    cli_substruct.add_argument('output', nargs='?', default=sys.stdout, type=argparse.FileType('w'), 
                               help='Destination to write matching SMILES to [default: stdout]')
    cli_similarity = commands.add_parser('similarity')
    cli_similarity.add_argument('smiles', help='SMILES pattern to match')
    cli_similarity.add_argument('-t', '--threshold', type=float, 
                                help='Similarity threshold to match')
    cli_similarity.add_argument('-c', '--coefficient', nargs='?', choices=COEFFICIENTS.keys(), default='tanimoto',
                                help='Similarity measure to use [default: %(default)s]')
    cli_similarity.add_argument('-d', '--descriptor', nargs='?', choices=DESCRIPTORS.keys(), default='path',
                                help='Fingerprint (descriptor) to use [default: %(default)s]')
    cli_similarity.add_argument('input', nargs='?', default=sys.stdin, type=argparse.FileType('r'),
                                help='Source SMILES to search [default: stdin]')
    cli_similarity.add_argument('output', nargs='?', default=sys.stdout, type=argparse.FileType('w'), 
                               help='Destination to write matching SMILES to [default: stdout]')
    params = parser.parse_args()
    if params.command == 'serve':
        sys.exit(app.run(host='0.0.0.0', port=8081, debug=DEBUG))
    else:
        sys.exit(main(params))

