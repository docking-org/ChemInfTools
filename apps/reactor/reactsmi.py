#!/usr/bin/env python
from __future__ import print_function

import itertools
import sys

from rdkit.Chem import (
    MolFromSmiles,
    MolToSmiles,
)
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
try:
    from flask import (
        abort,
        Flask, 
        redirect,
        render_template,
        request,
        Response,
        url_for,
    )
except ImportError:
    Flask = None


DEBUG = True


def load_smiles_file(it):
    for line in it:
        smiles, cid = str(line).strip().split()[:2]
        mol = MolFromSmiles(smiles)
        if mol is not None:
            mol.SetProp('_Name', cid)
            yield mol


def reaction_matrix(rxn, *molsLists):
    for mols in itertools.product(*molsLists):
        products = rxn.RunReactants(mols)
        for product in products:
            yield mols, product


def smiles_reaction_matrix(smarts, *sources, **kwargs):
    sep = kwargs.setdefault('sep', ' ')
    reaction = ReactionFromSmarts(smarts)
    smilesLists = [load_smiles_file(source) for source in sources]
    products = reaction_matrix(reaction, *smilesLists)
    for reactants, product in products:
        cids = [r.GetProp("_Name") for r in reactants]
        product_id = '.'.join(cids)
        for mol in product:
            smiles = MolToSmiles(mol, isomericSmiles=True)
            yield sep.join((smiles, product_id))

if Flask is not None:
    app = Flask(__name__)
    app.config['DEBUG'] = DEBUG
    app.config['UPLOAD_FOLDER'] = '/tmp'
    
    
    @app.route('/')
    @app.route('/index')
    def index():
        return render_template('index.html')
    
    
    @app.route('/react', methods=['GET', 'POST'])
    def react():
        reaction = request.form.get('reaction')
        smiles1 = request.files.get('reactant1')
        smiles2 = request.files.get('reactant2')
        if reaction and smiles1 and smiles2:
            reaction = str(reaction)
            smiles1 = list(smiles1)
            smiles2 = list(smiles2)
            products = smiles_reaction_matrix(reaction, smiles1, smiles2)
            return Response(products, mimetype='chemical/x-daylight-smiles')
        else:
            return redirect(url_for('index'))
else:
    app = None


if __name__ == '__main__':
    if '--serve' in sys.argv:
        if app is None:
            print("The flask package is missing. Cannot run reaction server.", file=sys.stderr)
            print("Try `pip install flask`", file=sys.stderr)
            sys.exit(-1)
        else:
            sys.exit(app.run(host='0.0.0.0', port=8081, debug=DEBUG))
    else:
        rxn = sys.argv[1]
        smiles_files = sys.argv[2:]
        fs = []
        try:
            for smiles_file in smiles_files:
                fs.append(open(smiles_file))
            for line in smiles_reaction_matrix(rxn, *fs):
                print(line)
        finally:
            for f in fs:
                if f and not f.closed:
                    f.close()
            fs = []


