#!/usr/bin/env python2

"""
Script to annotate MS2LDA documents using MAGMa

Writes json to stdout, according to following schema:
 
doc_annotations = [{
    'name': 'document name',            # The MS file name
    'scan': scanid,                     # The scanid in the magma DB, this can be ignored!
    'mol': 'molblock',                  # The chemical structure of the parent (reference) molecule
    'features': [{                      # The list of features
        'name': 'feature name',         # Feature name
        'intensity': intensity,         # Feature intensity?
        'type': 'loss' or 'fragment',   # Feature type
        'matches': [{                   # List of MS/MS peaks matching the feature 
            'smiles': 'smiles',         # Smiles of the substructure suggested by MAGMa
            'mol': 'molblock',          # Chemical structure (of the MAGMa substructure).
                                        # Only present for fragments, if possible to generate from smiles
            'mz': mz,                   # m/z of the peak
            'fragatoms': [atom indices] # mapping of the substructure on the parent molecule
            }]
        }]
    }]
"""

import sys
import magma
from magma.models import Molecule, Fragment
from magma.pars import Hmass, elmass
from sqlalchemy import create_engine, between #collate
from sqlalchemy.orm import sessionmaker #, exc
from rdkit import Chem
from rdkit.Chem import AllChem
import json
sys.path.insert(0, './ms2ldaviz/lda/code')
from ms2lda_feature_extraction import *
import requests
from argparse import ArgumentParser, RawDescriptionHelpFormatter


# --- globals ---
version = '0.1.0'

server_address = "http://ms2lda.org/"

spectrum = {} # dict containing de spectrum name for a scanid
scan = {} # dict containing de scanid for a spectrum name
scans_for_molid = {} # dict containing lists of scanids for each molid
motifspectra = {}
doc_annotations = []
loaders = {'GNPS': LoadGNPS, 'MSP': LoadMSP, 'MGF': LoadMGF}


def loss2smiles(molblock, atomlist):
    """
    Create smiles of the loss(es)
    from molblock and list of fragment atoms
    """
    atoms = [int(a) for a in atomlist.split(',')]
    mol = Chem.MolFromMolBlock(molblock)
    emol = Chem.EditableMol(mol)
    for atom in reversed(range(mol.GetNumAtoms())):
        if atom in atoms:
            emol.RemoveAtom(atom)
    frag = emol.GetMol()
    return Chem.MolToSmiles(frag)


def read_spectra():
    rt = 1
    peaklists={}
    for i in ms1:
        peaklists[rt] = []
        meta = metadata[i.name]
        if args.molname_field:
            compound = meta[args.molname_field]
        else:
            compound = i.name
        molid = struct_engine.add_smiles(meta['smiles'], compound)
        if molid in scans_for_molid:
            scans_for_molid[molid].append(rt)
        else:
            scans_for_molid[molid] = [rt]
        magma_session.commit()
        mim = db_session.query(Molecule.mim).filter_by(molid=molid).one()[0]
        meta['parentmass'] = mim + args.mode * (Hmass-elmass)
        scan[i.name] = rt
        spectrum[rt] = i.name
        rt += 2
    if not db_exists:
        for i in ms2:
            peaklists[scan[i[3].name]].append([i[0], i[2]])
        for rt in peaklists:
            if len(peaklists[rt]) > 0:
                maxint = sorted(peaklists[rt], key = lambda p: p[1])[-1][1]
                peaklist = [p for p in peaklists[rt]]
                ms_data_engine.store_peak_list(
                    rt,
                    rt,
                    metadata[spectrum[rt]]['parentmass'],
                    [metadata[spectrum[rt]]['parentmass'], 10000],
                    peaklist)


def run_magma_annotate():
    for molid in scans_for_molid:
        annotate_engine.scans = []
        annotate_engine.build_spectra(scans_for_molid[molid])
        annotate_engine.search_structures(molids = [molid], fast=True)
        magma_session.commit()


def create_doc_annotations():
    url = server_address + 'basicviz/get_all_doc_data/{}'.format(args.experiment_id)
    response = requests.get(url)
    doc_features = response.json()
    # json.dump(doc_features, open('tmp.json', 'w'))
    # doc_features = list(json.load(open('tmp.json', 'r')))

    for doc, feature_list, motifs in doc_features:
        if doc not in scan:
            continue
        results = db_session.query(Fragment.mz, Molecule.mol).\
                filter(Fragment.scanid == scan[doc]).\
                join(Molecule, Fragment.molid==Molecule.molid).all()
        sys.stderr.write('doc: ' + doc + '; scan: ' + str(scan[doc]) + '\n')
        if len(results) != 1:
            sys.stderr.write('---> No molecule !\n')
            continue
        parent_mz, parent_mol = results[0]
        ddoc = {'name': doc, 'scan': scan[doc], 'mol': parent_mol, 'features': []}
        # read fragments and losses
        frags = db_session.query(Fragment).\
                filter(Fragment.scanid == scan[doc] + 1)
        for f, i in feature_list:
            feature_type, massbin = f.split('_')
            feature = {'name': f, 'intensity': i, 'type': feature_type, 'matches': []}
            massbin = float(massbin)
            if feature_type == 'fragment':
                matching_frags = frags.filter(between(Fragment.mz, massbin-args.binsize/2, massbin+args.binsize/2)).all()
            else:
                matching_frags = frags.filter(between(parent_mz - Fragment.mz, massbin-args.binsize/2, massbin+args.binsize/2)).all()
            for frag in matching_frags:
                if frag.atoms in [m['fragatoms'] for m in feature['matches']]:
                    continue # avoid redundant matches
                smiles = frag.smiles if feature_type == 'fragment' else loss2smiles(parent_mol, frag.atoms)
                match = {'smiles': smiles, 'mz': frag.mz, 'fragatoms': frag.atoms}
                feature['matches'].append(match)
            ddoc['features'].append(feature)
        doc_annotations.append(ddoc)

class stdout2stderr(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = sys.stderr
        return self
    def __exit__(self, *args):
        sys.stdout = self._stdout

def arg_parser():
    ap = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    ap.add_argument('-v', '--version', action='version', version='%(prog)s ' + str(version))
    ap.add_argument('-l', '--loader', help="Specify which loader to use (default: %(default)s)",
                    default='GNPS', choices=['GNPS', 'MSP', 'MGF'])
    ap.add_argument('-a', '--loader_arguments', help="Dictionary (as JSON) with argements for the MS loader (default: %(default)s)",
                    default='{}', type=str)
    ap.add_argument('-b', '--binsize', help="Bin size used in MS2LDA experiment (Da) (default: %(default)s)",
                    default=0.005,type=float)
    ap.add_argument('-m', '--mode', help="Ionisation mode (default: %(default)s)",
                    default=1, choices=[-1,1], type=float)
    ap.add_argument('-n', '--molname_field', help="Meta data field containing the compound name", type=str)
    ap.add_argument('-p', '--mz_precision', help="Maximum relative m/z error (ppm) (default: %(default)s)",
                    default=5,type=float)
    ap.add_argument('-q', '--mz_precision_abs', help="Maximum absolute m/z error (Da) (default: %(default)s)",
                    default=0.001,type=float)
    ap.add_argument('experiment_id', help="Experiment ID", type=int)
    ap.add_argument('magma_db', help="(Non-)existing magma database", type=str)
    ap.add_argument('spectra_path', help="path to spectra", type=str)
    return ap


# MAIN
import glob

args = arg_parser().parse_args(sys.argv[1:])

# if magma_db exists it is reused, skipping the parts that create and run the magma job
db_exists = os.path.isfile(args.magma_db)
magma_session =   magma.MagmaSession(args.magma_db, 'ms2lda_dataset', 'debug')
struct_engine =   magma_session.get_structure_engine(pubchem_names=True)
if not db_exists:
    ms_data_engine =  magma_session.get_ms_data_engine(
                      ionisation_mode=args.mode,
                      abs_peak_cutoff=0,
                      mz_precision=args.mz_precision,
                      mz_precision_abs=args.mz_precision_abs)
    annotate_engine = magma_session.get_annotate_engine(
                      ms_intensity_cutoff=0,
                      msms_intensity_cutoff=0)

# create sqlalchemy connection to magma database
engine = create_engine('sqlite:///' + args.magma_db)
session = sessionmaker()
session.configure(bind=engine)
db_session = session()

if os.path.isdir(args.spectra_path):
    files = glob.glob(args.spectra_path + '/*')
else:
    files = [args.spectra_path]
loader = loaders[args.loader](**json.loads(args.loader_arguments))
with stdout2stderr():
    ms1, ms2, metadata = loader.load_spectra(files)
read_spectra()

if not db_exists:
    run_magma_annotate()
create_doc_annotations()
json.dump(doc_annotations, sys.stdout, indent=2)


