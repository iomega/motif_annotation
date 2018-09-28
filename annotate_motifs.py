"""
Script to annotate MS2LDA motifs using MAGMa

Writes json to stdout, according to following schema:
 
motifannotations = [{
    'name': 'motif name',
    'comments': 'comments'
    'features': [{
        'name': 'feature name',
        'type': 'loss' or 'fragment'
        'substr' [{
            'smiles': 'smiles',
            'mol': 'molblock',
            'matches': [{
                'spectrum': 'spectrum name',
                'mz': mz,
                'scan': scanid,
                'mol': 'molblock',
                'fragatoms': [atom indices]
                }]
            }]
        }]
    }]
"""

import sys, os
import magma
from sqlalchemy import create_engine, desc, func, between #collate
from sqlalchemy.orm import sessionmaker, aliased #, exc
from magma.models import Base, Molecule, Reaction, fill_molecules_reactions, Scan, Peak, Fragment, Run
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import groupby
from chemspipy import ChemSpider
import json
sys.path.insert(0, './ms2ldaviz/lda/code')
from ms2lda_feature_extraction import *
import requests
from argparse import ArgumentParser, RawDescriptionHelpFormatter


# --- globals ---
version = '0.1.0'

cs = ChemSpider('b07b7eb2-0ba7-40db-abc3-2a77a7544a3d')
server_address = "http://ms2lda.org/"

spectrum = {} # dict containing de spectrum name for a scanid
scan = {} # dict containing de scanid for a spectrum name
scans_for_molid = {} # dict containing lists of scanids for each molid
motifspectra = {}
motifannotations = []


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
        compound = meta['compound'].replace('"','')
        if not db_exists:
            if args.chemspider:
                csresults = cs.search(meta['InChIKey'])
                if csresults:
                    molid = struct_engine.add_structure(csresults[0].mol_2d, compound, 0.0, 0)
                else:
                    sys.stderr.write('--> No Chemspider for ' + compound + '\n')
                    molid = struct_engine.add_smiles(meta['smiles'], compound)
            else:
                molid = struct_engine.add_smiles(meta['smiles'], compound)
            if molid in scans_for_molid:
                scans_for_molid[molid].append(rt)
            else:
                scans_for_molid[molid] = [rt]
        scan[i.name] = rt
        spectrum[rt] = i.name
        rt += 2
    if not db_exists:
        for i in ms2:
            peaklists[scan[i[3].name]].append([i[0], i[2]])
        for rt in peaklists:
            if len(peaklists[rt]) > 0:
                maxint = sorted(peaklists[rt], key = lambda p: p[1])[-1][1]
                peaklist = [p for p in peaklists[rt] if p[1] > args.minrelint * maxint]
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


# read motif in spectra mapping
def read_motif_mapping():
    """
    motifspectra = {motifname: [spectra]}
    """
    url = server_address + 'basicviz/get_doc_m2m/{}'.format(args.experiment_id)
    response = requests.get(url)
    for m, f, p, o in response.json():
        if m in motifspectra:
            motifspectra[m].append(f)
        else:
            motifspectra[m] = [f]


def create_motif_annotations():
    # create sqlalchemy connection to magma database
    engine = create_engine('sqlite:///' + args.magma_db)
    session = sessionmaker()
    session.configure(bind=engine)
    db_session = session()

    # parse MS2LDA motifs
    for motif in sorted(motifspectra, key=lambda k: int(k[6:])):
        try:
            motiffile = open(args.motif_path + motif + '.m2m')
        except:
            continue
        ma = {'name': motif, 'comments': '', 'features': []}
        # read fragments and losses
        for line in motiffile:
            line = line[:-2]
            if line[0] == '#':
                ma['comments'] += '<b>'+ line + '</b><br>'
            if line.startswith('fragment_'):
                feature = {'name': line, 'type': 'fragment', 'substr': []}
                massbin, weight = line[9:-1].split(',')
                massbin = float(massbin)
                matching_frags = db_session.query(Fragment, Molecule.mol).\
                        filter(Fragment.scanid.in_([scan[s]+1 for s in motifspectra[motif] if s in scan])).\
                        filter(between(Fragment.mz, massbin-0.0025, massbin+0.0025)).\
                        join(Molecule, Fragment.molid==Molecule.molid).all()
            if line.startswith('loss_'):
                feature = {'name': line, 'type': 'loss', 'substr': []}
                massbin, weight = line[5:-1].split(',')
                massbin = float(massbin)
                frag_alias = aliased(Fragment)
                matching_frags = db_session.query(Fragment, Molecule.mol).\
                        filter(Fragment.scanid.in_([scan[s]+1 for s in motifspectra[motif] if s in scan])).\
                        join(frag_alias, Fragment.parentfragid==frag_alias.fragid).\
                        filter(between(frag_alias.mz - Fragment.mz, massbin-0.0025, massbin+0.0025)).\
                        join(Molecule, Fragment.molid==Molecule.molid).all()
            if line.startswith('fragment_') or line.startswith('loss_'):
                frags = {}
                for frag, molblock in matching_frags:
                    smiles = frag.smiles if feature['type'] == 'fragment' else loss2smiles(molblock, frag.atoms)
                    if smiles in frags:
                        frags[smiles].append([frag, molblock])
                    else:
                        frags[smiles] = [[frag, molblock]]
                for smiles, matches in sorted(frags.iteritems(), key = lambda (k,v): len(v), reverse = True):
                    substr = {'smiles': smiles, 'matches': []}
                    if feature['type'] == 'fragment':
                        try:
                            mol = Chem.MolFromSmiles(substr['smiles'])
                            AllChem.Compute2DCoords(mol)
                            substr['mol'] = Chem.MolToMolBlock(mol)
                        except:
                            pass
                    for frag, molblock in matches:
                        substr['matches'].append({
                            'spectrum': spectrum[frag.scanid-1],
                            'mz': frag.mz,
                            'scan': frag.scanid - 1,
                            'mol': molblock,
                            'fragatoms': frag.atoms
                            })
                    feature['substr'].append(substr)
                ma['features'].append(feature)
        motifannotations.append(ma)

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
    ap.add_argument('-i', '--minrelint', help="Threshold on intensity relative to basepeak (default: %(default)s)", default=0.05, type=float)
    ap.add_argument('-c', '--chemspider', help="Get structures from ChemSpider (default: %(default)s)", default=False, action='store_true')
    ap.add_argument('experiment_id', help="Experiment ID", type=int)
    ap.add_argument('magma_db', help="(Non-)existing magma database", type=str)
    ap.add_argument('spectra_path', help="path to spectra", type=str)
    ap.add_argument('motif_path', help="path to motifs", type=str)
    return ap


# MAIN
import glob

args = arg_parser().parse_args(sys.argv[1:])

# if magma_db exists it is reused, skipping the parts that create and run the magma job
db_exists = os.path.isfile(args.magma_db)
if not db_exists:
    magma_session =   magma.MagmaSession(args.magma_db, 'ms2lda_dataset', 'debug')
    ms_data_engine =  magma_session.get_ms_data_engine(
                      abs_peak_cutoff=0,
                      mz_precision=20,
                      mz_precision_abs=0.005)
    struct_engine =   magma_session.get_structure_engine(pubchem_names=False)
    annotate_engine = magma_session.get_annotate_engine(
                      ms_intensity_cutoff=0,
                      msms_intensity_cutoff=0)

if os.path.isdir(args.spectra_path):
    files = glob.glob(args.spectra_path + '/*')
else:
    files = [args.spectra_path]
loader = LoadGNPS()
with stdout2stderr():
    ms1, ms2, metadata = loader.load_spectra(files)
read_spectra()

if not db_exists:
    run_magma_annotate()
read_motif_mapping()
create_motif_annotations()
json.dump(motifannotations, sys.stdout, indent=2)


