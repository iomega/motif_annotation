"""
Script to annotate MS2LDA motifs using MAGMa

Usage: python annotate_motifs.py mspfile.mps magmadb.db field_containing_spectrumid ms2lda.csv path_to_motifs_including_prefix > annotation.json


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


# --- globals ---
mspfile = open(sys.argv[1], 'r')
magmadb = sys.argv[2]
spectrumid = sys.argv[3]
mappingfile = sys.argv[4]
motifpath = sys.argv[5]

# if magmadb exists it is reused, skipping the parts that create and run the magma job 
db_exists = os.path.isfile(magmadb)
if not db_exists:
    magma_session =   magma.MagmaSession(magmadb, 'ms2lda_dataset', 'debug')
    ms_data_engine =  magma_session.get_ms_data_engine(
                      abs_peak_cutoff=0,
                      mz_precision=20,
                      mz_precision_abs=0.005)
    struct_engine =   magma_session.get_structure_engine(pubchem_names=False)
    annotate_engine = magma_session.get_annotate_engine(
                      ms_intensity_cutoff=0,
                      msms_intensity_cutoff=0)

cs = ChemSpider('b07b7eb2-0ba7-40db-abc3-2a77a7544a3d')

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


def parse_mspfile():
    "Parse mspfile store spectra in a MAGMa database"

    rt = 1
    while True:
        l = mspfile.readline()
        if l == '':
            break
        l = l.rstrip('\n\r')
        if l[:5] == "NAME:":
            name = l[6:]
        if l[:12] == 'PRECURSORMZ:':
            pmz = float(l[13:])
        if l[:7] == 'SMILES:':
            smiles = l[8:]
        if l[:9] == 'INCHIKEY:' and not db_exists:
            results = cs.search(l[10:])
            if results:
                molid = struct_engine.add_structure(results[0].mol_2d, name, 0.0, 0)
            else:
                sys.stderr.write('--> No Chemspider for ' + name + '\n')
                molid = struct_engine.add_smiles(smiles, name)
        if l[:len(spectrumid)] == spectrumid:
            scan[l[len(spectrumid)+1:]] = rt
            spectrum[rt] = l[len(spectrumid)+1:]
        if l[:10] == 'Num Peaks:' and pmz < 1000:
            if not db_exists:
                peaklist = []
                for i in range(int(l[11:])):
                    mz, intensity = mspfile.readline().rstrip('\n\r').split()
                    peak = [float(mz), float(intensity)]
                    if peak not in peaklist:
                        peaklist.append(peak)
                ms_data_engine.store_peak_list(
                    rt,
                    rt,
                    pmz,
                    [pmz, 10000],
                    peaklist
                    )
                if molid in scans_for_molid:
                    scans_for_molid[molid].append(rt)
                else:
                    scans_for_molid[molid] = [rt]
            rt += 2


def run_magma_annotate():
    for molid in scans_for_molid:
        annotate_engine.scans = []
        annotate_engine.build_spectra(scans_for_molid[molid])
        annotate_engine.search_structures(molids = [molid], fast=True, ncpus=3)
        magma_session.commit()


# read motif in spectra mapping
def read_motif_mapping():
    """
    motifspectra = {motifname: [spectra]}
    """
    for line in open(mappingfile, 'r'):
        if line[:10] == '"Document"':
            continue
        f,m,p,null,null,null,name = [a.replace('"','') for a in line[:-1].split('","')]
        if m in motifspectra:
            motifspectra[m].append(f[:-3])
        else:
            motifspectra[m] = [f[:-3]]


def create_motif_annotations():
    # create sqlalchemy connection to magma database
    engine = create_engine('sqlite:///' + magmadb)
    session = sessionmaker()
    session.configure(bind=engine)
    db_session = session()

    # parse MS2LDA motifs
    for motif in sorted(motifspectra, key=lambda k: int(k[6:])):
        try:
            motiffile = open(motifpath + motif + '.m2m')
        except:
            continue
        ma = {'name': motif, 'comments': '', 'features': []}
        # read fragments and losses
        for line in motiffile:
            line = line[:-2]
            if line[0] == '#':
                ma['comments'] += '<b>'+ line + '</b><br>'
            if line[:9] == 'fragment_':
                feature = {'name': line, 'type': 'fragment', 'substr': []}
                massbin, weight = line[9:-1].split(',')
                massbin = float(massbin)
                matching_frags = db_session.query(Fragment, Molecule.mol).\
                        filter(Fragment.scanid.in_([scan[s]+1 for s in motifspectra[motif] if s in scan])).\
                        filter(between(Fragment.mz, massbin-0.0025, massbin+0.0025)).\
                        join(Molecule, Fragment.molid==Molecule.molid).all()
            if line[:5] == 'loss_':
                feature = {'name': line, 'type': 'loss', 'substr': []}
                massbin, weight = line[5:-1].split(',')
                massbin = float(massbin)
                frag_alias = aliased(Fragment)
                mol_alias = aliased(Molecule)
                matching_frags = db_session.query(Fragment, mol_alias.mol).\
                        filter(Fragment.scanid.in_([scan[s]+1 for s in motifspectra[motif] if s in scan])).\
                        join(frag_alias, Fragment.parentfragid==frag_alias.fragid).\
                        filter(between(frag_alias.mz - Fragment.mz, massbin-0.0025, massbin+0.0025)).\
                        join(mol_alias, Fragment.molid==mol_alias.molid).all()
            if line[:9] == 'fragment_' or line[:5] == 'loss_':
                frags = {}
                for frag, molblock in matching_frags:
                    smiles = frag.smiles if feature['type'] == 'fragment' else loss2smiles(molblock, frag.atoms)
                    if smiles in frags:
                        frags[smiles].append([frag, molblock.replace('\n','\\n')])
                    else:
                        frags[smiles] = [[frag, molblock.replace('\n','\\n')]]
                for smiles, matches in sorted(frags.iteritems(), key = lambda (k,v): len(v), reverse = True):
                    substr = {'smiles': smiles, 'matches': []}
                    if feature['type'] == 'fragment':
                        try:
                            mol = Chem.MolFromSmiles(substr['smiles'])
                            AllChem.Compute2DCoords(mol)
                            substr['mol'] = Chem.MolToMolBlock(mol)
                        except:
                            pass
                    for match in matches:
                        substr['matches'].append({
                            'spectrum': spectrum[frag.scanid-1],
                            'mz': frag.mz,
                            'scan': frag.scanid - 1,
                            'mol': molblock, #.replace('\n','\\n'),
                            'fragatoms': frag.atoms
                            })
                    feature['substr'].append(substr)
                ma['features'].append(feature)
        motifannotations.append(ma)


# MAIN
parse_mspfile()
if not db_exists:
    run_magma_annotate()
read_motif_mapping()
create_motif_annotations()
json.dump(motifannotations, sys.stdout, indent=2)


