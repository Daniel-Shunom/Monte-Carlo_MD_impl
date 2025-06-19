import sys
from rdkit import Chem

def get_logP(peptide):
    logP = None
    norm_min = len(peptide.split('-'))*-2.09
    norm_max = len(peptide.split('-'))*3.64

    for pep in peptide.split('-'):
        if logP is None:
            logP = peptide_logp[pep]
        else:
            logP += peptide_logp[pep]

    norm_logP = (logP - norm_min)/(norm_max - norm_min)
    return logP, norm_logP


def get_norm_AP(AP, norm_min=0.97, norm_max=2.7):
    norm_AP = (AP - norm_min)/(norm_max - norm_min)
    return norm_AP


def get_short_sequence(peptide):
    short_pep = ''
    for p in peptide:
        short_pep = short_pep + peptide_shortSeq[p]
    return short_pep


def get_long_sequence(pep):
    peptide = []
    for p in pep:
        peptide.append(peptide_longSeq[p])
    return '-'.join(peptide)


def get_smiles(pep):
    m = Chem.rdmolfiles.MolFromSequence(pep)
    return Chem.MolToSmiles(m)    

'''
def get_single_smile_pg_fp(smile):

    try:
        m = Chem.MolFromSmiles(smile)

        if m is not None:
            afp = fp.v2_atom(smiles_string=smile, debug=1, polymer_fp_type='aT', ismolecule=1)
            mfp = fp.v2_moleculardescriptor(smiles_string=smile, debug=1, ismolecule=1)
            efp = fp.v2_extendeddescriptor(smile, debug=1, ismolecule=1)

            smile_fp = {'smiles': smile, **afp,**mfp,**efp}
            #print('Success, Fingerprinting ', smile)

    except:
        print('Failed Fingerprinting for ', smile)
        smile_fp = None

    return smile_fp
'''

peptide_logp = {'ILE': -1.12,
                'LEU': -1.25,
                'PHE': -1.71,
                'VAL': -0.46,
                'MET': -0.67,
                'PRO': 0.14,
                'TRP': -2.09,
                'THR': 0.25,
                'GLN': 0.77,
                'CYS': -0.02,
                'TYR': -0.71,
                'ALA': 0.5,
                'SER': 0.46,
                'ASN': 0.85,
                'ARG': 1.81,
                'GLY': 1.15,
                'LYS': 2.8,
                'GLU': 3.63,
                'HIS+': 2.33,
                'ASP': 3.64,
                'GLU0': 0.11,
                'HIS': 0.11}


peptide_shortSeq = {'ARG': 'R',
                'HIS': 'H',
                'LYS': 'K',
                'ASP': 'D',
                'GLU': 'E',
                'SER': 'S',
                'THR': 'T',
                'ASN': 'N',
                'GLN': 'Q',
                'CYS': 'C',
                'GLY': 'G',
                'PRO': 'P',
                'ALA': 'A',
                'ILE': 'I',
                'LEU': 'L',
                'MET': 'M',
                'PHE': 'F',
                'TRP': 'W',
                'TYR': 'Y',
                'VAL': 'V'}

peptide_longSeq = {'R': 'ARG',
                 'H': 'HIS',
                 'K': 'LYS',
                 'D': 'ASP',
                 'E': 'GLU',
                 'S': 'SER',
                 'T': 'THR',
                 'N': 'ASN',
                 'Q': 'GLN',
                 'C': 'CYS',
                 'G': 'GLY',
                 'P': 'PRO',
                 'A': 'ALA',
                 'I': 'ILE',
                 'L': 'LEU',
                 'M': 'MET',
                 'F': 'PHE',
                 'W': 'TRP',
                 'Y': 'TYR',
                 'V': 'VAL'}    
