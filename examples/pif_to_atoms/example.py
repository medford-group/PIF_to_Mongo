import sys
import json
sys.path.insert(0,'../..')

from ase_utils import pif_to_atoms

pif = json.load(open('TiO2_N2.pif'))

atoms = pif_to_atoms(pif)

print atoms.get_potential_energy()
