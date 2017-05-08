from ase import Atoms, Atom

def pif_to_atoms(pif):
    props = pif['properties']
    pos = [p for p in props if p['name'] == 'Positions']
    if len(pos) < 1:
        raise ValueError('No Positions found. Cannot reconstruct atoms')
    elif len(pos) > 1:
        raise ValueError('Multiple positions found in PIF. Ambiguous atoms definitions')

    pos = pos[0]
    atoms = [c for c in pos['conditions'] if c['name'] == 'ASE atoms']

    if len(atoms) < 1:
        raise ValueError('No ASE atoms found in conditions. Cannot reconstruct atoms')
    elif len(atoms) > 1:
        raise ValueError('Multiple ASE atoms found in PIF. Ambiguous atoms definitions')
    else:
        doc = atoms[0]['scalars'][0]

    atoms = dict_to_atoms(doc)
    return atoms

def dict_to_atoms(doc):
    """
    Takes in a PIF dictionary and creates an atoms object. Lightly revised from
    from Kitchin group repo: https://github.com/jkitchin/vasp/blob/master/vasp/mongo.py
    """
    atoms = Atoms([Atom(atom['symbol'],
                            atom['position'],
                            tag=atom['tag'],
                            momentum=atom['momentum'],
                            magmom=atom['magmom'],
                            charge=atom['charge'])
                       for atom in doc['atoms']],
                      cell=doc['cell'],
                      pbc=doc['pbc'],
                      info=doc['info'],
                      constraint=[dict2constraint(c) for c in doc['constraints']])

    from ase.calculators.singlepoint import SinglePointCalculator
    calc = SinglePointCalculator(energy=doc.get('calculator.energy', None),
                                 forces=doc.get('calculator.forces', None),
                                 stress=doc.get('calculator.stress', None),
                                 atoms=atoms)
    atoms.set_calculator(calc)
    return atoms
