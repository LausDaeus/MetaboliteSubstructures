import copy

available_rings = None
organic_set = ['B','C','N','O','S','P','F','Cl','Br','I','*']
organic_aromatic = ['b','c','n','o','s','p']
other_atoms = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg',
               'Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
               'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
               'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
               'Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','Hf','Ta','W',
               'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
               'Fr','Ra','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Fl',
               'Lv','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
               'Er','Tm','Yb','Lu','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk',
               'Cf','Es','Fm','Md','No','Lr']
wild_card = '*'

def convert_to_smiles(original):
    molecule = copy.deepcopy(original)
    vertices = molecule.vertices()
    smiles = ""

    if not vertices:
        return smiles

    starting_vertices = molecule.create_spanning_tree()

    global available_rings
    available_rings = [True] * len(vertices)

    for root in starting_vertices:
        smiles +=  parse_atom(root, molecule)+"."

    smiles = smiles[:-1] # Remove extra "." separator

    molecule.clean_spanning_tree()
    clean_smiles_data(molecule)
    return smiles


def parse_atom(atom, molecule):
    #Create atom string

    #If the atom is in the organic subset and has no special features.
    if not atom.isotope and not atom.hydrogen and not atom.charge and \
            (atom.label in organic_set or atom.label in organic_aromatic):
        if atom.aromatic:
            atom_string = atom.label.lower()
        else:
            atom_string = atom.label

    #Else use square bracket notation and properties as defined by SMILES.
    else:
        atom_string = "["

        if atom.isotope:
            atom_string += str(atom.isotope)

        if atom.hydrogen:
            if atom.hydrogen == 1:
                atom_string += "H"
            else:
                atom_string += "H"+str(atom.hydrogen)

        atom_string += atom.label

        if atom.charge:
            atom_string += get_charge_string(atom.charge)

        atom_string += "]"

    #Sets the atom as parsed
    atom.parsed = True

    #Counter of number of branches.
    chain_count = 0

    #Scans for rings and counts branches.
    for neighbour in molecule.neighbours(atom):
        bond =  molecule.adjacency_dictionary[atom][neighbour]
        bondOpp = molecule.adjacency_dictionary[neighbour][atom]
        #Checks if bond is not included in spanning tree.
        #If it is not, add a ring label.
        if not bond.visited:
            atom_string += choose_bond_string(bond)
            #If there already is a label for this ring.
            if hasattr(bond, "ringID"):
                atom_string += get_ring_ID_syntax(bond.ringID)
                free_ring_ID(bond.ringID)
                del bond.ringID
            #Else, create a label for this ring.
            else:
                ringID = get_next_ring_ID()
                atom_string+= get_ring_ID_syntax(ringID)
                bond.ringID = ringID
        #All non-ring bonds are branches and are counted.
        elif not (hasattr(neighbour, "parsed") and neighbour.parsed):
            chain_count += 1

    #Add branches
    for neighbour in molecule.neighbours(atom):
        #Makes sure atom is only parsed once
        if not (hasattr(neighbour, "parsed") and neighbour.parsed):
            bond =  molecule.adjacency_dictionary[atom][neighbour]
            #If the bond is a branch.
            if bond.visited:
                #print "Next atom:", neighbour
                #If it is the only/last branch, don't use brackets.
                if chain_count == 1:
                    atom_string += choose_bond_string(bond)
                    #Recursive step for next atom.
                    atom_string += parse_atom(neighbour, molecule)
                #Otherwise, use brackets.
                else:
                    chain_count -= 1
                    atom_string += "("+choose_bond_string(bond)
                    #Recursive step for next atom.
                    atom_string += parse_atom(neighbour, molecule)+")"

    return atom_string

def choose_bond_string(bond):
    if bond.label == "-":
        return ""
    else:
        return bond.label


def get_next_ring_ID():
    global  available_rings
    try:
        next_val = available_rings.index(True)
        available_rings[next_val] = False
        return next_val+1
    except ValueError:
        available_rings = available_rings + len(available_rings) * True
        return get_next_ring_ID()

def get_ring_ID_syntax(ringID):
    if ringID > 9:
        return "%"+str(ringID)
    else:
        return str(ringID)

def free_ring_ID(value):
    available_rings[value-1] = True

def get_charge_string(charge):
    if charge == -1:
        return '-'
    if charge == 1:
        return '+'
    return str(charge)

def clean_smiles_data(molecule):
    for atom in molecule.vertices():
        if hasattr(atom, "parsed"):
            del atom.parsed
    for bond in molecule.edges():
        if hasattr(bond, "ringID"):
            del bond.ringID
