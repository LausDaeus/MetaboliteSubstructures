from enum import Enum

import src.molecule_graphs.graph as graph
from src.parsing import molecule_smiles_parser


class BondType(Enum):
    __order__ = 'single double triple quadruple aromatic'
    single = 'single'
    double = 'double'
    triple = 'triple'
    quadruple = 'quadruple'
    aromatic = 'aromatic'


# Represents an atom within the molecule, includes the atoms chemical properties
class Atom(graph.Vertex):
    """
    A single atom of a molecule which includes the chemical information

    :param label: a string representing the atomic element
    :param isotope: the isotopic number if specified
    :param hydrogen: the number of explicit hydrogens
    :param charge: the charge of the atom
    :param aromatic: sets if the atom is part of an aromatic ring
    """
    def __init__(self, label, isotope=None, hydrogen=None, charge=None, aromatic=False):
        graph.Vertex.__init__(self, label)
        self.isotope = isotope
        self.hydrogen = hydrogen
        self.charge = charge
        self.aromatic = aromatic
        self.ring_break = False

    def __str__(self):
        """
        String representation based on atom type
        :return: string representation
        """
        if self.aromatic:
            return 'aromatic %s atom' % self.label
        else:
            return '%s atom' % self.label

    def __repr__(self):
        """
        Unique representation based on atom type
        :return: None
        """

        # if self.aromatic:
        #     return '<Aromatic Atom %s at %s>' % (self.label, id(self))
        # else:
        #     return '<Atom %s at %s>' % (self.label, id(self))
        return '<Atom '+self.label+' at '+str(self.ID)+'>'

    def clone(self):
        new_atom = Atom(self.label, self.isotope, self.hydrogen, self.charge, self.aromatic)
        return new_atom

class Bond(graph.Edge):
    """
    A single bond of the molecule which has the bond type

    :param origin: one of the endpoints
    :param destination: one of the endpoints
    """
    def __init__(self, origin, destination, type, label):
        graph.Edge.__init__(self, origin, destination, label)
        self.type = type

    def __str__(self):
        """
        String representation that is different depending on the type of bond
        :return: string representation
        """
        return self.type.name + " bond"

    def __repr__(self):
        """
        Unique representation that is different depending on the type of bond
        :return: unique representation of the bond
        """
        return '<'+self.type.name+' bond at %s>' % id(self)

    def clone(self, origin=None, destination=None):
       return Bond(origin, destination, self.type, self.label)

class Molecule(graph.Graph):
    """
    A Molecule which inherits from the Graph class and contains atoms and bonds

    -smiles_string: the smiles string which is associated with the molecule
    """
    ID_counter = 0

    def __init__(self, smiles):
        graph.Graph.__init__(self)
        self._smiles_string = smiles
        self.ID = Molecule.ID_counter
        Molecule.ID_counter += 1


    # Return the original SMILES string
    def __str__(self):
        return self._smiles_string

    def __repr__(self):
        return '<Molecule %s at %s>' % (str(self), id(self))

    # def copy_self(self):
    #     duplicate = Molecule(str(self.characteristic_substructure))
    #     for vertex in self.vertices():
    #         possible_location.vertex_to_graph(vertex)
    #         for neighbour in self.characteristic_substructure.adjacency_dictionary[vertex]:
    #             possible_location.adjacency_dictionary[vertex][neighbour] = copy(self.characteristic_substructure.
    #                                                                              adjacency_dictionary
    #                                                                              [vertex][neighbour])

    def add_atom(self, label, isotope=None, hydrogen=None, charge=None, aromatic=False):
        """
        Adds a new atom to the molecule
        :param label: the atomic element
        :param isotope: the isotope number
        :param hydrogen: the number of explicit hydrogens
        :param charge: the charge of the atom
        :return: the new atom object that has been added
        """
        new_atom = Atom(label, isotope, hydrogen, charge, aromatic)
        graph.Graph.vertex_to_graph(self, new_atom)
        return new_atom

    def add_bond(self, origin, destination, type, label):
        new_bond = Bond(origin, destination, type, label)
        graph.Graph.edge_to_graph(self, origin, destination, new_bond)
        return new_bond

    def set_list_membership(self, list_number):
        """
        When given lists of molecules, sets this molecule's list ID.
        :param list_number: listID
        :return: nothing
        """
        self.list_member = list_number

    def molecule_info(self, reverse=False):
        info = ""
        for idx,atom in enumerate(self.vertices()):
            atom.identifier = idx

        for atom in self.vertices():
            info += "\n["+str(atom.identifier)+"]"+atom.label+": "
            for neighbour in self.neighbours(atom):
                if reverse:
                    info += "("+str(self.get_edge_if_exists(neighbour, atom)) + "[" + str(neighbour.identifier) + "]) "
                else:
                    info += "("+str(self.get_edge_if_exists(atom, neighbour)) + "[" + str(neighbour.identifier) + "]) "

        for atom in self.vertices():
            del atom.identifier

        return info

    def to_smiles(self):
        return molecule_smiles_parser.convert_to_smiles(self)