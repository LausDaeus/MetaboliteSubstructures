from collections import  Counter

import src.cs_algorithms.nx_structures as nx
import src.molecule_graphs.molecule as m

class PathStructure(m.Molecule):

    def __init__(self, label, path_vertices, parent_molecule):

        m.Molecule.__init__(self, label)
        self.label = label
        self.parent_struct_map = {} #map vertices from parent molecule to struct's
        self.struct_parent_map = {} #map vertices from struct to parents's


        for atom in path_vertices:
            # if not atom.aromatic:
            new_atom = self.add_atom(atom.label, atom.aromatic)
            self.parent_struct_map[atom] = new_atom
            self.struct_parent_map[new_atom] = atom

        for atom in path_vertices:
            for neighbour in path_vertices:
                # Test if there is an edge to any of the other atoms in the path
                edge = parent_molecule.get_edge_if_exists(atom, neighbour)
                if isinstance(edge, m.Bond):
                    self.add_bond(self.parent_struct_map[atom], self.parent_struct_map[neighbour], edge.type,edge.label)

        # Return the original SMILES string
    def __str__(self):
        return self.label

    def __repr__(self):
        return '<Path structure %s at %s>' % (str(self), id(self))

    def check_if_duplicate(self, molecule, vertices, path_structures, multi_dict, nx_structures):
        """
        Adds the other_struct graph to path_structures or multiple_structures depending on if there is a duplicate of it.

        :param other_struct: the graph object which is to be tested against the other graphs
        :param molecule: the molecule which contained the subgraph isomorphic to the structure
        :param vertices: the dictionary mapping the other_struct vertices to the molecule vertices
        :return: the graph object which has been newly added or updated
        """
        # A copy of path_structures is made so that it can be altered while also looping over the label

        for structure in path_structures:

            nx_mapping = nx.isomorphism(self, structure,nx_structures)
            if not nx_mapping:
                # Continues to next structure for testing as they are not isomorphic
                continue
            if nx_mapping:

                # A multiple structure has been made which matches a structure that has already been added to the CS

                # Maps the vertices from the other_struct to the target
                isomorphic_mapping = {}
                for target_position in nx_mapping:
                    target_match = structure.vertex_from_position(target_position)
                    other_struct_match = self.vertex_from_position(nx_mapping[target_position])
                    mole_vertex = vertices[other_struct_match]
                    isomorphic_mapping[target_match] = mole_vertex
                if molecule in path_structures[structure]:
                    # If the molecule vertices are different to the vertices currently in dictionary
                    # then there are multiple subgraphs in molecule isomorphic to other_struct
                    if Counter(isomorphic_mapping.values()) != Counter(
                            path_structures[structure][molecule].values()):
                        self.add_structure_to_multiple_dictionary(structure, molecule,
                                                                  isomorphic_mapping, path_structures, multi_dict)
                else:
                    path_structures[structure][molecule] = isomorphic_mapping
                unique_structure = structure
                break
        else:
            path_structures[self] = {molecule: vertices}
            unique_structure = self
        return unique_structure

    def add_structure_to_multiple_dictionary(self, structure, molecule, mapping, path_structures, multi_dict):
        """
        Adds a structure which has appeared more than once in a molecule to the multiple_structures dictionary

        :param structure: the structure which has appeared more than once in a molecule
        :param molecule: the molecule that has repeated instances of the structure
        :param mapping: a mapping from the structure to the molecule vertices which comes from the isomorphism check
        :return: None
        """
        if structure in multi_dict and molecule in multi_dict[structure]:
            # Creates a list of lists which sets out the vertices that are mapped to each instance of the structure
            # within the molecule, this ensures that an instance isn't duplicated within the dictionary.
            multiples = [[vertex[index] for vertex in multi_dict[structure][molecule].values()]
                         for index in range(len(multi_dict[structure][molecule].values()[0]))]
            for v in multiples:
                if Counter(mapping.values()) == Counter(v):
                    return
            # Store the molecule vertex which maps to each of the structure vertices
            for key in mapping:
                multi_dict[structure][molecule][key].append(mapping[key])
        else:
            if structure in multi_dict:
                multi_dict[structure][molecule] = {}
            elif structure not in multi_dict:
                multi_dict[structure] = {molecule: {}}
            # When the structure is encountered for the second time in the molecule the mapping from path_structures
            # is also stored into the multiple_structures dictionary so it includes all instances of the structure
            for key in mapping:
                multi_dict[structure][molecule][key] = [path_structures[structure][molecule][key],
                                                                         mapping[key]]

