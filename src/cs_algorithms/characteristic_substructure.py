# coding=utf-8
from collections import OrderedDict
import src.utils.verbose as verbose
import src.cs_algorithms.nx_structures as nx
from src.parsing.smiles_molecule_parser import Parser
import src.cs_algorithms.path_structure as ps

import src.molecule_graphs.molecule as m



# Implementation of Finding Characteristic Substructures for Metabolite Classes
# Ludwig, Hufsky, Elshamy, Böcker


class CSAlgorithm(object):
    def __init__(self, **kwargs):
        # The initial parameters for the algorithm
        self.length_start = None
        self.length_end = 5
        self.step = 4
        self.threshold = 0.8
        self.isomorphism_factor = 0.8
        self.list_threshold = 0.8
        self.list_number = -1
        self.using_lists = False
        self.not_random = False
        self.name = "CS"

        for key, value in kwargs.iteritems():
            if hasattr(self, key):
                if value:
                    setattr(self, key, value)
                    verbose.info("show_setting", key, value)


        # The structure which has been created through the combination of representative path structures
        self.characteristic_substructure = m.Molecule('Characteristic Substructure')
        # Indicate if the characteristic substructure contains a path structure yet
        self.cs_begun = False
        # List holding structures which have been added to the characteristic substructure
        self.cs_structures = []
        # Dictionary holding the locations of the molecules which map to the characteristic substructure
        # {molecule: {molecule vertex: cs vertex}}
        self.cs_locations = {}
        # All the given molecules
        self.molecules = []
        # Dictionary of path structures with the isomorphic molecules and mapping vertices from structure to molecule
        # {path structure: {molecule: {structure vertex: molecule vertex}}}
        self.path_structures = {}
        # Dictionary of structures that appear multiple times in molecules and list of molecule vertices that map them
        # {single structure: {molecule: {structure vertex: [molecule vertices]}}}
        self.multiple_dict = {}

        self.nx_structures = {}

    def clear_old_structures(self):
        self.multiple_dict.clear()
        self.path_structures.clear()
        self.nx_structures.clear()


    def clear_all_structures(self):
        self.cs_structures[:] = []
        self.cs_locations.clear()
        self.clear_old_structures()
        self.cs_begun = False

    def add_molecules(self, smiles_set, list_membership=None):
        if list_membership and self.using_lists:
            for smiles, membership in zip(smiles_set, list_membership):
                molecule = Parser().parse_smiles(smiles)
                molecule.set_list_membership(membership)
                self.molecules.append(molecule)
        else:
            for smiles in smiles_set:
                molecule = Parser().parse_smiles(smiles)
                self.molecules.append(molecule)

    def find_characteristic_substructure(self, paths):
        """
        Find the characteristic substructure for a set of molecules

        Calls the representative paths method for each of the lengths between the start length and the end
        Creates the CS with the most frequent path structure
        Each of the subsequent path structures is added to the CS in order of frequency
        Swaps the vertices in the molecules which map to the CS with the CS vertices

        :return: a molecule object that is the characteristic substructure of the list of molecules
        """
        self.clear_all_structures()

        if not paths:
            return None

        self.length_start = self.get_starting_length(paths)
        length = self.length_start
        verbose.info("starting_length",length)
        while length >= self.length_end:
            verbose.info("current_length",length, self.name)
            representative_paths = self._find_representative_paths(paths, length)
            sorted_dictionary = self._find_representative_structures(representative_paths)
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            #print [struct.label for struct in sorted_dictionary.keys()]
            for structure in sorted_dictionary.keys():
                self._add_structure_to_characteristic(structure)

            self.clear_old_structures()

            # The step size only increases if the characteristic substructure has been started
            if self.cs_begun:
                length -= self.step
            else:
                length -= 1

        return self.characteristic_substructure

    def find_all_representative_structures(self, paths, name="CS"):
        """
        Creates a list of all the structure of different lengths which are representative sorted in terms of frequency

        :return: list of structures which appear frequently in molecules
        """
        #print paths

        self.clear_all_structures()
        all_structures = {}

        if not paths:
            return None

        length = self.get_starting_length(paths)
        verbose.info("starting_struct_length", length)
        while length >= self.length_end:
            verbose.info("current_length",length, name)
            representative_paths = self._find_representative_paths(paths, length)
            #print representative_paths
            sorted_dictionary = self._find_representative_structures(representative_paths)
            #print sorted_dictionary
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_dictionary:
                all_structures.update(sorted_dictionary)

            self.clear_old_structures()
            # To get the structures of all lengths the step does not alter
            length -= 1

        representative_structures = OrderedDict(sorted(all_structures.items(), key=lambda x: x[1], reverse=True))
        return representative_structures

    def find_graphs_paths(self):
        """
        For each SMILES string a molecule object is created and all of its paths found

        :param smiles_set: a list of SMILES strings
        :return: a dictionary containing all of the paths as strings with their length as value
        """
        paths = []
        for molecule in self.molecules:
            molecule.find_all_paths(paths)

        verbose.info("show_path_info",paths)

        return paths

    def _find_representative_paths(self, paths, length):
        """
        Find the paths which occur with a high enough frequency
        
        :param length: an integer which sets the length of the paths that should be considered
        :return: a list of the paths as strings which are representative 
        """
        verbose.info("rep_path_start")
        representative_paths = []
        #print paths
        if(self.using_lists):
            for path in paths[length]:
            # Check all the paths which are of the chosen length
                counter = 0
                list_membership = [False] * self.list_number
                for molecule in self.molecules:
                    # Search the tuples for each path which consists of path string and vertices present in path
                    if path in molecule.paths:
                        counter += 1
                        list_membership[molecule.list_member] = True

                if float(counter) / len(self.molecules) >= self.threshold and \
                        float(list_membership.count(True)) / self.list_number >= self.list_threshold:

                    representative_paths.append(path)

        else:
            for path in paths[length]:
                # Check all the paths which are of the chosen length
                counter = 0
                for molecule in self.molecules:
                    # Search the tuples for each path which consists of path string and vertices present in path
                    if path in molecule.paths:
                        counter += 1

                if float(counter) / len(self.molecules) >= self.threshold:
                    representative_paths.append(path)


        verbose.info("rep_path_finish", len(representative_paths))
        return representative_paths

    def _find_representative_structures(self, rep_paths):
        """
        Find the path structures which appear with a high enough frequency.

        Ensures there are no duplicates in the dictionary of path structures so that the frequency count is accurate
        :param rep_paths: a list of the paths which are representative
        :return: a list of path structures which are representative sorted into order of frequency, highest to lowest
        """
        #start = datetime.now()
        #print rep_paths

        rep_structures = {}
        verbose.info("rep_struct_start")
        for path in rep_paths:
            for molecule in self.molecules:
                if path in molecule.paths:
                    #print len(molecule.paths[path])
                    path_vertices_groups = molecule.paths[path]  # Lists of path vertices
                    for vertices in path_vertices_groups:
                        # Test if the structure has already been encountered
                        path_structure = ps.PathStructure(path, vertices, molecule)
                        #self._check_structure_duplicates(path_structure, molecule, path_structure.struct_parent_map)
                        path_structure.check_if_duplicate(molecule, path_structure.struct_parent_map,
                                                          self.path_structures, self.multiple_dict, self.nx_structures)

        for structure in self.path_structures:
            relative_frequency = len(self.path_structures[structure].keys()) / float(len(self.molecules))
            if float(relative_frequency) >= self.threshold and self.check_list_frequency(self.path_structures[structure].keys()):
                # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
                rep_structures[structure] = relative_frequency
        # Sort dictionary based on frequency highest to lowest and return the path structures only

        if self.not_random:
            temp = OrderedDict(sorted(rep_structures.items(), key=lambda x: (x[1], x[0].ID), reverse=True))
        else:
            temp = OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True))
        verbose.info("rep_struct_finish",temp)
        return temp

    def check_list_frequency(self, molecules):
        if not self.using_lists:
            return True
        list_membership = [False] * self.list_number

        for molecule in molecules:
            list_membership[molecule.list_member] = True

        return float(list_membership.count(True))/self.list_number >= self.list_threshold

    def _add_structure_to_characteristic(self, structure):
        if not self.cs_begun:
            self.add_CS_base(structure)
            verbose.info("add_to_CS", structure.label)

        elif self.can_add_multiple_times(structure):
            possible_locations = self.find_multiple_addable_locations(structure)
            if possible_locations:
                molecule_map = {}
                verbose.info("add_to_CS", structure.label, len(possible_locations))
                for location in possible_locations:
                    molecule_map.update(self.merge_to_CS(location, structure))
                self._add_cs_locations(structure, molecule_map)

        else:
            possible_locations = []
            #Iterates through all molecules this path struct is in.
            for molecule in self.path_structures[structure]:
                structure_mapping = self.path_structures[structure][molecule]
                self.find_addable_location(structure, molecule, possible_locations, structure_mapping)
            if possible_locations:
                frequency = self.get_location_frequency(possible_locations)
                verbose.info("add_to_CS", structure.label)
                molecule_map = self.merge_to_CS(self.get_k_most_frequent_locations(
                    frequency[0], frequency[1],1)[0], structure)

                self._add_cs_locations(structure, molecule_map)

    def can_add_multiple_times(self, structure):
        if structure not in self.multiple_dict:
            return False
        molecule_dicts = self.multiple_dict[structure]
        return self.isomorphism_factor * len(self.molecules) < len(molecule_dicts)

    def find_multiple_addable_locations(self, structure):

        #Determine how many times the path struct can be added
        molecule_dicts = self.multiple_dict[structure]
        locations = []
        counts = [0]
        for molecule in molecule_dicts:
            struct_mol_map = molecule_dicts[molecule]
            occurrences = len(struct_mol_map.itervalues().next())

            if occurrences >= len(counts):
                counts.extend([0]*(occurrences-len(counts)+1))
            counts[occurrences] +=1
            for idx in xrange(0, occurrences):
                struct_map = {st_vertex:struct_mol_map[st_vertex][idx] for st_vertex in struct_mol_map}
                self.find_addable_location(structure, molecule, locations, struct_map)

        k_value = 0
        for idx, struct_count in enumerate(reversed(counts)):
            if struct_count >= len(self.molecules)*self.isomorphism_factor:
                k_value = idx
                break

        if not locations:
            return None

        #print "loc:",locations
        frequency = self.get_location_frequency(locations)
        #print "freq:",frequency

        return self.get_k_most_frequent_locations(frequency[0], frequency[1], k_value)


    def find_addable_location(self, structure, molecule, locations, structure_mapping):
        if molecule not in self.cs_locations:
            return
        #structure_mapping = ps.Path_Structure.all_structures[structure][molecule]
        cs_mapping = self.cs_locations[molecule]
        cs_path_mapping = {}
        for struct_vertex in structure.vertices():
            # The vertex maps to a molecule vertex which is already represented in the characteristic substructure
            molecule_vertex = structure_mapping[struct_vertex]
            if molecule_vertex in cs_mapping:
                #Maps CS vertices to PS vertices
                #Test if spot is occupied by another structure
                cs_vertex = cs_mapping[molecule_vertex]
                # if molecule.degree(molecule_vertex) != self.characteristic_substructure.degree(cs_vertex):
                #     return
                cs_path_mapping[cs_vertex] = struct_vertex

        if cs_path_mapping:
            locations.append(cs_path_mapping)

    def merge_to_CS(self, location, structure):
        copy_mapping = {}
        connectors = set(location.values())

        #Copy new vertices
        for vertex in structure.vertices():
            if vertex not in connectors:
                cloned_vertex = vertex.clone()
                self.characteristic_substructure.vertex_to_graph(cloned_vertex)
                copy_mapping[vertex] = cloned_vertex

        #Copy non-connecting edges.
        for vertex in structure.vertices():
            if vertex not in connectors:
                for neighbour in structure.neighbours(vertex):
                    if neighbour not in connectors:
                        cloned_edge = structure.get_edge_if_exists(vertex, neighbour).\
                            clone(copy_mapping[vertex],copy_mapping[neighbour])
                        self.characteristic_substructure.edge_to_graph(copy_mapping[vertex],
                                                                       copy_mapping[neighbour], cloned_edge)

        #Add edges from connecting vertices
        for cs_vertex in location:
            struct_vertex = location[cs_vertex]
            for neighbour in structure.neighbours(struct_vertex):
                if neighbour not in connectors:
                    cloned_edge = structure.get_edge_if_exists(struct_vertex, neighbour).\
                        clone(cs_vertex,copy_mapping[neighbour])
                    self.characteristic_substructure.edge_to_graph(cs_vertex,
                                                                   copy_mapping[neighbour], cloned_edge)


        return copy_mapping


    def _add_cs_locations(self, structure, cs_mapping):
        """
        Remembers the locations where molecules map to the characteristic substructure and the structure which is added.

        :param structure: the graph has been recently added to the characteristic substructure
        :return: None
        """

        for molecule in self.path_structures[structure]:
            if molecule not in self.cs_locations:
                self.cs_locations[molecule] = {}
            for vertex in self.path_structures[structure][molecule]:
                # If this vertex of the structure has been added to the CS
                # then add the association to the cs_locations
                if vertex in cs_mapping:
                    self.cs_locations[molecule][self.path_structures[structure][molecule][vertex]] = cs_mapping[vertex]

    def get_location_frequency(self, possible_locations):
        """
        Finds the most frequent location of a structure given a list of graphs displaying possible locations.

        :param possible_locations: list of graphs which combine the CS and a structure
        :return: the graph which is chosen as most frequent
        # """
        buckets = []
        counters = []

        for location in possible_locations:
            not_in_bucket = True
            #Check if bucket contains vertices in question
            for idx, bucket in enumerate(buckets):
                match = True
                for vertex in location.keys():
                    if vertex not in bucket:
                        match = False

                #If so, increment bucket
                if match:
                    counters[idx] +=1
                    not_in_bucket = False

            #Add new bucket with vertex group
            if not_in_bucket:
                new_bucket = {}
                counters.append(1)
                new_bucket.update(location)
                buckets.append(new_bucket)

        return (buckets, counters)

    def get_k_most_frequent_locations(self, buckets, counters, k):
        counts = []
        for idx,count in enumerate(counters):
            counts.append([count, idx*-1]) #Force ordering

        counts = sorted(counts, key=lambda x: (x[0], x[1]), reverse=True)

        bucks = []
        for idx in xrange(0,min(k,len(counts))):
            bucks.append(buckets[counts[idx][1]*-1])

        return bucks

    def get_starting_length(self,paths):
        longest_path_length = len(paths)-1
        if self.length_start and self.length_start <= longest_path_length:
            return self.length_start
        else:
            return longest_path_length

    def add_CS_base(self, structure):
        copy_mapping = {}
        for vertex in structure.vertices():
            cloned_vertex = vertex.clone()
            self.characteristic_substructure.vertex_to_graph(cloned_vertex)
            copy_mapping[vertex] = cloned_vertex
        for vertex in structure.vertices():
            for neighbour in structure.neighbours(vertex):
                cloned_edge = structure.get_edge_if_exists(vertex, neighbour).\
                    clone(copy_mapping[vertex], copy_mapping[neighbour])
                self.characteristic_substructure.edge_to_graph(copy_mapping[vertex],
                                                               copy_mapping[neighbour], cloned_edge)

        self.cs_begun = True
        self._add_cs_locations(structure, copy_mapping)
