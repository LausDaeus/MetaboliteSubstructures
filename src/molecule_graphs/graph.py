#from draw_molecule import draw_molecule as draw
import Queue
"""
This module contains a graph implemented with an adjacency dictionary
and the edge and vertex objects which can be added to the graph.
Based on the Graph ADT in 'Algorithms and Data Structures in Python'
by  Michael T. Goodrich, Roberto Tamassia, Michael H. Goldwasser
"""


class Vertex(object):
    """
    A single vertex of a graph which includes its label.

    The vertex does not hold information about its connections; this information is stored in the graph
    :param label: a string which is the label of the vertex
    position: an integer which is used to signify the order in which the vertex was added to the graph
            it can be used as a simple unique identifier for the vertex within the graph
    visited: a boolean flag which is used in the find_all_paths method of a graph
    """
    ID_Counter = 0
    def __init__(self, label):
        self.label = label
        self.visited = False
        self.ID = Vertex.ID_Counter
        Vertex.ID_Counter +=1

    def __hash__(self):
        return hash(id(self))


class Edge(object):
    """
    A single Edge of a Graph which includes the two vertices that it connects and its label.

    :param origin: one of the vertex objects which the edge connects
    :param destination: the other vertex which the edge connects
    :param label: a string which is the label of the edge
    """
    def __init__(self, origin, destination, label=None):
        self._origin = origin
        self._destination = destination
        self.label = label
        self.visited = False

    @property
    def endpoints(self):
        """
        Gives the endpoints of the edge
        :return: a tuple of the two endpoints
        """
        return self._origin, self._destination

    def opposite(self, vertex):
        """
        Gives the vertex object which lies opposite to the given vertex on the edge
        :param vertex: the vertex object whose opposite is to be found
        :return: the vertex object that is the other endpoint in the edge
        """
        if vertex == self._origin:
            return self._destination
        if vertex == self._destination:
            return self._origin
        else:
            return None

    def __hash__(self):
        return hash((self._origin, self._destination))

    def __str__(self):
        return 'join %s and %s' % (self._origin, self._destination)


class Graph(object):
    """
    A Graph which contains an adjacency dictionary that contains the information about its vertices and edges.

    adjacency_dictionary: a dictionary that takes the form of {vertex: {adjacent vertex: connecting edge}}
    size: an integer gives the number of individual vertices in the graph
    paths: a mapping of a path (string) to a list of vertices in the graph.
    """
    ID_counter = 0
    def __init__(self):
        self.adjacency_dictionary = {}
        self.size = 0
        self.paths = {}
        # list of vertices, index in list acts as position indicator for isomorphism checking
        self.positions = []
        self.ID = Graph.ID_counter
        Graph.ID_counter += 1

    def vertices(self):
        """
        Returns all the vertices of the graph

        Performance: O(n)
        :return: a list of all the vertex objects present in the graph
        """
        #return self.adjacency_dictionary.keys()
        return self.positions

    def vertex_in_graph(self, vertex):
        return vertex in self.adjacency_dictionary

    def edges(self):
        """
        Returns all the edges of the graph in a set to remove duplicates

        Performance: O(m)
        :return: a set of all the edge objects present in the graph
        """
        edges = set()
        for adjacentVertex in self.adjacency_dictionary.values():
            edges.update(adjacentVertex.values())
        return edges

    def clear(self):
        """
        Clears all the parameters of the graph

        Performance: O(1)
        :return: None
        """
        self.adjacency_dictionary = {}
        self.size = 0
        self.paths = {}
        self.positions = []

    def add_vertex(self, label):
        """
        Creates a new vertex object and adds it to the graph using the vertex_to_graph method.

        Performance: O(1)
        :param label: a string representing the label of the vertex
        :return: the vertex object which has been newly created
        """
        new_vertex = Vertex(label)
        self.vertex_to_graph(new_vertex)
        return new_vertex

    def vertex_to_graph(self, vertex):
        """
        Adds a vertex object to the graph by assigning a dictionary which will contain all adjacent vertices and edges.

        Performance: O(1)
        :param vertex: the vertex object that is to be added to the graph
        :return: None
        """
        self.adjacency_dictionary[vertex] = {}
        self.positions.append(vertex)
        self.size += 1

    def remove_vertex(self, vertex):
        """
        Delete the vertex object from the graph by removing it from the dictionaries of vertices which are adjacent.

        Performance: O(degree of vertex)
        :param vertex: the vertex object that is to be removed
        :return: None
        """
        for neighbour in self.neighbours(vertex):
            del self.adjacency_dictionary[neighbour][vertex]
        del self.adjacency_dictionary[vertex]
        self.positions.remove(vertex)

    def add_edge(self, first_vertex, second_vertex, label=None):
        """
        Create an edge object and add it to the graph using the edge_to_graph method.

        Performance: O(1)
        :param first_vertex: a vertex object that is to be an endpoint of the edge
        :param second_vertex: a vertex object that is to be an endpoint of the edge
        :param label: the label of the edge object
        :return: the edge object which has been newly added to the graph
        """
        new_edge = Edge(first_vertex, second_vertex, label)
        self.edge_to_graph(first_vertex, second_vertex, new_edge)
        return new_edge

    def edge_to_graph(self, first_vertex, second_vertex, edge):
        """
        Adds an edge object to the graph by adding it to the adjacency dictionary entries for both of its endpoints.

        Performance: O(1)
        :param first_vertex: one of the endpoints of the edge that is being added
        :param second_vertex: one of the endpoints of the edge that is being added
        :param edge: the edge that is to be added to the graph
        :return: None
        """
        # Try for key error
        self.adjacency_dictionary[first_vertex][second_vertex] = edge
        self.adjacency_dictionary[second_vertex][first_vertex] = edge

    def remove_edge(self, first_vertex, second_vertex):
        """
        Remove the edge that is found between the two given vertices.

        Performance: O(1)
        :param first_vertex:
        :param second_vertex:
        :return:
        """
        del self.adjacency_dictionary[first_vertex][second_vertex]
        del self.adjacency_dictionary[second_vertex][first_vertex]

    def neighbours(self, vertex):
        """
        Returns the vertices adjacent to the given vertex in the graph

        :param vertex: the vertex whose neighbours are to be found
        :return: a list of vertex objects adjacent to the vertex
        """
        return self.adjacency_dictionary[vertex].keys()

    def connecting_edges(self, vertex):
        """
        Returns the edges which are attached to the given vertex in the graph

        :param vertex: the vertex whose attached edges are to be found
        :return: a list of edge objects which are attached to the given vertex
        """
        return self.adjacency_dictionary[vertex].values()

    def degree(self, vertex):
        """
        Returns the degree of the given vertex which is the number of edges that are attached to the vertex

        :param vertex: the vertex object whose degree is to be returned
        :return: an integer which is the number of edges attached to the vertex
        """
        return len(self.adjacency_dictionary[vertex])

    def position_of_vertex(self, vertex):
        """
        Returns the position of the given vertex

        :param vertex: the vertex object whose position in the positions list is to be returned
        :return: an integer which is the position of the vertex in the graph
        """
        return self.positions.index(vertex)

    def vertex_from_position(self, position):
        """
        Tests if there is a vertex at that position and return it if there is

        :param position: an integer which is to be tested
        :return: a vertex if it had has the given position
        """
        if len(self.positions) > position:
            return self.positions[position]
        else:
            return False

    def endpoints_position(self, edge):
        """
        Gives the positions of the endpoints of the edge
        :return: a tuple of the positions of the endpoints
        """
        endpoints = edge.endpoints
        return self.position_of_vertex(endpoints[0]), self.position_of_vertex(endpoints[1])

    def get_edge_if_exists(self, first_vertex, second_vertex):
        """
        Test if there is an edge joining two vertices in the graph

        :param first_vertex: a vertex that is to be tested for a connection
        :param second_vertex: a vertex that is to be tested for a connection
        :return: an edge object which joins the two given vertices
        """
        if second_vertex in self.adjacency_dictionary[first_vertex]:
            return self.adjacency_dictionary[first_vertex][second_vertex]
        else:
            return False

    def find_all_paths(self, path_list):
        """
        Finds all the possible paths in a graph using a depth first search

        The depth first search is carried out starting from each vertex of the graph so all path combinations are found
        Algorithm structure from Handbook of Graph Theory, Gross & Yellen
        :return: a dictionary containing the strings representing the paths and their lengths as the value
        """

        for v in self.vertices():
            path_stack = []             # Used to create the string of the path
            vertex_stack = []         # Used to create the list of vertices
            self._path_visit(v, path_stack, vertex_stack, path_list)

        #print [len(all_paths[bucket]) for bucket in all_paths]
        #return all_paths

    def _path_visit(self, v, path_stack, vertex_stack, path_list):
        """
        The recursive label of the find_all_paths method.
        As each new path is found the string representing it is created using the path_stack.
        The vertices which are present in the path are stored along with the path as a tuple in self.paths
        :param v: the vertex object which is to be added to the path
        :param path_stack: The stack which is used to construct the path string
        :param vertex_stack: The stack which is used to store the vertices encountered in order
        :param all_paths: A dictionary which is used to store a non-duplicated list of paths along with their lengths
        :return: None
        """
        v.visited = True
        path_stack.append(v.label)
        vertex_stack.append(v)
        for w in self.neighbours(v):
            if not w.visited:
                e = self.get_edge_if_exists(v, w)
                path_stack.append(e.label)
                self._path_visit(w, path_stack, vertex_stack, path_list)

        path = ''.join(path_stack)
        length = len(path)/2 + 1

        while len(path_list) <= length:
            path_list.append(set())

        positions = list(vertex_stack)
        #print length, len(path_list)
        if not path[::-1] in path_list[length] or path[::-1] == path:
            path_list[length].add(path)

            if not path in self.paths:
                self.paths[path] = []
                self.paths[path].append(positions)
            else:
                if not self._is_same_path(self.paths[path], positions):
                    self.paths[path].append(positions)

        # Once the vertex has been processed it is popped from the stack to create the string and to store the vertices
        path_stack.pop()
        if path_stack:
            path_stack.pop() #Remove edge label if non-empty
        vertex_stack.pop()

        v.visited = False

    def _is_same_path(self, paths, given_path):
        for vertices in paths:

            #Check if same path
            same_path = True
            for vertex, position in zip(vertices, given_path):
                if not (vertex is position):
                    same_path = False
                    break

            if same_path:
                return True

            #Check if reverse path
            reverse_path = True
            for vertex, position in zip(reversed(vertices),given_path):
                if not (vertex is position):
                    reverse_path = False
                    break

            if reverse_path:
                return True

        return False
    def create_spanning_tree(self):
        """
        Creates a spanning tree by traversing the molecule in a BFS and setting edge attributes.
        When finished with the tree, call "clear_spanning_tree" to remove any excess data.
        """


        #If there are no vertices in the graph is empty
        if not self.positions:
            return None

        for vertex in self.positions:
            vertex.visited = False

        for edge in self.edges():
            edge.visited = False

        starting_vertices = []
        queue = Queue.Queue()
        root = self.next_unvisited_vertex()

        while root:
            queue.put(root)
            starting_vertices.append(root)
            root.visited = True

            while not queue.empty():
                current = queue.get()

                for neighbour in self.neighbours(current):
                    edge = self.adjacency_dictionary[current][neighbour]
                    if not neighbour.visited:
                        neighbour.visited = True
                        edge.visited = True
                        queue.put(neighbour)

            root = self.next_unvisited_vertex()

        return starting_vertices

    def clean_spanning_tree(self):
        for vertex in self.positions:
            vertex.visited = False

        for edge in self.edges():
            if hasattr(edge, "visited"):
                del edge.visited

    def next_unvisited_vertex(self):
        for vertex in self.vertices():
            if not vertex.visited:
                return vertex
        return None