import networkx as nx
import networkx.algorithms.isomorphism as iso

import src.molecule_graphs.molecule as m


def isomorphism(pattern, target, nx_structures):

    """
    Uses the NetworkX isomorphism algorithm to check if the pattern graph and the target graph are isomorphic.

    The faster_could_be_isomorphic method is used to discount two structures if they could not be isomorphic.
    :param pattern: a molecule object which is to be tested for isomorphism
    :param target: a molecule object which pattern graph is to be compared against
    :return: None if the graphs are not isomorphic
    :return: a dictionary which maps the indices of the two NetworkX graphs together if they are isomorphic
    """
    if pattern not in nx_structures:
        nx_structures[pattern] = create_nx_graph(pattern, nx_structures)
    if target not in nx_structures:
        nx_structures[target] = create_nx_graph(target, nx_structures)
    if not nx.faster_could_be_isomorphic(nx_structures[pattern], nx_structures[target]):
        # Graphs are definitely not isomorphic
        return None

    #print pattern, target
    # Ensures the isomorphism considers the vertex label and edge type
    matcher = iso.GraphMatcher(nx_structures[pattern], nx_structures[target],
                               node_match=iso.categorical_node_match('label', 'C'),
                               edge_match=iso.categorical_edge_match('type', 'single'))
    if matcher.is_isomorphic():
        return matcher.mapping


def create_nx_graph(structure, nx_structures):
    """
    Creates a copy of the molecule object as an NetworkX graph
    The position attirbute of each structure vertices is used as the index of the NetworkX graph

    :param structure: the molecule object which will be turned into a NetworkX graph
    :return: None
    """
    g = nx.Graph()
    # For each vertex and edge in molecule graph add node and edge in NetworkX graph
    for n in structure.vertices():
        g.add_node(structure.position_of_vertex(n), label=n.label)
    for e in structure.edges():
        if isinstance(e, m.Bond):
            g.add_edge(structure.endpoints_position(e)[0], structure.endpoints_position(e)[1], type=e.type.name)

    return g