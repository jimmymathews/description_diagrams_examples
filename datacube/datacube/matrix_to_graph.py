
import networkx as nx

from .data_structures import DescriptionDiagram

def bipartite_diagram_representation_of_feature_matrix(points):
    graph = nx.DiGraph()
    dimension = points.shape[1]
    feature_index = lambda index, value: index + (1-int(value))*dimension
    signatures_unsorted = [[feature_index(i, v) for i, v in enumerate(row)] for row in points]
    signatures = list(set([tuple(sorted(signature)) for signature in signatures_unsorted]))
    diagram = DescriptionDiagram(dimension = dimension)
    point_facets = [diagram.facet(signature) for signature in signatures]
    feature_facets = [diagram.facet((j,)) for j in range(2*dimension)]
    facets = point_facets + feature_facets
    graph.add_nodes_from([facet.signature for facet in facets])
    for point in point_facets:
        for entry in point.signature:
            feature_signature = tuple(sorted(set([entry])))
            feature = diagram.facet(feature_signature)
            if point < feature:
                p = point.signature
                f = feature.signature
                graph.add_edge(p, f)
                p_is_in_f = (p, f)
                graph.edges[p, f]['supports feature inference'] = [p_is_in_f]
    diagram.set_graph(graph)
    return diagram
