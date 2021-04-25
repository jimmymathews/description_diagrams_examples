import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plot

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

def color_by_signature(signature, dimension):
    if len(signature) == dimension:
        return 'green'
    if len(signature) == 1:
        return 'blue'
    return 'red'

def save_preview_visualization(filename, diagram):
    graph = diagram.graph
    signatures = [node for node in graph]
    color_map = [color_by_signature(signature, diagram.dimension) for signature in signatures]
    figure = plot.figure(dpi=200)
    nx.draw(
        graph,
        node_color=color_map,
        ax=figure.add_subplot(111),
        with_labels=True,
        arrowsize=7,
        node_size=125,
        font_size=5,
        pos=nx.fruchterman_reingold_layout(graph, scale=2, k=0.5)
    )
    figure.savefig(filename)

