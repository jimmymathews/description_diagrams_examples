#!/usr/bin/env python3

import igraph
from igraph import Graph

# import matplotlib
# matplotlib.use('agg') # no UI backend

from log_formats import colorized_logger
logger = colorized_logger(__name__)


class Facet:
    def __init__(self, signature, dimension):
        signature = tuple(signature)
        if tuple(sorted(signature)) != signature:
            logger.warning('Use sorted integer lists as signatures for facets. Found: %s', str(signature))
        if any([feature_index >= 2*dimension or feature_index < 0 for feature_index in signature]):
            logger.error('Signature (%s) is malformed for the given dimension (%s)', str(signature), str(dimension))
            raise ValueError()
        self.signature = tuple(signature)
        self.dimension = dimension

    def __hash__(self):
        return self.signature

    def __eq__(self, other):
        return self.signature == other.signature

    def __lt__(self, other):
        return [set(other.signature).issubset(self.signature)]

    def __gt__(self, other):
        return [set(self.signature).issubset(other.signature)]


class DescriptionDiagram:
    def __init__(self, points):
        self.graph = Graph(directed=True)
        self.create_bipartite_feature_matrix_graph(points)

    def create_bipartite_feature_matrix_graph(self, points):
        dimension = points.shape[1]
        self.dimension = dimension
        signatures = [[i + (1-int(value))*dimension for i, value in enumerate(row)] for row in points]
        signatures = list(set([tuple(sorted(signature)) for signature in signatures]))
        f = self.facet
        s = self.signature
        point_facets = [f(signature) for signature in signatures]
        feature_facets = [f([i]) for i in range(2*dimension)]
        facets = point_facets +  feature_facets
        self.graph.add_vertices(len(facets))
        self.graph.vs['signature'] = [s(facet) for facet in facets]
        self.graph.vs['color'] = ['blue' for point in point_facets] + ['red' for feature in feature_facets]
        for point in point_facets:
            signature = s(point)
            for entry in signature:
                feature_signature = tuple(sorted(set([entry])))
                feature = f(feature_signature)
                if point < feature:
                    V_point = self.vertex(signature, read_only=True)
                    V_feature = self.vertex(feature_signature, read_only=True)
                    if (self.graph.get_eid(V_point, V_feature, directed=True, error=False) == -1):
                        self.graph.add_edge(V_point, V_feature)

        # labels = [str(s(facet)) for facet in facets]
        # igraph.plot(self.graph, target='plot.pdf', layout=self.graph.layout_kamada_kawai(), vertex_label=labels)

    def vertex(self, signature, read_only=False):
        """
        Args:
            signature (tuple):
                Tuple of integers representing the signature of a facet which is a node in the graph.
            read_only (bool):
                If set, will not attempt to add the vertex in case it does not yet exist in the graph.
                Otherwise, and this is the default behavior, the vertex will be automatically created.

        Returns:
            v (igraph.Vertex):
                The vertex with the given signature.
        """
        subset = self.graph.vs.select(signature=signature)
        if len(subset) > 0:
            return subset[0]
        else:
            if read_only:
                logger.debug(
                    'You requested vertex by signature (%s) in read-only mode, but the vertex was not found.',
                    str(signature),
                )
                return None
            L = len(self.graph.vs)
            self.graph.add_vertex(1)
            V = self.graph.vs[L]
            V['signature'] = signature
            return V

    def signature(self, facet):
        return facet.signature

    def facet(self, signature):
        return Facet(signature, self.dimension)

    def minimize_cubical_length_in_homology_class(self):
        self.regularize()
        self.do_local_discrete_gradient_descent()
        self.package_results_for_presentation()

    def regularize(self):
        """
        Create maximal chains out of each edge of the original bipartite graph.
        """

    def do_local_discrete_gradient_descent(self):
        best_option = self.compute_best_basis_homology()
        while best_option:
            self.apply(best_option)
            best_option = self.compute_best_basis_homology()

    def compute_best_basis_homology(self):
        return None

    def apply(self, basis_option):
        pass

    def package_results_for_presentation(self):
        self.calculate_maximal_directed_spanning_tree()
        


import numpy as np
points = np.array([[1,0,1,1], [0,1,1,1], [0,0,1,1], [1,0,1,1], [1,0,1,1], [1,0,1,1]])
d = DescriptionDiagram(points)

# g = Graph(directed=True)
# g.add_vertices(3)
# g.add_edges([(0,1), (1,2), (0,2)])
# f1 = Facet([0,1], 3)
# f2 = Facet([0,1,4], 3)
# f3 = Facet([1,3], 3)
# g.vs['facet'] = [f1.signature, f2.signature, f3.signature]
# V1 = g.vs.select(facet = f1.signature)[0]
# V2 = g.vs.select(facet = f2.signature)[0]
# V3 = g.vs.select(facet = f3.signature)[0]

# V1.degree(mode='in')
# V1.degree(mode='out')
# V1.degree()

# V1.in_edges()
# V1.out_edges()

# for V in g.vs:
#     x = V.index
#     print(str(x))
#     for E in V.in_edges():
#         W = E.source_vertex
#         y = W.index
#         print('   ' + str(x) + ' <== ' + str(y))
#         for Ep in W.in_edges():
#             Wp = Ep.source_vertex
#             z = Wp.index
#             if x == z:
#                 continue
#             print('      ' + str(x) + ' <== ' + str(y) + ' <== ' + str(z))
#         for Ep in W.out_edges():
#             Wp = Ep.target_vertex
#             z = Wp.index
#             if x == z:
#                 continue
#             print('      ' + str(x) + ' <== ' + str(y) + ' ==> ' + str(z))

#     for E in V.out_edges():
#         W = E.target_vertex
#         y = W.index
#         print('   ' + str(x) + ' ==> ' + str(y))
#         for Ep in W.in_edges():
#             Wp = Ep.source_vertex
#             z = Wp.index
#             if x == z:
#                 continue
#             print('      ' + str(x) + ' ==> ' + str(y) + ' <== ' + str(z))
#         for Ep in W.out_edges():
#             Wp = Ep.target_vertex
#             z = Wp.index
#             if x == z:
#                 continue
#             print('      ' + str(x) + ' ==> ' + str(y) + ' ==> ' + str(z))





