#!/usr/bin/env python3
import random
import math

import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plot
import numpy as np

from .progress_bar import ProgressBar
from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class DescriptionDiagram:
    def __init__(self, points, progress_bar: ProgressBar=None):
        self.progress_bar = progress_bar
        self.graph = self.create_bipartite_feature_matrix_graph(points)
        self.minimize_cubical_length_in_spanning_class()
        self.save_preview_visualization()

    def create_bipartite_feature_matrix_graph(self, points):
        graph = nx.DiGraph()
        dimension = points.shape[1]
        self.dimension = dimension
        signatures = [[i + (1-int(value))*dimension for i, value in enumerate(row)] for row in points]
        signatures = list(set([tuple(sorted(signature)) for signature in signatures]))
        point_facets = [self.facet(signature) for signature in signatures]
        feature_facets = [self.facet([i]) for i in range(2*dimension)]
        facets = point_facets +  feature_facets
        graph.add_nodes_from([facet.signature for facet in facets])
        for point in point_facets:
            for entry in point.signature:
                feature_signature = tuple(sorted(set([entry])))
                feature = self.facet(feature_signature)
                if point < feature:
                    graph.add_edge(point.signature, feature.signature, weight=1)
                    graph.edges[point.signature, feature.signature]['supports feature inference'] = [(point.signature, feature.signature)]
        return graph

    def facet(self, signature):
        return Facet(sorted(signature), self.dimension)

    def minimize_cubical_length_in_spanning_class(self):
        self.regularize()
        self.do_local_discrete_gradient_descent()
        self.package_results_for_presentation()

    def regularize(self):
        """
        Create maximal chains out of each edge of the original bipartite graph.
        """
        existing_edges = list(self.graph.edges)
        for point_signature, feature_signature in existing_edges:
            if self.progress_bar:
                self.progress_bar.report(
                    expected_total=len(existing_edges),
                    task_description='Creating maximal chains.'
                )
            p = point_signature
            f = feature_signature
            reduced_signature = list(p)
            reduced_signature.remove(f[0])
            permutation = list(f) + list(np.random.permutation(reduced_signature))
            self.graph.remove_edge(p, f)
            for i in range(len(permutation)-1):
                s1 = self.facet(permutation[0:(i+1)]).signature
                s2 = self.facet(permutation[0:(i+2)]).signature
                if self.graph.has_edge(s2, s1):
                    support = self.graph.edges[s2, s1]['supports feature inference']
                    self.graph.edges[s2, s1]['supports feature inference'] = support + [(p, f)]
                else:
                    self.graph.add_edge(s2, s1)
                    self.graph.edges[s2, s1]['supports feature inference'] = [(p, f)]
        logger.info(
            'Weights: %s',
            (list(set([len(self.graph.edges[e[0], e[1]]['supports feature inference']) for e in self.graph.edges]))),
        )

    def do_local_discrete_gradient_descent(self):
        best_option = self.compute_best_basis_homology()
        while best_option['Homology'] is not None:
            self.apply(best_option)
            best_option = self.compute_best_basis_homology()

        isolated = list(nx.isolates(self.graph))
        self.graph.remove_nodes_from(isolated)

    def compute_best_basis_homology(self):
        best_length_change = 0
        best_strand_aggregation_change = 0
        best_homology = None
        null_homologies = []
        ties = []
        triples = self.get_triples()
        for f1, f2, f3 in triples:
            if self.progress_bar:
                self.progress_bar.report(
                    expected_total = len(triples),
                    task_description = 'Considering each basis homology.')
            ss12 = self.get_spanning_support(f1.signature, f2.signature)
            ss23 = self.get_spanning_support(f2.signature, f3.signature)
            ss13 = self.get_spanning_support(f1.signature, f3.signature)
            if len(ss12) + len(ss23) + len(ss13) == 0:
                logger.error(
                    'Triangle considered has no support overlap with the current 1-chain. Weights 12, 23, 13: [%s, %s, %s]',
                    ss12,
                    ss23,
                    ss13,
                )

            # should ensure that only [>1, >1, *]  occurs?
            if len(ss12) == 0 or len(ss23) == 0:
                logger.error(
                    'Wrong triangle type: %s %s %s',
                    ss12,
                    ss23,
                    ss13,
                )

            homology = [
                (f1.signature, f2.signature),
                (f2.signature, f3.signature),
                (f1.signature, f3.signature),
            ]
            length_change, strand_aggregation_change = self.compute_effect(homology)

            if length_change == 0 and strand_aggregation_change != 0:
                null_homologies.append(homology)
            if length_change < best_length_change:
                best_length_change = length_change
                best_strand_aggregation_change = strand_aggregation_change
                ties = [homology]
            elif length_change == best_length_change:
                if strand_aggregation_change < best_strand_aggregation_change:
                    best_strand_aggregation_change = strand_aggregation_change
                    ties = [homology]
                elif strand_aggregation_change == best_strand_aggregation_change:
                    ties.append(homology)

        if best_length_change == 0:
            if len(null_homologies) > 0:
                logger.info(
                    'Using one of the %s null homologies at random.',
                    str(len(null_homologies)),
                )
                best_homology = null_homologies[random.randint(0, len(null_homologies)-1)]
        else:
            logger.info('Success, got better length change: %s', best_length_change)

            if len(ties) > 1:
                best_homology = ties[random.randint(0, len(ties)-1)]
                logger.info('Still, got ties. Choosing from among %s ties randomly.', len(ties))
            elif len(ties) == 1:
                best_homology = ties[0]
            elif len(ties) == 0:
                logger.error('No option found')

        return {
            'Length improvement' : abs(best_length_change),
            'Homology' : best_homology,
        }

    def compute_effect(self, homology):
        length_change = 0

        e12 = homology[0]
        e23 = homology[1]
        e13 = homology[2]
        f1 = self.facet(e12[0])
        f2 = self.facet(e12[1])
        f3 = self.facet(e13[1])

        segment_length_12 = math.sqrt(len(f1.signature) - len(f2.signature))
        segment_length_23 = math.sqrt(len(f2.signature) - len(f3.signature))
        segment_length_13 = math.sqrt(len(f1.signature) - len(f3.signature))

        ss12 = self.get_spanning_support(f1.signature, f2.signature)
        ss23 = self.get_spanning_support(f2.signature, f3.signature)
        ss13 = self.get_spanning_support(f1.signature, f3.signature)

        movable_strands = list(set(ss12).intersection(set(ss23)))
        if len(movable_strands) == 0:
            length_change = 0
            strand_aggregation_change = 0
            return length_change, strand_aggregation_change
        else:
            new_ss12 = list(set(ss12).difference(movable_strands))
            new_ss23 = list(set(ss23).difference(movable_strands))
            new_ss13 = list(set(ss13).union(movable_strands))

        length_change = 0
        if len(new_ss12) == 0:
            length_change += -1*segment_length_12
        if len(new_ss23) == 0:
            length_change += -1*segment_length_23
        if len(ss13) == 0:
            length_change += 1*segment_length_13

        strand_aggregation_before = max([len(ss12), len(ss23), len(ss13)])
        strand_aggregation_after = max([len(new_ss12), len(new_ss23), len(new_ss13)])
        strand_aggregation_change = strand_aggregation_after - strand_aggregation_before

        return length_change, strand_aggregation_change

    # def compute_norm(self):
    #     length = 0
    #     for signature_i, signature_j in self.graph.edges:
    #         codimension = len(signature_i) - len(signature_j)
    #         weight_ij = self.get_weight(signature_i, signature_j)
    #         length += pow(weight_ij, 1/2) * codimension
    #     return length

    def compute_length(self):
        length = 0
        for signature_i, signature_j in self.graph.edges:
            codimension = len(signature_i) - len(signature_j)
            length += math.sqrt(codimension)
        return length

    def get_triples(self):
        triples = []
        edges = list(self.graph.edges)
        for signature1, signature2 in edges:
            f1 = self.facet(signature1)
            f2 = self.facet(signature2)
            if self.progress_bar:
                self.progress_bar.report(expected_total = len(edges), task_description='Building list of basis homology options.')

            # 1-supported case; guided by dim0 support
            # for signature3 in self.graph.nodes:
            #     f3 = self.facet(signature3)
            #     if f1 < f3 and f3 < f2:
            #         F1 = f1
            #         F2 = f3
            #         F3 = f2
            #         triples.append([F1, F2, F3])
            #         logger.info('1-case, going to consider %s %s %s', F1.signature, F2.signature, F3.signature)

            # consistent pair type
            for signature3 in self.graph.successors(signature2):
                f3 = self.facet(signature3)
                ss12 = self.get_spanning_support(f1.signature, f2.signature)
                ss23 = self.get_spanning_support(f2.signature, f3.signature)
                movable_strands = list(set(ss12).intersection(set(ss23)))
                if len(movable_strands) > 0:
                    triples.append([f1, f2, f3])

            # inconsistent sink-center type
            # for signature3 in self.graph.predecessors(signature2):
            #     f3 = self.facet(signature3)
            #     swap = f2.signature
            #     F2 = self.facet(f3.signature)
            #     F3 = self.facet(swap)
            #     F1 = f1
            #     if not (F1 < F3):
            #         logger.error(
            #             'Unexpected facet non-containment: [%s, %s]',
            #             str(F1.signature),
            #             str(F2.signature),
            #         )
            #     if not (F2 < F3):
            #         logger.error(
            #             'Unexpected facet non-containment: [%s, %s]',
            #             str(F2.signature),
            #             str(F3.signature),
            #         )
            #     if not (F2 < F1) and not (F1 < F2):
            #         continue
            #     if F2 < F1:
            #         swap = F2.signature
            #         F2 = self.facet(F1.signature)
            #         F1 = self.facet(swap)
            #     triples.append([F1, F2, F3])

            # inconsistent source-center type
            # for signature3 in self.graph.successors(signature1):
            #     f3 = self.facet(signature3)
            #     if not (f2 < f3) and not (f3 < f2):
            #         continue
            #     if f3 < f2:
            #         swap = f3.signature
            #         F3 = self.facet(f2.signature)
            #         F2 = self.facet(swap)
            #         F1 = f1
            #     else:
            #         F2 = f2
            #         F3 = f3
            #         F1 = f1
            #     triples.append([F1, F2, F3])
        return triples

    def get_spanning_support(self, source, target):
        if not self.graph.has_edge(source, target):
            return []
        else:
            return self.graph.edges[source, target]['supports feature inference']

    def apply(self, basis_option):
        homology = basis_option['Homology']
        e12 = homology[0]
        e23 = homology[1]
        e13 = homology[2]
        f1 = self.facet(e12[0])
        f2 = self.facet(e12[1])
        f3 = self.facet(e13[1])
        s1 = f1.signature
        s2 = f2.signature
        s3 = f3.signature
        ss12 = self.get_spanning_support(f1.signature, f2.signature)
        ss23 = self.get_spanning_support(f2.signature, f3.signature)
        ss13 = self.get_spanning_support(f1.signature, f3.signature)

        movable_strands = list(set(ss12).intersection(set(ss23)))
        if len(movable_strands) == 0:
            logger.error('This homology can not be applied, no strands in common.')

        new_ss12 = list(set(ss12).difference(movable_strands))
        new_ss23 = list(set(ss23).difference(movable_strands))
        new_ss13 = list(set(ss13).union(movable_strands))

        if len(new_ss12) == 0:
            self.graph.remove_edge(s1, s2)
        if len(new_ss23) == 0:
            self.graph.remove_edge(s2, s3)
        if len(new_ss13) == 0:
            self.graph.remove_edge(s1, s3)

        if len(new_ss12) > 0:
            if not self.graph.has_edge(s1, s2):
                self.graph.add_edge(s1, s2)
            self.graph.edges[s1, s2]['supports feature inference'] = new_ss12
        if len(new_ss23) > 0:
            if not self.graph.has_edge(s2, s3):
                self.graph.add_edge(s2, s3)
            self.graph.edges[s2, s3]['supports feature inference'] = new_ss23
        if len(new_ss13) > 0:
            if not self.graph.has_edge(s1, s3):
                self.graph.add_edge(s1, s3)
            self.graph.edges[s1, s3]['supports feature inference'] = new_ss13

        logger.info('Total length: %s', self.compute_length())
        logger.info('Number of edges: %s', str(len(self.graph.edges)))
        logger.info('Number of nodes: %s', str(len(self.graph.nodes)))

    def package_results_for_presentation(self):
        """
        Output.
        """
        # self.calculate_maximal_directed_spanning_tree()
        
    def save_preview_visualization(self):
        signatures = [node for node in self.graph]
        assignment = lambda signature: 'red' if len(signature) > 1 and len(signature) < self.dimension else 'blue' if len(signature) == 1 else 'green' 
        # color_map = ['red' if len(signature) > 1 else 'blue' for signature in signatures]
        color_map = [assignment(signature) for signature in signatures]
        figure = plot.figure(dpi=200)
        # nx.draw_spring(self.graph, node_color=color_map, ax=figure.add_subplot(111), with_labels=True, arrowsize=7, node_size=15, font_size=9)

        nx.draw(
            self.graph,
            node_color=color_map,
            ax=figure.add_subplot(111),
            with_labels=True,
            arrowsize=7,
            node_size=125,
            font_size=5,
            pos=nx.fruchterman_reingold_layout(self.graph, scale=2, k=0.5)
        )
        figure.savefig('graph.png')


import numpy as np
points = np.array([
    [1,0,1,1],
    [0,1,1,1],
    [0,0,1,1],
    [1,1,1,1],
    [0,0,0,0],
    [1,0,0,0],
    [0,1,0,0],
    [1,1,0,0],
])
# points = np.array([
#     [1,1,0],
#     [1,1,1],
# ])
points = np.loadtxt('colon_test_binary_10.csv', delimiter=',')
progress_bar = ProgressBar()
d = DescriptionDiagram(points, progress_bar=progress_bar)

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





