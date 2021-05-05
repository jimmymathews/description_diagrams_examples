import random
import math

import networkx as nx
import numpy as np

from .diagram_mutation import DescriptionDiagramMutator
from .progress_bar import ProgressingTask
from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class StrandEnrichedHomologyAssessor(ProgressingTask):
    def __init__(self, diagram):
        self.diagram = diagram
        super().__init__()

    def find_best_homology(self, priority='aggregation'):
        best_length_change = 0
        best_strand_aggregation_change = 0
        best_homology = None
        null_homologies = []
        ties = []
        triples = self.get_triples()
        for f1, f2, f3 in triples:
            # self.record_progress(
            #     expected_total = len(triples),
            #     task_description = 'Considering each basis homology.',
            # )
            ss12 = self.diagram.get_spanning_support(f1.signature, f2.signature)
            ss23 = self.diagram.get_spanning_support(f2.signature, f3.signature)
            ss13 = self.diagram.get_spanning_support(f1.signature, f3.signature)
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

            if priority == 'length':
                pair = (length_change, -1*strand_aggregation_change)
                best_pair = (best_length_change, -1*best_strand_aggregation_change)

            if priority == 'aggregation':
                pair = (-1*strand_aggregation_change, length_change)
                best_pair = (-1*best_strand_aggregation_change, best_length_change)

            if pair < best_pair:
                best_homology = homology
                best_length_change = length_change
                best_strand_aggregation_change = strand_aggregation_change

            # if length_change < best_length_change:
            #     best_length_change = length_change
            #     best_strand_aggregation_change = strand_aggregation_change
            #     ties = [homology]
            # elif length_change == best_length_change:
            #     if strand_aggregation_change < best_strand_aggregation_change:
            #         best_strand_aggregation_change = strand_aggregation_change
            #         ties = [homology]
            #     elif strand_aggregation_change == best_strand_aggregation_change:
            #         ties.append(homology)

        if best_length_change == 0:
            if len(null_homologies) > 0:
                logger.info(
                    'Using one of the %s null homologies at random.',
                    str(len(null_homologies)),
                )
                best_homology = null_homologies[random.randint(0, len(null_homologies)-1)]
        # else:
        #     logger.info('Success, got better length change: %s', best_length_change)
            # if len(ties) > 1:
            #     best_homology = ties[random.randint(0, len(ties)-1)]
            #     logger.info('Still, got ties. Choosing from among %s ties randomly.', len(ties))
            # elif len(ties) == 1:
            #     best_homology = ties[0]
            # elif len(ties) == 0:
            #     logger.error('No option found')

        return {
            'Length improvement' : abs(best_length_change),
            'Homology' : best_homology,
        }

    def get_triples(self):
        triples = []
        graph = self.diagram.graph
        edges = list(graph.edges)
        for signature1, signature2 in edges:
            f1 = self.diagram.facet(signature1)
            f2 = self.diagram.facet(signature2)
            # self.record_progress(
            #     expected_total = len(edges),
            #     task_description='Building list of basis homology options.',
            # )
            for signature3 in graph.successors(signature2):
                f3 = self.diagram.facet(signature3)
                ss12 = self.diagram.get_spanning_support(f1.signature, f2.signature)
                ss23 = self.diagram.get_spanning_support(f2.signature, f3.signature)
                movable_strands = list(set(ss12).intersection(set(ss23)))
                if len(movable_strands) > 0:
                    triples.append([f1, f2, f3])
        return triples

    def compute_effect(self, homology):
        length_change = 0
        e12 = homology[0]
        e23 = homology[1]
        e13 = homology[2]
        f1 = self.diagram.facet(e12[0])
        f2 = self.diagram.facet(e12[1])
        f3 = self.diagram.facet(e13[1])
        segment_length_12 = math.sqrt(len(f1.signature) - len(f2.signature))
        segment_length_23 = math.sqrt(len(f2.signature) - len(f3.signature))
        segment_length_13 = math.sqrt(len(f1.signature) - len(f3.signature))
        ss12 = self.diagram.get_spanning_support(f1.signature, f2.signature)
        ss23 = self.diagram.get_spanning_support(f2.signature, f3.signature)
        ss13 = self.diagram.get_spanning_support(f1.signature, f3.signature)
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


class SteinerReduction(DescriptionDiagramMutator):
    def __init__(self):
        super().__init__()

    def mutate(self, diagram):
        self.L0 = diagram.compute_length()
        self.E0 = len(diagram.graph.edges)
        self.minimize_cubical_length_in_spanning_class(diagram)
        isolated = list(nx.isolates(diagram.graph))
        diagram.graph.remove_nodes_from(isolated)

    def minimize_cubical_length_in_spanning_class(self, diagram):
        self.do_local_discrete_gradient_descent(diagram)

    def do_local_discrete_gradient_descent(self, diagram):
        diagram_crawler = StrandEnrichedHomologyAssessor(diagram)
        for listener in self.get_progress_listeners():
            diagram_crawler.add_progress_listener(listener)

        best_option = diagram_crawler.find_best_homology(priority='aggregation')
        while best_option['Homology'] is not None:
            self.record_progress(
                expected_total = None,
                task_description = str(diagram.compute_length()),
            )
            self.apply(best_option, diagram)
            best_option = diagram_crawler.find_best_homology(priority='aggregation')

        best_option = diagram_crawler.find_best_homology(priority='length')
        while best_option['Homology'] is not None:
            self.apply(best_option, diagram)
            best_option = diagram_crawler.find_best_homology(priority='length')

    def apply(self, basis_option, diagram):
        homology = basis_option['Homology']
        e12 = homology[0]
        e23 = homology[1]
        e13 = homology[2]
        f1 = diagram.facet(e12[0])
        f2 = diagram.facet(e12[1])
        f3 = diagram.facet(e13[1])
        s1 = f1.signature
        s2 = f2.signature
        s3 = f3.signature
        ss12 = diagram.get_spanning_support(f1.signature, f2.signature)
        ss23 = diagram.get_spanning_support(f2.signature, f3.signature)
        ss13 = diagram.get_spanning_support(f1.signature, f3.signature)
        movable_strands = list(set(ss12).intersection(set(ss23)))
        if len(movable_strands) == 0:
            logger.error('This homology can not be applied, no strands in common.')
        new_ss12 = list(set(ss12).difference(movable_strands))
        new_ss23 = list(set(ss23).difference(movable_strands))
        new_ss13 = list(set(ss13).union(movable_strands))
        graph = diagram.graph
        if len(new_ss12) == 0:
            graph.remove_edge(s1, s2)
        if len(new_ss23) == 0:
            graph.remove_edge(s2, s3)
        if len(new_ss13) == 0:
            graph.remove_edge(s1, s3)
        if len(new_ss12) > 0:
            if not graph.has_edge(s1, s2):
                graph.add_edge(s1, s2)
            graph.edges[s1, s2]['supports feature inference'] = new_ss12
        if len(new_ss23) > 0:
            if not graph.has_edge(s2, s3):
                graph.add_edge(s2, s3)
            graph.edges[s2, s3]['supports feature inference'] = new_ss23
        if len(new_ss13) > 0:
            if not graph.has_edge(s1, s3):
                graph.add_edge(s1, s3)
            graph.edges[s1, s3]['supports feature inference'] = new_ss13
