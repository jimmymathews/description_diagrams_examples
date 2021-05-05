import random
import math

from .strand_homology import StrandHomologyPriority, StrandHomologyEffectCalculator
from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class StrandEnrichedHomologyAssessor:
    def __init__(self, diagram):
        self.diagram = diagram
        super().__init__()

    def find_best_homology(self, priority: StrandHomologyPriority=None):
        best_length_change = 0
        best_strand_aggregation_change = 0
        best_homology = None
        null_homologies = []
        ties = []
        triples = self.get_triples()
        for f1, f2, f3 in triples:
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

            if priority == StrandHomologyPriority.TOTAL_LENGTH:
                pair = (length_change, -1*strand_aggregation_change)
                best_pair = (best_length_change, -1*best_strand_aggregation_change)

            if priority == StrandHomologyPriority.AGGREGATION:
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

    def edge_length(self, signature1, signature2):
        return math.sqrt(len(signature1) - len(signature2))

    def compute_effect(self, homology):
        supp, new_supp, s, D = StrandHomologyEffectCalculator.compute_effect(homology, self.diagram)

        I = lambda input_set: 0 if len(input_set) == 0 else 1
        length_changes = [
           (I(new_supp[(i,j)]) - I(supp[(i,j)])) * self.edge_length(s[i],s[j]) for i,j in supp.keys()
        ]
        length_change = sum(length_changes)

        strand_aggregation_before = max([len(entry) for entry in supp.values()])
        strand_aggregation_after = max([len(entry) for entry in new_supp.values()])
        strand_aggregation_change = strand_aggregation_after - strand_aggregation_before

        return length_change, strand_aggregation_change
