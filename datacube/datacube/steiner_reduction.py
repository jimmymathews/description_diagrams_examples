import math

import networkx as nx

from .diagram_mutation import DescriptionDiagramMutator
from .homology_assessor import StrandEnrichedHomologyAssessor
from .strand_homology import StrandHomologyPriority, StrandHomologyEffectCalculator, InapplicableHomology
from .lazy_heap import LazyHeap
from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class SteinerReduction(DescriptionDiagramMutator):
    def __init__(self):
        super().__init__()

    def mutate(self, diagram):
        self.do_local_discrete_gradient_descent(diagram)

    def do_local_discrete_gradient_descent(self, diagram):
        """
        Optimizes the geometric length of a description diagram with respect to
        its realization on the data hypercube.

        The implementation takes place in two phases.

        1. A "coarse" phase intended mainly to aggregate strands together.
        (Optimize a strand aggregation metric.)

        2. A "fine" phase that directly aims for length reduction. (Optimize
        length.)

        In both phases, the metric used on the *other* phase is used to break
        ties when deciding between the elementary moves.

        Args:
            diagram (DescriptionDiagram):
                The positive 1-chain on the barycentric subdivision of the
                hypercube to be modified in-place.
        """
        logger.info('Optimizing for strand aggregation.')
        self.optimize(diagram, priority = StrandHomologyPriority.AGGREGATION)
        logger.info('Optimizing for length.')
        self.optimize(diagram, priority = StrandHomologyPriority.TOTAL_LENGTH)
        logger.info('Finished optimizing.')
        self.declare_task_finished()

    def optimize(self, diagram, priority: StrandHomologyPriority=None):
        self.initialize_elementary_moves(diagram, priority=priority)
        while True:
            continuing = self.take_best_option(diagram, priority=priority)
            if continuing is False:
                break
            self.record_progress(
                expected_total = None,
                task_description = 'Length: ' + str(diagram.compute_length()),
            )

    def initialize_elementary_moves(self, diagram, priority: StrandHomologyPriority=None):
        graph = diagram.graph
        edges = list(graph.edges)
        self.homologies_supported_by_edge = {}

        self.elementary_moves = LazyHeap()
        self.register_triples(self.get_triples_along(diagram), diagram, priority)

    def register_triples(self, triples, diagram, priority, silent=False):
        graph = diagram.graph
        for f1, f2, f3 in triples:
            if not silent:
                self.record_progress(
                    expected_total = len(triples),
                    task_description = 'Considering each 2-simplex for effectivity.',
                )

            ss12 = diagram.get_spanning_support(f1.signature, f2.signature)
            ss23 = diagram.get_spanning_support(f2.signature, f3.signature)
            ss13 = diagram.get_spanning_support(f1.signature, f3.signature)
            if len(ss12) + len(ss23) + len(ss13) == 0:
                logger.error(
                    'Triangle considered has no support overlap with the current 1-chain. Weights 12, 23, 13: [%s, %s, %s]',
                    ss12,
                    ss23,
                    ss13,
                )
            if len(ss12) == 0 or len(ss23) == 0:
                logger.error(
                    'Wrong triangle type: %s %s %s',
                    ss12,
                    ss23,
                    ss13,
                )
            s1 = f1.signature
            s2 = f2.signature
            s3 = f3.signature

            homology = (
                (s1, s2),
                (s2, s3),
                (s1, s3),
            )

            try:
                length_change, strand_aggregation_change = self.compute_metric_changes(diagram, homology)
            except InapplicableHomology as ih:
                continue

            if priority == StrandHomologyPriority.TOTAL_LENGTH:
                priority_object = (length_change, -1*strand_aggregation_change)
            elif priority == StrandHomologyPriority.AGGREGATION:
                priority_object = (-1*strand_aggregation_change, length_change)
            else:
                logger.error('Strand homology priority type %s not supported.', priority.name)
                return


            if length_change >= 0 and priority == StrandHomologyPriority.TOTAL_LENGTH:
                continue
            elif strand_aggregation_change <= 0 and priority == StrandHomologyPriority.AGGREGATION:
                continue

            self.elementary_moves.push_item(homology, priority_object=priority_object)

            if (s1,s2) in self.homologies_supported_by_edge:
                supported = self.homologies_supported_by_edge[(s1,s2)]
            else:
                supported = set([])
            self.homologies_supported_by_edge[(s1,s2)] = supported.union(set([homology]))

            if (s2,s3) in self.homologies_supported_by_edge:
                supported = self.homologies_supported_by_edge[(s2,s3)]
            else:
                supported = set([])
            self.homologies_supported_by_edge[(s2,s3)] = supported.union(set([homology]))

            if (s1,s3) in self.homologies_supported_by_edge:
                supported = self.homologies_supported_by_edge[(s1,s3)]
            else:
                supported = set([])
            self.homologies_supported_by_edge[(s1,s3)] = supported.union(set([homology]))

    def get_triples_along(self, diagram):
        triples = []
        graph = diagram.graph
        edges = list(graph.edges)
        for signature1, signature2 in edges:
            self.record_progress(
                expected_total = len(edges),
                task_description = 'Searching for 2-simplices along the given 1-chain. ',
            )
            f1 = diagram.facet(signature1)
            f2 = diagram.facet(signature2)
            for signature3 in graph.successors(signature2):
                f3 = diagram.facet(signature3)
                ss12 = diagram.get_spanning_support(f1.signature, f2.signature)
                ss23 = diagram.get_spanning_support(f2.signature, f3.signature)
                movable_strands = list(set(ss12).intersection(set(ss23)))
                if len(movable_strands) > 0:
                    triples.append([f1, f2, f3])
        return triples

    def compute_metric_changes(self, diagram, homology):
        supp, new_supp, s, D = StrandHomologyEffectCalculator.compute_effect(homology, diagram)
        I = lambda input_set: 0 if len(input_set) == 0 else 1
        length_changes = [
           (I(new_supp[(i,j)]) - I(supp[(i,j)])) * self.edge_length(s[i],s[j]) for i,j in supp.keys()
        ]
        length_change = sum(length_changes)
        strand_aggregation_before = max([len(entry) for entry in supp.values()])
        strand_aggregation_after = max([len(entry) for entry in new_supp.values()])
        strand_aggregation_change = strand_aggregation_after - strand_aggregation_before
        return length_change, strand_aggregation_change

    def edge_length(self, signature1, signature2):
        return math.sqrt(len(signature1) - len(signature2))

    def get_best_option(self):
        homology = self.elementary_moves.pop()
        return homology

    def take_best_option(self, diagram, priority):
        homology = self.get_best_option()
        if homology is None:
            return False

        try:
            length_change, strand_aggregation_change = self.compute_metric_changes(diagram, homology)
            if priority == StrandHomologyPriority.TOTAL_LENGTH:
                if length_change > 0:
                    logger.debug('Met last good option by length.')
                    return False
            if priority == StrandHomologyPriority.AGGREGATION:
                if strand_aggregation_change < 0:
                    logger.debug('Met last good option by aggregation.')
                    return False
        except InapplicableHomology as ih:
            logger.debug('Best option was no good: %s', homology)
            return False

        supp, new_supp, s, D = StrandHomologyEffectCalculator.compute_effect(homology, diagram)

        graph = diagram.graph
        for i, j in D.values():
            if len(new_supp[(i,j)]) == 0:
                if graph.has_edge(s[i], s[j]):
                    graph.remove_edge(s[i], s[j])

        for i, j in D.values():
            if len(new_supp[(i,j)]) > 0:
                if not graph.has_edge(s[i], s[j]):
                    graph.add_edge(s[i], s[j])
                graph.edges[s[i], s[j]]['supports feature inference'] = new_supp[(i,j)]
        isolated = list(nx.isolates(diagram.graph))
        diagram.graph.remove_nodes_from(isolated)

        self.update_heap(homology, diagram, priority)

        return True

    def update_heap(self, homology, diagram, priority):
        graph = diagram.graph
        triples = []
        for si, sj in homology:
            if (si,sj) in self.homologies_supported_by_edge:
                homologies = self.homologies_supported_by_edge[(si,sj)]

                for homology in homologies:
                    self.elementary_moves.remove_item(homology)

                s1 = si
                s2 = sj
                f1 = diagram.facet(s1)
                f2 = diagram.facet(s2)
                if s2 in graph.nodes:
                    for s3 in graph.successors(s2):
                        f3 = diagram.facet(s3)
                        ss12 = diagram.get_spanning_support(f1.signature, f2.signature)
                        ss23 = diagram.get_spanning_support(f2.signature, f3.signature)
                        movable_strands = list(set(ss12).intersection(set(ss23)))
                        if len(movable_strands) > 0:
                            triples.append([f1, f2, f3])
                if s1 in graph.nodes:
                    for s0 in graph.predecessors(s1):
                        f0 = diagram.facet(s0)
                        ss01 = diagram.get_spanning_support(f0.signature, f1.signature)
                        ss12 = diagram.get_spanning_support(f1.signature, f2.signature)
                        movable_strands = list(set(ss01).intersection(set(ss12)))
                        if len(movable_strands) > 0:
                            triples.append([f0, f1, f2])

        self.register_triples(triples, diagram, priority, silent=True)



        # consider each of the 3 legs of support of the homology
        # given a leg,
        # look up all heap objects that contain that leg in their support
        # for each heap object,
        # look at the homology it has
        # calculate the metric changes it would have *now*
        #    catch if it is no longer applicable, and remove it from the queue
        # update the priority ojbect based on the new metric changes

        # then:
        # still given a leg
        # find all new supported homologies with consecutive legs in the diagram's support
        # push them into the queue









