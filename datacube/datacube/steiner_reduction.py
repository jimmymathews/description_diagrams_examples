
import networkx as nx

from .diagram_mutation import DescriptionDiagramMutator
from .homology_assessor import StrandEnrichedHomologyAssessor
from .strand_homology import StrandHomologyPriority, StrandHomologyEffectCalculator
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

        2. A "fine" phase the directly aims for length reduction. (Optimize
        length.)

        In both phases, the metric used on the *other* phase is used to break
        ties when deciding between the elementary moves.

        Args:
            diagram (DescriptionDiagram):
                The positive 1-chain on the barycentric subdivision of the
                hypercube to be modified in-place.
        """
        self.optimize(diagram, priority = StrandHomologyPriority.AGGREGATION)
        self.optimize(diagram, priority = StrandHomologyPriority.TOTAL_LENGTH)
        self.declare_task_finished()

    def optimize(self, diagram, priority: StrandHomologyPriority=None):
        diagram_crawler = StrandEnrichedHomologyAssessor(diagram)
        while True:
            best_option = diagram_crawler.find_best_homology(priority=priority)
            if best_option['Homology'] is None:
                break
            self.record_progress(
                expected_total = None,
                task_description = 'Length: ' + str(diagram.compute_length()),
            )
            self.apply_homology(best_option, diagram)
            isolated = list(nx.isolates(diagram.graph))
            diagram.graph.remove_nodes_from(isolated)

    def apply_homology(self, basis_option, diagram):
        homology = basis_option['Homology']
        supp, new_supp, s, D = StrandHomologyEffectCalculator.compute_effect(homology, diagram)
        graph = diagram.graph
        for i, j in D.values():
            if len(new_supp[(i,j)]) == 0:
                graph.remove_edge(s[i], s[j])

        for i, j in D.values():
            if len(new_supp[(i,j)]) > 0:
                if not graph.has_edge(s[i], s[j]):
                    graph.add_edge(s[i], s[j])
                graph.edges[s[i], s[j]]['supports feature inference'] = new_supp[(i,j)]
