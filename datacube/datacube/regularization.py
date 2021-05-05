
import numpy as np

from .diagram_mutation import DescriptionDiagramMutator
from .progress_bar import ProgressingTask
from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class RandomResolutionOfBipartite(DescriptionDiagramMutator, ProgressingTask):
    def __init__(self):
        super().__init__()

    def mutate(self, diagram):
        self.do_regularization(diagram)

    def do_regularization(self, diagram):
        """
        Create maximal chains out of each edge of the original bipartite graph.
        """
        graph = diagram.graph
        existing_edges = list(graph.edges)
        for point_signature, feature_signature in existing_edges:
            self.record_progress(
                expected_total = len(existing_edges),
                task_description = 'Creating maximal chains.',
            )
            p = point_signature
            f = feature_signature
            reduced_signature = list(p)
            reduced_signature.remove(f[0])
            permutation = list(f) + list(np.random.permutation(reduced_signature))
            graph.remove_edge(p, f)
            for i in range(len(permutation)-1):
                s1 = diagram.facet(permutation[0:(i+2)]).signature
                s2 = diagram.facet(permutation[0:(i+1)]).signature
                if graph.has_edge(s1, s2):
                    support = graph.edges[s1, s2]['supports feature inference']
                    graph.edges[s1, s2]['supports feature inference'] = support + [(p, f)]
                else:
                    graph.add_edge(s1, s2)
                    graph.edges[s1, s2]['supports feature inference'] = [(p, f)]
        logger.info(
            'Weights: %s',
            (list(set([len(graph.edges[e[0], e[1]]['supports feature inference']) for e in graph.edges]))),
        )
