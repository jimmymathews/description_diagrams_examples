import math

import networkx as nx

from .log_formats import colorized_logger
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

    def __eq__(self, other):
        return self.signature == other.signature

    def __lt__(self, other):
        return set(other.signature) < set(self.signature)

    def __gt__(self, other):
        return set(self.signature) < set(other.signature)


class DescriptionDiagram:
    def __init__(self, dimension):
        self.dimension = dimension
        self.graph = None

    def set_graph(self, graph):
        self.graph = graph

    def facet(self, signature):
        return Facet(sorted(signature), self.dimension)

    def get_spanning_support(self, source, target):
        if not self.graph.has_edge(source, target):
            return []
        else:
            return self.graph.edges[source, target]['supports feature inference']

    def compute_length(self):
        length = 0
        for signature_i, signature_j in self.graph.edges:
            codimension = len(signature_i) - len(signature_j)
            length += math.sqrt(codimension)
        return length
