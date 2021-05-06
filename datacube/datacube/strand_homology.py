import enum
from enum import Enum, auto

from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class StrandHomologyPriority(Enum):
    """
    These options specify whether a strand homology (typically an atomic
    one, consisting of a single 2-simplex boundary), considered as applying
    to a "strand-enriched" 1-chain, should be selected based on amount of
    increase of strand aggregation, or based on the the amount of reduction
    of geometric length. 
    """
    AGGREGATION = auto()
    TOTAL_LENGTH = auto()


class StrandHomologyEffectCalculator:
    def compute_effect(homology, diagram):
        """
        Computes the effect of applying an elementary 1-homology on strand
        routing in a strand-enriched description diagram.

        Args:
            homology (tuple):
                A 3-element tuple, each element of which is a pair of facet
                signatures. Interpreted as a 2-simplex boundary.

            diagram (DescriptionDiagram):
                The diagram context in which the homology application is to be
                considered.

        Returns:
            support (dict):
                For each of the 3 edges of the homology, the set of strands
                (i.e. (point, feature-plane) pairs) which are currently routed
                through that edge in the diagram.

            new_support (dict):
                The support (as described above) that would result from
                application of the elementary homology.

            signatures (dict):
                For each vertex (indicated with an integer 0, 1, 2 in reference
                to the vertices in the standard form of the given 2-simplex
                homology), the facet signature it corresponds to. This is
                provided for convenience.

            edge_order (dict):
                The (standard) order of the edges in the homology. This is provided
                for convenience.
        """
        V = [0,1,2]
        D = {
            0 : (0,1),
            1 : (1,2),
            2 : (0,2),
        }
        e = {D[i] : homology[i] for i in V}
        f = {i : diagram.facet(StrandHomologyEffectCalculator.get_simplex_vertex(i, e)) for i in V}
        s = {i : f[i].signature for i in V}
        supp = {(i,j) : diagram.get_spanning_support(s[i], s[j]) for i, j in D.values()}
        movable_strands = list(
            set(supp[(0,1)]).intersection(set(supp[(1,2)]))
        )
        if len(movable_strands) == 0:
            logger.debug('The homology %s can not be applied, no strands in common.', [[diagram.facet(v).get_binary_vector_string() for v in e] for e in homology])
            logger.debug('Domain sizes of would-be application: %s %s %s', len(supp[(0,1)]), len(supp[(1,2)]), len(supp[(0,2)]))
            raise InapplicableHomology()
        strand_alteration = {
            (0,1) : (lambda x, y: x.difference(y)),
            (1,2) : (lambda x, y: x.difference(y)),
            (0,2) : (lambda x, y: x.union(y)),
        }
        new_supp = {
            (i,j) : list(strand_alteration[(i,j)](set(supp[(i,j)]), set(movable_strands))) for i, j in D.values()
        }
        return [supp, new_supp, s, D]

    def get_simplex_vertex(index, edges):
        """
        The information identifying vertices (the signatures) is contained with
        the edges themselves. This helper function picks out the vertex
        signature by index by traversing the edges.

        Args:
            index (int):
                The index of the desired vertex (0, 1, or 2).

            edges (dict):
                The edge signatures. Keys are (0,1), (1,2), (0,2), and each
                value is a pair of signatures.
        """
        for i, j in edges.keys():
            if i == index:
                return edges[(i,j)][0]
            if j == index:
                return edges[(i,j)][1]


class InapplicableHomology(Exception):
    pass
