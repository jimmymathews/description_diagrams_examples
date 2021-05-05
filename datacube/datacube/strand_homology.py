import enum
from enum import Enum, auto


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
            logger.error('This homology can not be applied, no strands in common.')
        strand_alteration = {
            (0,1) : (lambda x, y: x.difference(y)),
            (1,2) : (lambda x, y: x.difference(y)),
            (0,2) : (lambda x, y: x.union(y)),
        }
        new_supp = {
            (i,j) : list(strand_alteration[(i,j)](set(supp[(i,j)]), set(movable_strands))) for i, j in D.values()
        }
        return [supp, new_supp, s, D]

    def get_simplex_vertex(index, e):
        for i, j in e.keys():
            if i == index:
                return e[(i,j)][0]
            if j == index:
                return e[(i,j)][1]
