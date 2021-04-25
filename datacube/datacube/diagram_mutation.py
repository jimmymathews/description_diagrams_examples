from .data_structures import DescriptionDiagram


class DescriptionDiagramMutator:
    """
    An interface for a generic unit of alteration to an existing description
    diagram, i.e. a subgraph of the directed graph of standard hypercube
    facet containments which describes a binary feature matrix.
    """
    def mutate(self, diagram: DescriptionDiagram=None):
        pass

