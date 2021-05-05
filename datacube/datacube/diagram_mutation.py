
from .data_structures import DescriptionDiagram
from .progress_bar import ProgressingTask


class DescriptionDiagramMutator(ProgressingTask):
    """
    An interface for a generic unit of alteration to an existing description
    diagram, i.e. a subgraph of the directed graph of standard hypercube
    facet containments which describes a binary feature matrix.
    """
    def __init__(self):
        super().__init__()

    def mutate(self, diagram: DescriptionDiagram=None):
        pass
