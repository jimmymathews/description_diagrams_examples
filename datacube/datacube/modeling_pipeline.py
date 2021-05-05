import os
from os.path import basename
import enum
from enum import Enum, auto

import networkx as nx

from .io import binary_matrix_from_file
from .matrix_to_graph import bipartite_diagram_representation_of_feature_matrix
from .regularization import RandomResolutionOfBipartite
from .steiner_reduction import SteinerReduction
from .annotation import DiagramAnnotator
from .progress_bar import ProgressBar
from .io import save_preview_visualization


class OutputFormats(Enum):
    GRAPHML = auto()
    PNG = auto()


class ModelingPipeline:
    def __init__(self,
        input,
        output_formats=[OutputFormats.GRAPHML, OutputFormats.PNG],
        interactive=False
    ):
    """
    Args:
        input (str or matrix-like):
            If a string, this should be the filename of a file containing
            a binary feature matrix (possibly with column/row names).

        output_formats (datacube.modeling_pipeline.OutputFormats):
            Specifies where to send output.

        interactive (boolean):
            Indicates whether to allow printing to the terminal.
    """
        self.input = input
        self.output_formats = output_formats
        self.interactive = interactive

    def run(self):
        mutation_steps = self.get_mutation_steps()
        points = self.gather_input()
        diagram = bipartite_diagram_representation_of_feature_matrix(points)
        for step in mutation_steps:
            step.mutate(diagram)
        self.send_to_output(diagram)

    def get_mutation_steps(self):
        mutation_steps = [
            RandomResolutionOfBipartite(),
            SteinerReduction(),
            DiagramAnnotator(),
        ]
        if self.interactive:
            progress_bar = ProgressBar()
            for step in mutation_steps:
                step.add_progress_listener(listener=progress_bar)

    def gather_input(self):
        if type(self.input) == str:
            filename = self.input
            self.file_basename = basename(filename).rstrip('.csv').rstrip('.CSV')
            return binary_matrix_from_file(filename)
        else:
            self.file_basename = 'graph'
            return self.input

    def send_to_output(self, diagram):
        for f in self.output_formats:
            if f == OutputFormats.GRAPHML:
                nx.write_graphml(diagram.graph, self.file_basename + '.graphml')
            if f == OutputFormats.PNG:
                save_preview_visualization(self.file_basename + '.png', diagram)

