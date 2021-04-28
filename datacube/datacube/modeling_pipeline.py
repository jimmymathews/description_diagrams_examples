import os
from os.path import basename
import enum
from enum import Enum, auto

import numpy as np
import networkx as nx

from .io import binary_matrix_from_file
from .io import bipartite_diagram_representation_of_feature_matrix
from .steiner_reduction import SteinerReduction, RandomResolutionOfBipartite
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
        self.input = input
        self.output_formats = output_formats
        self.interactive = interactive

    def run(self):
        points = self.gather_input()
        random_resolver = RandomResolutionOfBipartite()
        reducer = SteinerReduction()
        mutation_steps = [random_resolver, reducer]
        if self.interactive:
            progress_bar = ProgressBar()
            for mutation_step in mutation_steps:
                mutation_step.add_progress_listener(progress_bar)
        diagram = bipartite_diagram_representation_of_feature_matrix(points)
        for mutation_step in mutation_steps:
            mutation_step.mutate(diagram)
        self.annotate(diagram)
        self.send_to_output(diagram)

    def gather_input(self):
        if type(self.input) == str:
            filename = self.input
            self.file_basename = basename(filename).rstrip('.csv').rstrip('.CSV')
            return binary_matrix_from_file(filename)
        else:
            self.file_basename = 'graph'
            return self.input

    def annotate(self, diagram):
        for s1, s2 in diagram.graph.edges:
            l = diagram.graph.edges[s1, s2]['supports feature inference']
            diagram.graph.edges[s1, s2]['supports feature inference'] = 0
            diagram.graph.edges[s1, s2]['number of feature inferences supported'] = len(l)
            diagram.graph.edges[s1, s2]['object set size'] = len(set([obj for obj, feature in l]))
            diagram.graph.edges[s1, s2]['feature set size'] = len(set([feature for obj, feature in l]))
        for s in diagram.graph.nodes:
            diagram.graph.nodes[s]['number of features'] = len(s)
            diagram.graph.nodes[s]['binary representation'] = diagram.facet(s).get_binary_vector_string()

    def send_to_output(self, diagram):
        for f in self.output_formats:
            if f == OutputFormats.GRAPHML:
                nx.write_graphml(diagram.graph, self.file_basename + '.graphml')
            if f == OutputFormats.PNG:
                save_preview_visualization(self.file_basename + '.png', diagram)

