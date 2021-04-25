#!/usr/bin/env python3
import numpy as np

import datacube
from datacube.io import bipartite_diagram_representation_of_feature_matrix
from datacube.steiner_reduction import SteinerReduction, RandomResolutionOfBipartite
from datacube.progress_bar import ProgressBar
from datacube.io import save_preview_visualization

points = np.array([
    [1,0,1,1],
    [0,1,1,1],
    [0,0,1,1],
    [1,1,1,1],
    [0,0,0,0],
    [1,0,0,0],
    [0,1,0,0],
    [1,1,0,0],
])
# points = np.loadtxt('colon_test_binary_10.csv', delimiter=',')
progress_bar = ProgressBar()
diagram = bipartite_diagram_representation_of_feature_matrix(points)
L0 = diagram.compute_length()
E0 = len(diagram.graph.edges)

random_resolver = RandomResolutionOfBipartite()
reducer = SteinerReduction()

mutation_steps = [random_resolver, reducer]
for mutation_step in mutation_steps:
    mutation_step.add_progress_listener(progress_bar)

for mutation_step in mutation_steps:
    mutation_step.mutate(diagram)

print('Length per feature-point: ' + str(diagram.compute_length() / (points.shape[0] * points.shape[1])))
print('Intermediate edge number: ' + str(len(diagram.graph.edges) - E0))

save_preview_visualization('graph.png', diagram)

