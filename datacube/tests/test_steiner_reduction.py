#!/usr/bin/env python3
import numpy as np

import datacube
from datacube.io import bipartite_diagram_representation_of_feature_matrix
from datacube.steiner_reduction import SteinerReduction, RandomResolutionOfBipartite
from datacube.progress_bar import ProgressBar
from datacube.io import save_preview_visualization

points = np.array([
    [0,0,1,1,1,1],
    [1,0,1,1,1,1],
    [1,1,1,1,1,1],
    [0,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,1,1,0,0,0],
])
# points = np.loadtxt('colon_test_binary_10.csv', delimiter=',')
# import pandas as pd
# m = pd.read_csv('colon_gtex_cell_cycle_genes_no_labels_50_samples_15_genes.csv', index_col=0)
# threshold = lambda x: 1 if x > 0.5 else 0
# threshold = np.vectorize(threshold)
# points = threshold(m.to_numpy())

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

for s1, s2 in diagram.graph.edges:
    l = diagram.graph.edges[s1, s2]['supports feature inference']
    diagram.graph.edges[s1, s2]['supports feature inference'] = 0
    diagram.graph.edges[s1, s2]['number of feature inferences supported'] = len(l)
    diagram.graph.edges[s1, s2]['object set size'] = len(set([obj for obj, feature in l]))
    diagram.graph.edges[s1, s2]['feature set size'] = len(set([feature for obj, feature in l]))

def binarize(s, dim=0):
    ss = []
    for i in range(dim):
        if i in s:
            ss.append('1')
        elif i+dim in s:
            ss.append('0')
        else:
            ss.append('-')
    return ''.join(ss)

for s in diagram.graph.nodes:
    diagram.graph.nodes[s]['number of features'] = len(s)
    diagram.graph.nodes[s]['binary representation'] = binarize(s, dim=6)

import networkx as nx
nx.write_graphml(diagram.graph, 'graph_example.graphml')


# att_list = 'supports feature inference'
# G = diagram.graph
# for n1, n2, d in G.edges(data=True):
#     for att in att_list:
#         d.pop(att, None)
