#!/usr/bin/env python3

import datacube
from datacube.modeling_pipeline import ModelingPipeline, OutputFormats

points = np.array([
    [0,0,1,1,1,1],
    [1,0,1,1,1,1],
    [1,1,1,1,1,1],
    [0,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,1,1,0,0,0],
])
# points = np.loadtxt('colon_test_binary_10.csv', delimiter=',')
# m = pd.read_csv('colon_gtex_cell_cycle_genes_no_labels_50_samples_15_genes.csv', index_col=0)

pipeline = ModelingPipeline(
    points,
    output_formats=[OutputFormats.GRAPHML, OutputFormats.PNG],
    interactive=True,
)
pipeline.run()
