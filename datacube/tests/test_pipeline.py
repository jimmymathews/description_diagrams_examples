#!/usr/bin/env python3

import numpy as np

import datacube
from datacube.modeling_pipeline import ModelingPipeline, OutputFormats

def test_direct_input():
    points = np.array([
        [0,0,1,1,1,1],
        [1,0,1,1,1,1],
        [1,1,1,1,1,1],
        [0,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,1,1,0,0,0],
    ])
    pipeline = ModelingPipeline(
        points,
        output_formats=[OutputFormats.GRAPHML, OutputFormats.PNG],
        interactive=True,
    )
    pipeline.run()

def test_from_file():
    # input_filename = 'colon_gtex_cell_cycle_genes_no_labels_50_samples_15_genes.csv'
    input_filename = 'colon_gtex_cell_cycle_genes_no_labels_20_samples.csv'
    pipeline = ModelingPipeline(
        input_filename,
        output_formats=[OutputFormats.GRAPHML, OutputFormats.PNG],
        interactive=True,
    )
    pipeline.run()


if __name__=='__main__':
    test_direct_input()
    test_from_file()
