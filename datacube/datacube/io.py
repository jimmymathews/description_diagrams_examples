
import networkx as nx
import matplotlib
matplotlib.use('agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plot

from .log_formats import colorized_logger
logger = colorized_logger(__name__)

def binary_matrix_from_file(filename):
    """
    Will treat first row as column names if there is at least one
    non-numeric entry in this row.
    Only in case column names are present, this function will then treat the
    first column as row names if there is at least one non-numeric entry in
    this column. (That is, column names must be present if row names are present).

    Args:
        filename (str):
            File containing a CSV-formatted matrix.

    Returns:
        matrix (numpy.ndarray):
            A binary Numpy matrix.
    """
    df = pd.read_csv(filename, header=None, keep_default_na=False)
    if any([not is_floatable(entry) for entry in df.iloc[0]]):
        row0 = list([str(entry) for entry in df.iloc[0]])
        if len(row0) == len(set(row0)):
            if row0[0] == '':
                row0[0] = 'unnamed column possibly row names'
            df.columns = row0
            df = df.iloc[[i for i in range(df.shape[0]) if i != 0]]
        else:
            logger.error(
                'First row contains non-numeric data, but values are not column names because of non-uniqueness.'
            )
            return None
    colname0 = df.columns[0]
    if any([not is_floatable(entry) for entry in df[colname0]]):
        col0 = list([str(entry) for entry in df[colname0]])
        if len(col0) == len(set(col0)):
            df.index = col0
            df = df[[df.columns[i] for i in range(df.shape[1]) if i != 0]]
        else:
            logger.error(
                'First column contains non-numeric data, but values are not row names because of non-uniqueness.'
            )
            return None
    threshold = np.vectorize(lambda x: 1 if x >= 0.5 else 0)
    matrix = threshold(df.astype(float).to_numpy())
    return matrix

def is_floatable(value):
    try:
        f_value = float(value)
        return True
    except ValueError:
        return False

def color_by_signature(signature, dimension):
    if len(signature) == dimension:
        return 'green'
    if len(signature) == 1:
        return 'blue'
    return 'red'

def save_preview_visualization(filename, diagram):
    graph = diagram.graph
    signatures = [node for node in graph]
    color_map = [color_by_signature(signature, diagram.dimension) for signature in signatures]
    figure = plot.figure(dpi=200)
    graph_prettier = nx.relabel.relabel_nodes(
        graph,
        {signature: diagram.facet(signature).get_binary_vector_string() for signature in signatures},
    )
    nx.draw(
        graph_prettier,
        node_color=color_map,
        ax=figure.add_subplot(111),
        with_labels=True,
        arrowsize=7,
        node_size=125,
        font_size=5,
        pos=nx.fruchterman_reingold_layout(graph_prettier, scale=2, k=0.5)
    )
    figure.savefig(filename)

