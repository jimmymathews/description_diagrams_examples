#!/usr/bin/env python3
import numpy as np

import datacube
from datacube.io import binary_matrix_from_file

def test_matrix_import():
    f1 = 'fully_annotated.csv'
    f2 = 'only_column_names.csv'
    f3 = 'bare_matrix.csv'

    m1 = binary_matrix_from_file(f1)
    m2 = binary_matrix_from_file(f2)
    m3 = binary_matrix_from_file(f3)

    assert((m1 == m2).all())
    assert((m1 == m3).all())
    assert((m1 == np.matrix([[1,0,1], [0,1,0]])).all())


if __name__=='__main__':
    test_matrix_import()
