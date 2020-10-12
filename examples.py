#!/usr/bin/env python3

from homological_data_cube_diagrams import *

class Examples:
    def single_merge():
        feature_matrix = np.array(
            [[0,0,0],
             [0,0,1],
            ]
        )
        calculator = HomologicalMinimumCalculator(feature_matrix=feature_matrix)
        d = calculator.cube.get_d()
        p = calculator.calculate_projection()
        print(p.representation(emphasize_non_vertices=True))
        calculator.projected_c_raw.plot_histogram()
        return calculator

    def higher_dimension_single_merge(dimension=4):
        feature_matrix = np.array(
            [[0 for i in range(dimension-1)] + [0],
             [0 for i in range(dimension-1)] + [1],
            ]
        )
        calculator = HomologicalMinimumCalculator(feature_matrix=feature_matrix)
        d = calculator.cube.get_d()
        p = calculator.calculate_projection()
        print(p.representation(emphasize_non_vertices=True))
        calculator.projected_c_raw.plot_histogram()
        return calculator

    def more_points():
        feature_matrix = np.array(
            [[0,0,0,0],
             [0,0,0,1],
             [0,0,1,0],
             [0,0,1,1],
             [0,1,0,0],
             [1,0,0,0],
             [1,1,0,0],
            ]
        )
        calculator = HomologicalMinimumCalculator(feature_matrix=feature_matrix)
        d = calculator.cube.get_d()
        p = calculator.calculate_projection()
        print(p.representation(emphasize_non_vertices=True))
        calculator.projected_c_raw.plot_histogram()
        return calculator

    def cluster():
        feature_matrix = np.array(
            [[0,0,1,0,1],
             [0,0,1,1,0],
             [0,0,0,1,1],
            ]
        )
        calculator = HomologicalMinimumCalculator(feature_matrix=feature_matrix)
        d = calculator.cube.get_d()
        p = calculator.calculate_projection()
        print(p.representation(emphasize_non_vertices=True))
        calculator.projected_c_raw.plot_histogram()
        return calculator

if __name__=='__main__':
    Ex = Examples
    # Ex.single_merge()
    # Ex.higher_dimension_single_merge()
    Ex.more_points()
    # Ex.cluster()
