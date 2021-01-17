#!/usr/bin/env python3

from homological_data_cube_diagrams import *
from stochastic_description_diagrams import *

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
        return calculator

    def stochastic_cluster():
        feature_matrix = np.array(
            [[0,0,1,0,1],
             [0,0,1,1,0],
             [0,0,0,1,1],
            ]
        )
        sim = StochasticDescriptionDiagrams(feature_matrix = feature_matrix)

    def stochastic_higher_dimension(dimension = 4):
        feature_matrix = np.array(
            [[0 for i in range(dimension-1)] + [0],
             [0 for i in range(dimension-1)] + [1],
            ]
        )
        sim = StochasticDescriptionDiagrams(feature_matrix = feature_matrix)



if __name__=='__main__':
    Ex = Examples
    # Ex.single_merge()
    # Ex.higher_dimension_single_merge()
    Ex.more_points()
    # Ex.cluster()
    # Ex.stochastic_higher_dimension(dimension=6)
    # Ex.stochastic_cluster()


