## Tools for attempted solutions of the data-cube Steiner problem; homological and stochastic variants

### Examples

```
python3 examples.py
```

Out of the box this will run the simplest example, "single merge". To run the others, uncomment the corresponding lines in `examples.py`:

```

Ex.single_merge()
# Ex.higher_dimension_single_merge()
# Ex.more_points()
# Ex.cluster()
# Ex.stochastic_higher_dimension(dimension=6)
# Ex.stochastic_cluster()
...
```

What's going on is:

1. An initial set of vertices of a cube of a given dimension is chosen, for the purpose of the example ('single merge' is 2 adjacent vertices on a standard, 3-dimensional cube; 'more points' is two groups of 4 vertices, each group the set of 4 vertices on a square in a 4-cube.).
2. The "real-coefficients homological" approach is used to determine a heuristically-optimal graph-ical description of the vertex set.
3. In this approach, the initial vertex set determines a 1-chain `c_raw` in the barycentric subdivision of the cube.
4. The geometric length on the cube surface defines a norm on the 1-chains space.
5. We calculate the member of the homology class of `c_raw` of minimal norm (like a harmonic form).

Details can be found in `linear_quadratic_optimization_homology_class.pdf`. Background on the Formal Concept Analysis viewpoint on the problem is helpful, and some of this background can be gleaned from `presentation_fca.pdf` (for an overview) and `fca_notes_for_cubical_steiner.pdf` (for some more details; very rough, only read the parts that make sense).

Two major arms of this research project remain to be completed:
1. Implement the straightforward version of the homologically-minimal description involving positive integer coefficients. See the open github issue "Enumerating atomic homology cases". As I noted there, this is the low-hanging fruit!
2. Investigate a possible formulation of optimal concept binding via matroids. Nissim and I began this but did not go very far. This would almost certainly lead to a purely combinatorial and efficiently computable notion of minimality. This is a very exciting direction!



