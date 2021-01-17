#!/usr/bin/env python3

from homological_data_cube_diagrams import *
import random
import time

class FacetTraversal:
    def __init__(self, cube, initial_facet):
        self.cube = cube
        self.traversed_facet = initial_facet
        self.linked_to = {}
        self.next = None
        self.previous = None

    def get_traversed_facet(self):
        return self.traversed_facet

    def link_to(self, other):
        self.linked_to[other] = other

    def set_next(self, other):
        self.next = other

    def set_previous(self, other):
        self.previous = other

    def get_next(self):
        return self.next

    def get_previous(self):
        return self.previous

    def mutate(self):
        sig_prev = self.get_previous().get_traversed_facet().get_signature()
        sig = self.get_traversed_facet().get_signature()
        sig_next = self.get_next().get_traversed_facet().get_signature()
        removed_by_self = set(sig_prev).difference(sig)
        removed_by_next = set(sig).difference(sig_next)

        mutated_sig = tuple(sorted(set(sig).difference(removed_by_next).union(removed_by_self)))
        self.traversed_facet = self.cube.get_facet(mutated_sig, forgetful=True)

    def cover_same_ground(traversal1, traversal2):
        return traversal1.get_traversed_facet() == traversal2.get_traversed_facet()

    def fuse(traversal1, traversal2):
        traversal1.link_to(traversal2)
        traversal2.link_to(traversal1)

    def __repr__(self):
        return repr(self.get_traversed_facet())

    def get_emphasis(self):
        sig = self.get_traversed_facet().get_signature()

        if self.get_previous() != None:
            sig_prev = self.get_previous().get_traversed_facet().get_signature()
            removed_by_self = set(sig_prev).difference(sig)
        else:
            removed_by_self = []

        if self.get_next() != None:
            sig_next = self.get_next().get_traversed_facet().get_signature()
            removed_by_next = set(sig).difference(sig_next)
        else:
            removed_by_next = []

        entries = list(removed_by_next) + list(removed_by_self)
        n = self.cube.dimension
        return [entry if entry < n else entry - n for entry in entries]

class Strand:
    def __init__(self, cube):
        self.cube = cube
        self.facet_traversals = []

    def random_initialization(self, vertex, face):
        n = self.cube.dimension
        other_face_integers = list(set(vertex.get_signature()).difference(set(face.get_signature())))
        face_order = random.sample(other_face_integers, k=len(other_face_integers))
        face_order.append(face.get_signature()[0])

        t = FacetTraversal(self.cube, initial_facet=vertex)
        self.facet_traversals.append(t)

        for i in range(1,len(face_order)):
            tp = FacetTraversal(self.cube, initial_facet=self.cube.get_facet(signature=face_order[i:n], forgetful=True))
            self.facet_traversals.append(tp)
            t.set_next(tp)
            tp.set_previous(t)
            t = tp

    def random_mutation(self):
        n = self.cube.dimension
        r = random.randint(1, n-2) # Ensures that (1) vertex is unchanged, (2) face is unchanged
        self.facet_traversals[r].mutate()

    # def __repr__(self):
    #     return '\n'.join([ft.get_traversed_facet().representation(emphasis = i) for i, ft in enumerate(self.facet_traversals)])


class Homotopy1Class:
    '''
    A homotopy class of (a union of) relative paths in the barycentric subdivision of a given cube.
    Restricted to maximally factored paths at the moment (must consist of 1-cells of Q norm 1).
    In this program it is used to capture the information of a given feature matrix.
    (i.e. homotopic to the raw data chain of the feature matrix).
    '''
    def __init__(self, cube):
        '''
        For each element (s,f) of the feature matrix, a path from s to f is stored in instance data here.
        Such a path is specified by an ordering f', f'', ..., f(n-1), f of the set of facets containing s, terminating in f.
        
        '''
        self.cube = cube
        self.strands = {}

    def initialize_random_raw_data_strands(feature_matrix, cube):
        n = cube.dimension
        assert(feature_matrix.shape[1] == n)
        tuples = [tuple(row) for row in feature_matrix]
        uniques = np.unique(tuples, axis=0)
        vertices = [cube.get_vertex([row[l] for l in range(n)]) for row in uniques]
        faces = [cube.get_face((l,)) for l in range(2*n)]
        paths = Homotopy1Class(cube)
        for vertex in vertices:
            for face in faces:
                if vertex.contained_in(face):
                    strand = Strand(cube)
                    strand.random_initialization(vertex, face)
                    paths.strands[strand] = strand
        return paths

    def __repr__(self):
        n = self.cube.dimension
        return '\n'.join([' '.join([strand.facet_traversals[i].get_traversed_facet().representation(emphasis=strand.facet_traversals[i].get_emphasis()) for strand in self.strands]) for i in range(n)])


class StochasticDescriptionDiagrams(Verbose):
    def __init__(self, dimension=None, feature_matrix=np.array([[]]), verbose=True):
        self.verbose = verbose
        assert(dimension == None or feature_matrix == np.array([[]]))
        if dimension != None:
            self.cube = DataCube(dimension=dimension, verbose=verbose)
        else:
            self.cube = DataCube(dimension=feature_matrix.shape[1], verbose=verbose)

        self.cube.current_coordinate_system = 'standard Euclidean'
        if feature_matrix != np.array([[]]):
            self.set_feature_matrix(feature_matrix)
        self.start_simulation()

    def set_feature_matrix(self, feature_matrix):
        self.paths = Homotopy1Class.initialize_random_raw_data_strands(feature_matrix, self.cube)

    def start_simulation(self):
        while(True):
            for strand in self.paths.strands:
                strand.random_mutation()
            self.display()
            try:
                time.sleep(0.01)
            except KeyboardInterrupt:
                print()
                return
            self.clear_display()

    def clear_display(self):
        print(clear_line, end='')
        for i in range(self.cube.dimension):
            print(up_line + clear_line, end='')

    def display(self):
        print(repr(self.paths))

