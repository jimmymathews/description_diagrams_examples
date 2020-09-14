import itertools
import numpy as np
from colors import *

class FacetSignature:
    '''
    An object to wrap the signature of a facet with behaviors.
    The signature is a sorted tuple of integers labelling the faces containing the facet.
    Faces are labelled as follows:
    - The 0-value face for the i-th feature is i
    - The 1-value face for the i-th features is i+n
    Indices start with 0.
    '''
    def __init__(self, signature, top=None):
        self.signature = signature
        self.top = top

    def get_values(self):
        if self.top == None:
            return ()
        n = self.top.dimension
        values = [None for j in range(n)]
        for i in range(n):
            if i in self.signature:
                values[i] = 0
            if i+n in self.signature:
                values[i] = 1
        return values

    def implies(self, other):
        return set(self.signature).issubset(other.signature)

    def __str__(self):
        return ''.join([str(element) if element!= None else '-' for element in self.get_values()])

    def __repr__(self):
        return cyan + str(self) + reset


class Facet:
    '''
    Represents a facet of a given cube.
    signature_object is an object of type FacetSignature.
    This object identifies which facet this is in the context of the given cube.
    The parents and children are direct only. They are not determined at constructor
    time because the facet containment data structure is not a tree, so it does
    not admit of one-shot recursive construction starting from the top element.
    '''
    def __init__(self, dimension=3, signature_object=None):
        self.dimension = dimension
        self.signature_object = signature_object
        self.parents = {}
        self.children = {}

    def register_direct_parent(self, other):
        self.parents[other] = other
        other._register_direct_child(self)

    def _register_direct_child(self, other):
        '''
        To be used only by register_direct_parent.
        This reduces potential for error where a parent relationship is recorded but not the inverse relation.
        '''
        self.children[other] = other

    def contained_in(self, other):
        return other.signature_object.implies(self.signature_object)

    def __str__(self):
        return str(self.signature_object)

    def __repr__(self):
        return repr(self.signature_object)


class FacetPair:
    '''
    Intended to be used as a basis element in the set of 1-chains of the barycentric subdivision of a given cube.
    '''
    def __init__(self, f1, f2):
        assert(f1.contained_in(f2))
        self.f1 = f1
        self.f2 = f2

    def __repr__(self):
        return repr(self.f1) + ' ' + yellow + '\u2192' + reset + '  ' + repr(self.f2)


class DataCube:
    '''
    Represents the standard cube complex structure on the cube of a given dimension.
    Currently it is memory intensive, storing objects for every facet.
    In the future, there may be a way to create them on the fly as needed.
    '''
    def __init__(self, dimension=3, verbose=True):
        self.verbose = verbose
        self.dimension = dimension
        n = self.dimension
        assert(n >= 1)
        self.facets = {}
        self.facets_of_dimension = {i:{} for i in range(n+1)}
        self.signature_objects = {}
        self.create_all_sub_facets()

    def create_all_sub_facets(self):
        n = self.dimension
        top = self.create_facet(dimension=n, signature=())
        for k in range(1, n+1):
            k_subsets_faces = itertools.combinations(range(2*n), k)
            for k_subset_faces in k_subsets_faces:
                if not self.check_valid_signature(k_subset_faces):
                    continue
                signature = k_subset_faces
                facet = self.create_facet(dimension=n-len(signature), signature=signature)
                if(self.verbose):
                    print(green + 'Creating facet with signature ' + reset + str(signature).ljust(25) + magenta + 'defined by: ' + reset + str(facet).ljust(25) + reset)
                self.register_all_direct_parents(facet, signature=signature)

    def register_all_direct_parents(self, facet, signature=None):
        for entry in signature:
            truncated = tuple(sorted(list(set(signature).difference(set([entry])))))
            assert(truncated in self.signature_objects.keys())
            facet.register_direct_parent(self.get_facet(truncated))

    def check_valid_signature(self, signature):
        n = self.dimension
        for i in range(n):
            if i in signature and i+n in signature:
                return False
        return True

    def create_facet(self, dimension=0, signature=()):
        top = None
        if signature != ():
            top = self.get_facet(())
        signature_object = FacetSignature(signature=signature, top=top)
        self.signature_objects[signature] = signature_object
        f = Facet(dimension=dimension, signature_object=signature_object)
        self.facets[signature_object] = f
        self.facets_of_dimension[dimension][signature_object] = f
        return f

    def get_facet(self, signature):
        return self.facets[self.signature_objects[signature]]

    def get_vertex(self, vector):
        n = self.dimension
        signature = tuple(sorted([i if vector[i] == 0 else i+n for i in range(len(vector))]))
        return self.get_facet(signature)

    def __getitem__(self, signature):
        if signature in self.signature_objects.keys():
            return self.get_facet(signature)
        else:
            return None

    def get_raw_data_chain(self, feature_matrix):
        n = self.dimension
        assert(feature_matrix.shape[1] == n)
        tuples = [tuple(row) for row in feature_matrix]
        uniques = np.unique(tuples, axis=0)
        vertices = [self.get_vertex([row[l] for l in range(n)]) for row in uniques]
        # chain = DescriptionChain(self, vertices)
        chain = DescriptionChain(self)
        for vertex in vertices:
            for face in self.facets_of_dimension[n-1].values():
                if vertex.contained_in(face):
                    pair = FacetPair(vertex, face)
                    chain.coefficients[pair] = 1
        return chain


class DescriptionChain:
    '''
    A 1-chain in the barycentric subdivision of a given cube, capturing the
    information of a given feature matrix (i.e. homologous to the raw data
    chain of the feature matrix).
    '''
    def __init__(self, cube):
        self.cube = cube
        self.coefficients = {}

    def get_vertices(self):
        return list(set([pair.f1 if pair.f1.dimension==0 else None for pair in self.coefficients]))

    def __repr__(self):
        display_width = max([len(str(coefficient)) for coefficient in self.coefficients.values()])
        return '\n'.join([green + str(self.coefficients[pair]).ljust(display_width+1) + reset + repr(pair) for pair in self.coefficients])


class HomologicalMinimumCalculator:
    '''
    For a given feature matrix, leading to the 1-chain c0, this calculator attempts to find the element of
    smallest geometric norm belonging to the integral homology class of c0 relative to the sample set
    (vertices of the cube) and the attribute set barycenters. Equivalently, c0 plus the closest element of
    the lattice of boundaries of relative 2-chains to the element -c0.
    Note that the latter lattice is independent of the feature matrix. Computations concerning this lattice,
    like short bases, can be performed in advance knowing only the dimension of the cube.
    '''
    def __init__(self, feature_matrix):
        self.cube = DataCube(dimension=feature_matrix.shape[1], verbose=False)
        self.c0 = self.cube.get_raw_data_chain(feature_matrix)

    def calculate_barycentric_1_skeleton(self):
        pass

    def calculate_search_lattice(self):
        pass

    def get_minimal_chain(self):
        pass


feature_matrix = np.array([[0,0,0,0,1], [0,0,0,1,1], [0,1,1,1,1]])
calculator = HomologicalMinimumCalculator(feature_matrix)
print(          'Raw data chain')
print(magenta + '--------------' + reset)
print(calculator.c0)

