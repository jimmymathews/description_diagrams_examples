import itertools
import numpy as np
from colors import *

class FacetSignature:
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
    def __init__(self, dimension=3, signature_object=None):
        self.dimension = dimension
        self.signature_object = signature_object
        self.parents = {}
        self.children = {}

    def register_direct_parent(self, other):
        self.parents[other] = other
        other._register_direct_child(self)

    def _register_direct_child(self, other):
        self.children[other] = other

    def contained_in(self, other):
        return other.signature_object.implies(self.signature_object)

    def __str__(self):
        return str(self.signature_object)

    def __repr__(self):
        return repr(self.signature_object)


class FacetPair:
    def __init__(self, f1, f2):
        assert(f1.contained_in(f2))
        self.f1 = f1
        self.f2 = f2

    def __repr__(self):
        return repr(self.f1) + ' ' + yellow + '\u2192' + reset + '  ' + repr(self.f2)


class DataCube:
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
        vertices = [cube.get_vertex([row[l] for l in range(n)]) for row in uniques]
        chain = DescriptionChain(self, vertices)
        for vertex in chain.vertices:
            for face in cube.facets_of_dimension[n-1].values():
                if vertex.contained_in(face):
                    pair = FacetPair(vertex, face)
                    chain.coefficients[pair] = 1
        return chain



class DescriptionChain:
    def __init__(self, cube, vertices):
        self.cube = cube
        self.vertices = vertices
        self.coefficients = {}

    def __repr__(self):
        display_width = max([len(str(coefficient)) for coefficient in self.coefficients.values()])
        return '\n'.join([magenta + str(self.coefficients[pair]).ljust(display_width+1) + reset + repr(pair) for pair in self.coefficients])


class HomologicalMinimumCalculator:
    def __init__(self, feature_matrix, dimension=3):
        self.dimension = dimension
        self.feature_matrix = feature_matrix

    def calculate_barycentric_1_skeleton(self):
        pass

    def calculate_search_lattice(self):
        pass

    def get_minimal_chain(self):
        pass


cube = DataCube(dimension=5, verbose=False)
feature_matrix = np.array([[0,0,0,0,1], [0,0,0,1,1], [0,1,1,1,1]])
c0 = cube.get_raw_data_chain(feature_matrix)
print(c0)
