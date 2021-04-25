#!/usr/bin/env python3
import itertools
import numpy as np
from scipy import sparse
from scipy import linalg
import matplotlib.pyplot as plt
from math import sqrt
import os
from terminal_codes import *


class Facet:
    '''
    Represents a facet of a given cube.
    '''
    def __init__(self, signature_object=None):
        '''
        'signature_object' is an object of type FacetSignature.
        It identifies which facet this is in the context of the given cube.
        The parents and children are direct only. They are not determined at constructor time because the facet
        containment data structure is not a tree, so it does not admit of one-shot recursive construction starting
        from the top element (that would be the usual design pattern for object composition).
        '''
        self.signature_object = signature_object
        self.parents = {}
        self.children = {}

    def get_dimension(self):
        return self.signature_object.top_dimension - len(self.signature_object.signature)

    def is_top(self):
        return self.signature_object.signature == ()

    def contained_in(self, other):
        return other.signature_object.implies(self.signature_object)

    def get_cube(self):
        return self.signature_object.cube

    def get_direct_parents(self):
        '''
        Constructs the list of other Facet objects representing the direct parents of this one.
        It works by going to the FacetSignature level, i.e. the list of hyperplanes containing
        this Facet. It removes each hyperplane one by one, and records the resulting Facet of dimension
        1 higher than this one.
        '''
        signature = self.signature_object.signature
        parents = []
        for entry in signature:
            truncated = tuple(sorted(list(set(signature).difference(set([entry])))))
            parents.append(self.get_cube().get_facet(truncated))
        return parents

    def get_all_parents(self, exclude_self=False):
        '''
        The more-than-direct version of get_direct_parents. Uses get_direct_parents recursively.
        '''
        if self.is_top():
            return []
        if exclude_self:
            return list(set(itertools.chain(*[direct_parent.get_all_parents() for direct_parent in self.get_direct_parents()])))
        else:
            return [self] + list(itertools.chain(*[direct_parent.get_all_parents() for direct_parent in self.get_direct_parents()]))

    def get_direct_children(self):
        '''
        Gets the direct children of the given Facet (i.e. Facet objects of one dimension less), using
        the opposite strategy as get_direct_parents.
        '''
        signature = self.signature_object.signature
        parents = []
        n = self.signature_object.top_dimension
        for extra in setdiff(range(2*n),signature):
            augmented = list(sorted(signature + [extra]))
            if not FacetSignature.check_valid_signature(augmented):
                continue
            parents.append(self.get_cube().get_facet(augmented))
        return parents

    def get_all_children(self, exclude_self=False):
        if self.is_vertex():
            return []
        if exclude_self:
            return list(set(itertools.chain(*[direct_parent.get_all_children() for direct_parent in self.get_direct_children()])))
        else:
            return [self] + list(itertools.chain(*[direct_parent.get_all_children() for direct_parent in self.get_direct_children()]))

    def is_vertex(self):
        '''
        Returns True or False depending on whether this Facet object represents one of the vertices of
        the given cube, or instead one of the higher-dimensional facets.
        '''
        return self.signature_object.is_vertex()

    def get_hashable(self):
        return self.signature_object.get_hashable()

    def __str__(self):
        return str(self.signature_object)

    def representation(self, emphasis = None):
        if emphasis == None:
            return repr(self.signature_object)
        else:
            no_char = ' '
            vals = self.signature_object.get_values()
            svals = [str(element) if element!=None else no_char for element in vals]
            tvals = [element if not i in emphasis else emph + element + backgroundhighlight for i, element in enumerate(svals)]
            ss = ''.join([str(element) if element!= None else no_char for i, element in enumerate(tvals)])
            return backgroundhighlight + ''.join(tvals) + resetcode

    def __repr__(self):
        return self.representation()

    def __eq__(self, other):
        return self.get_hashable() == other.get_hashable()

    def __hash__(self):
        return self.get_hashable().__hash__()

    def get_signature(self):
        return self.signature_object.signature


class FacetSignature:
    '''
    An object to wrap the signature of a facet with behaviors.
    '''
    def __init__(self, signature, cube):
        '''
        The 'signature' is a sorted tuple of integers labelling the faces containing the facet.
        Faces are labelled as follows:
        - The 0-value face for the i-th feature is i
        - The 1-value face for the i-th features is i+n
        Indices start with 0.
        'top_dimension' is the dimension of the larger cube, containing the facet described by
        the signature. This dimension is needed to interpret the signature.
        '''
        self.signature = signature
        self.cube = cube
        self.top_dimension = self.cube.dimension

    def get_dimension(self):
        return self.top_dimension - len(self.signature)

    def is_vertex(self):
        return self.get_dimension() == 0

    def get_values(self):
        '''
        A 'trinary' vector representing the facet. An entry is 0 or 1 if that coordinate of a vertex is determined
        by membership in the facet, and None if not.
        '''
        n = self.top_dimension
        values = [None for j in range(n)]
        for i in range(n):
            if i in self.signature:
                values[i] = 0
            if i+n in self.signature:
                values[i] = 1
        return values

    def implies(self, other):
        '''
        The 'signature' design as a tuple of integers has the benefit that it supports this straightforward
        method of determining the containment relation.
        '''
        return set(self.signature).issubset(other.signature)

    def get_hashable(self):
        return self.signature

    def __str__(self):
        no_char = ' '
        return ''.join([str(element) if element!= None else no_char for element in self.get_values()])

    def __repr__(self):
        return backgroundhighlight + str(self) + resetcode

    def check_valid_signature(signature, n):
        '''
        Used during comprehensive construction of all facets to select 'signature's that represent non-empty
        intersections of faces.
        '''
        for i in range(n):
            if i in signature and i+n in signature:
                return False
        return True


class FacetPair:
    '''
    Intended to be used as a basis element in the set of 1-chains of the barycentric
    subdivision of a given polygonal cell complex. (In the current context, an n-cube).
    '''
    def __init__(self, f1, f2):
        '''
        The 1-cells of the barycentric subdivision of a polygonal cell complex are labelled
        by pairs (f1,f2) of cells with f1 < f2 ('<' means 'contained in' here).
        Provide Facets f1 and f2.
        '''
        assert(f1.contained_in(f2))
        self.f1 = f1
        self.f2 = f2

    def get_hashable(self):
        return (self.f1.get_hashable(), self.f2.get_hashable())

    def __repr__(self):
        return repr(self.f1) + '' + yellow + arrow + resetcode + '' + repr(self.f2)


class FacetTriple:
    '''
    Intended to be used as a basis element in the set of 2-chains of the barycentric subdivision
    of a given polygonal cell complex. (In the current context, an n-cube).
    '''
    def __init__(self, f1, f2, f3):
        '''
        The 2-dimensional simplices of the barycentric subdivision of a polygonal cell complex are labelled by a
        triple (f1,f2, f3) of polygonal cells with f1 < f2 < f3 ('<' means 'contained in' here).
        Here the usage is for cubical facets. Provide Facets f1, f2, f3.
        '''
        assert(f1.contained_in(f2))
        assert(f2.contained_in(f3))
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3

    def get_facets(self):
        return [self.f1, self.f2, self.f3]

    def get_hashable(self):
        return (self.f1.get_hashable(), self.f2.get_hashable(), self.f3.get_hashable())

    def __repr__(self):
        return repr(self.f1) + '' + red + arrow + resetcode + '' + repr(self.f2) + '' + red + arrow + resetcode + '' + repr(self.f3)


class Verbose:
    def print(self, message, end='\n'):
        if self.verbose:
            print(message, end=end)


class DataCube(Verbose):
    '''
    Represents the standard cube complex structure on the cube of a given dimension.
    The full facet lattice need not be computed in order to consider a DataCube object valid.
    For high dimensions, this would be a computationally intensive task.
    FacetPair and FacetTriple objects (representing 1- and 2-simplices in the barycentric
    subdivision of the cube) are created as needed by 'get_facet_pair' and 'get_facet_triple'
    using a factory pattern.
    '''
    def __init__(self, dimension=3, verbose=False):
        self.dimension = dimension
        n = self.dimension
        self.verbose = verbose
        assert(n >= 2)
        self.signature_objects = {}
        self.facets = {}
        self.facets_of_dimension = {i:{} for i in range(n+1)}
        self.pairs = {}
        self.triples = {}
        self.C1_basis = {}
        self.C2_basis = {}
        self.filled_out = False

    def create_all_sub_facets(self):
        '''
        Creates objects for each facet, and registers all direct containments.
        Possibly very computationally intensive. Use only if necessary.
        '''
        if self.filled_out:
            return
        else:
            self.filled_out = True
        self.print('Creating ' + yellow + 'data cube' + resetcode + ' facet objects.')
        n = self.dimension
        top = self.create_facet(dimension=n, signature=())
        for k in range(1, n+1):
            k_subsets_faces = itertools.combinations(range(2*n), k)
            for k_subset_faces in k_subsets_faces:
                if not FacetSignature.check_valid_signature(k_subset_faces, n):
                    continue
                signature = k_subset_faces
                facet = self.create_facet(dimension=n-len(signature), signature=signature)
                self.print(green + 'Creating facet with signature ' + resetcode + str(signature).ljust(25) + magenta + 'defined by: ' + resetcode + str(facet).ljust(25) + resetcode)

    def create_facet(self, dimension=0, signature=()):
        top = None
        n = self.dimension
        if signature != ():
            top = self.get_facet(())
        signature_object = FacetSignature(signature=signature, cube=self)
        self.signature_objects[signature] = signature_object
        f = Facet(signature_object=signature_object)
        self.facets[signature_object] = f
        self.facets_of_dimension[dimension][signature_object] = f
        return f

    def get_facet(self, signature, forgetful=False):
        '''
        Factory pattern for creating a facet of the cube on the fly.
        '''
        if not forgetful:
            if not signature in self.signature_objects:
                self.create_facet(dimension=self.dimension-len(signature), signature=signature)
            return self.facets[self.signature_objects[signature]]            
        else:
            return Facet(signature_object=FacetSignature(signature=signature, cube=self))

    def get_vertex(self, vector):
        '''
        Retrieve Facet object representing a vertex labelled by a binary n-vector.
        '''
        n = self.dimension
        signature = tuple(sorted([i if vector[i] == 0 else i+n for i in range(len(vector))]))
        return self.get_facet(signature)

    def get_face(self, signature):
        '''
        Retrieve Facet object representing a face labelled by a 1-length signature.        
        '''
        assert(len(signature)==1)
        return self.get_facet(signature)

    def get_facet_pair(self, f1, f2):
        '''
        Factory pattern for creating a 1-cell in the barycentric subdivision on the fly.
        '''
        pair = FacetPair(f1, f2)
        key = pair.get_hashable()
        if key in self.pairs:
            return self.pairs[key]
        else:
            self.pairs[key] = pair
            return pair

    def get_facet_triple(self, f1, f2, f3):
        '''
        Factory pattern for creating a 2-simplex in the barycentric subdivision on the fly.
        '''
        triple = FacetTriple(f1, f2, f3)
        key = triple.get_hashable()
        if key in self.triples:
            return self.triples[key]
        else:
            self.triples[key] = triple
            return triple

    def calculate_C1_basis(self):
        self.create_all_sub_facets()
        self.C1_basis_indices = {}
        self.C1_basis = {}
        n = self.dimension
        count = 0
        for facet in self.facets.values():
            f1 = facet
            parents = f1.get_all_parents(exclude_self=True)
            for parent in parents:
                f2 = parent
                pair = self.get_facet_pair(f1,f2)
                self.C1_basis_indices[pair] = count
                self.C1_basis[count] = pair
                count = count + 1

    def calculate_C2_basis(self):
        self.create_all_sub_facets()
        self.C2_basis_indices = {}
        self.C2_basis = {}
        n = self.dimension
        count = 0
        for facet in self.facets.values():
            f1 = facet
            parents1 = f1.get_all_parents(exclude_self=True)
            for parent1 in parents1:
                f2 = parent1
                parents2 = f2.get_all_parents(exclude_self=True)
                for parent2 in parents2:
                    f3 = parent2
                    triple = self.get_facet_triple(f1, f2, f3)
                    self.C2_basis_indices[triple] = count
                    self.C2_basis[count] = triple
                    count = count + 1

    def get_C1_basis(self):
        if self.C1_basis == {}:
            self.calculate_C1_basis()
        return self.C1_basis

    def get_C2_basis(self):
        if self.C2_basis == {}:
            self.calculate_C2_basis()
        return self.C2_basis

    def get_C1_basis_indices(self):
        if self.C1_basis == {}:
            self.calculate_C1_basis()
        return self.C1_basis_indices

    def get_C2_basis_indices(self):
        if self.C2_basis == {}:
            self.calculate_C2_basis()
        return self.C2_basis_indices

    def display_C1_basis(self):
        for basis_element in self.get_C1_basis().values():
            self.print(basis_element)

    def display_C2_basis(self):
        for basis_element in self.get_C2_basis().values():
            self.print(basis_element)

    def get_d(self):
        if not hasattr(self, 'd'):
            self.calculate_boundary_operator()
        return self.d

    def calculate_boundary_operator(self):
        self.print('')
        self.print(yellow + 'B' + resetcode + ' = barycentric simplicial subdivision of boundary of standard ' + str(self.dimension) + '-cube')
        self.print('Calculating ' + yellow + 'C\u2081B' + resetcode + ' (degeneracies excluded) ... ')
        self.print(green + str(len(self.get_C1_basis())) + resetcode + ' elements')
        self.print('Calculating ' + yellow + 'C\u2082B' + resetcode + ' (degeneracies excluded) ... ')
        self.print(green + str(len(self.get_C2_basis())) + resetcode + ' elements')
        N = len(self.get_C1_basis())
        K = len(self.get_C2_basis())
        c2_index = []
        c1_index = []
        values = []
        for i in range(K):
            c2 = self.get_C2_basis()[i]
            f1, f2, f3 = c2.get_facets()

            f1f2 = self.get_C1_basis_indices()[self.get_facet_pair(f1,f2)]
            f2f3 = self.get_C1_basis_indices()[self.get_facet_pair(f2,f3)]
            f1f3 = self.get_C1_basis_indices()[self.get_facet_pair(f1,f3)]

            c2_index = c2_index + [   i,    i,    i]
            c1_index = c1_index + [f1f2, f2f3, f1f3]
            values =   values   + [   1,    1,   -1]

        self.print('Calculating ' + yellow + 'boundary operator ' + resetcode + magenta + roundd + yellow + ': C\u2082B ' + arrow + ' C\u2081B' + resetcode + ' ...')
        self.d_sparse = sparse.coo_matrix((np.array(values), (np.array(c1_index), np.array(c2_index))), dtype=float)
        self.d = self.d_sparse.toarray()

        for i in self.get_C2_basis():
            basis_element = self.get_C2_basis()[i]
            c1= DescriptionChain.from_vector(self.d[:,i], self)
            self.print('')
            self.print(magenta + roundd + resetcode + '(' + repr(basis_element) + ') =')
            self.print(repr(c1))
            if (i > 5):
                self.print('...')
                self.print('')
                break


class DescriptionChain:
    '''
    A 1-chain in the barycentric subdivision of a given cube.
    In this program it is used to capture the information of a given feature matrix
    (i.e. homologous to the raw data chain of the feature matrix).
    '''
    def __init__(self, cube):
        '''
        '''
        self.cube = cube
        self.coefficients = {}

    def to_vector(self):
        '''
        Create vector representation of the 1-chain, with respect to ordered basis.
        '''
        basis = self.cube.get_C1_basis()
        pairs = basis.values()
        v = np.zeros(len(pairs), dtype=float)
        for i, pair in enumerate(pairs):
            if pair in self.coefficients.keys():
                v[i] = self.coefficients[pair]
        return v

    def from_vector(v, cube):
        basis = cube.get_C1_basis()
        c = DescriptionChain(cube)
        pairs = basis.values()
        for i, pair in enumerate(pairs):
            if v[i] != 0:
                c.coefficients[pair] = v[i]
        return c

    def geometric_norm(self, only_support=False, p=2):
        l = 0
        for pair in self.coefficients:
            coefficient = self.coefficients[pair]
            if only_support:
                l = l + DescriptionChain.get_basis_norm(pair)
            else:
                if p == 2:
                    l = l + coefficient*coefficient * DescriptionChain.get_basis_norm(pair)
                elif p == 1:
                    l = l + abs(coefficient) * DescriptionChain.get_basis_norm(pair)
        if p == 2:
            return sqrt(l)
        elif p == 1:
            return l

    def euclidean_norm(self, only_support=False):
        l = 0
        for pair in self.coefficients:
            coefficient = self.coefficients[pair]
            if only_support:
                l = l + 1
            else:
                l = l + coefficient*coefficient
        return sqrt(l)

    def get_basis_norm(pair):
        return (pair.f2.get_dimension() - pair.f1.get_dimension())

    def color_line(self, pair):
        return green if pair.f1.get_dimension() > 0 else magenta

    def pretty_num(self, c):
        if c == int(c):
            return int(c)
        else:
            return c

    def representation(self, emphasize_non_vertices=False):
        display_width = max([len(str(self.pretty_num(coefficient))) for coefficient in self.coefficients.values()])
        kv = [item for item in self.coefficients.items()]
        kv0 = [item for item in self.coefficients.items() if item[0].f1.get_dimension() > 0]
        kvs = sorted(kv, key=lambda x: -1*abs(x[1]))
        kvs0 = sorted(kv0, key=lambda x: -1*abs(x[1]))
        # if not omit_vertices:
        lines = [self.color_line(pair) + str(self.pretty_num(coefficient)).rjust(display_width) + ' ' + resetcode + repr(pair) for pair, coefficient in kvs]
        # else:
            # lines = [self.color_line(pair) + str(coefficient).rjust(display_width) + ' ' + resetcode + repr(pair) for pair, coefficient in kvs0]
        # truncation = min(len(lines), 25)
        # truncation = len(lines)
        # if truncation != len(lines):
        #     lines = lines[0:truncation] + ['...', '\n']
        if emphasize_non_vertices:
            top10 = kvs0[0:(min(10, len(kvs0)))]
            top10 = [self.color_line(pair) + str(coefficient).rjust(display_width) + ' ' + resetcode + repr(pair) for pair, coefficient in top10]
            lines = lines + ['\n', '...', '\n', 'Non-trivial top 10 support 1-simplices in ' + yellow + 'Q' + resetcode + '-optimal representative of 1-homology class of raw-data 1-chain', '(original-vertex-sourced 1-simplices excluded; facet-center-sourced 1-simplices included)'] + top10
        return '\n'.join(lines)

    def __repr__(self):
        return self.representation()

    def plot_histogram(self):
        fig, axs = plt.subplots(1, 2)
        fig.suptitle('distribution of coefficients of 1-chain (resp. excl. vertex-sourced basis elements)')
        axs[0].hist(list(self.coefficients.values()), color = 'blue', edgecolor = 'black', bins = int(50))
        axs[1].hist([coefficient for pair, coefficient in self.coefficients.items() if pair.f1.get_dimension()>0], color = 'green', edgecolor = 'black', bins = int(50))
        plt.show()

    def get_raw_data_chain(feature_matrix, cube):
        '''
        'feature_matrix' a binary matrix.
        Returns the 1-chain in the barycentric subdivision of the cube which is the
        sum of one 1-cell for each entry of the feature matrix; for each vertex labelled by
        a column of 'feature_matrix', an edge connecting that vertex to the face center for each
        face to which it belongs.
        '''
        n = cube.dimension
        assert(feature_matrix.shape[1] == n)
        tuples = [tuple(row) for row in feature_matrix]
        uniques = np.unique(tuples, axis=0)
        vertices = [cube.get_vertex([row[l] for l in range(n)]) for row in uniques]
        faces = [cube.get_face((l,)) for l in range(2*n)]
        chain = DescriptionChain(cube)
        for vertex in vertices:
            for face in faces:
                if vertex.contained_in(face):
                    pair = cube.get_facet_pair(vertex, face)
                    chain.coefficients[pair] = 1
        return chain


class GeometricWeightedCoordinateChanger(object):
    '''
    A setup-and-tear-down-able class to allow conceptually straightforward
    performance of operations on the calculator's geometric items in a coordinate
    system for which the Q norm becomes the Euclidean norm.
    '''
    def __init__(self, calculator):
        self.calculator = calculator

    def __enter__(self):
        self.calculator.change_C1_coordinates('geometric weighted')

    def __exit__(self, exception_type, exception_value, traceback):
        self.calculator.change_C1_coordinates('standard Euclidean')


class HomologicalMinimumCalculator(Verbose):
    '''
    For a given feature matrix, leading to the 1-chain c_raw, this calculator attempts to find the element of
    smallest geometric norm (Q) belonging to the integral homology class of c_raw relative to the union of the
    sample set (certain vertices of the cube) and the attribute set barycenters.
    Equivalently, c_raw plus: the closest element of the lattice of boundaries of relative 2-chains to the element -c_raw.

    Note that this lattice is independent of the feature matrix. Computations concerning this lattice,
    like short bases, can be performed in advance knowing only the dimension of the cube.

    Note also that there is a shortcut heuristic using basic linear algebra, namely to find the intersection
    of the translation of the plane P spanned by the lattice to c_raw and the Q-orthogonal complement of subspace P. The
    resulting point will perhaps be close to the sought-after lattice element, but will not typically have integer coefficients.

    To do this last operation, one still needs a basis for the boundaries subspace.
    This can be done as follows. The matrix of the mapping C2(B)->C1(B) consists of rows of the following format.
    A basis element of C2(B) is given by an ordered triple of facets, one properly contained in the next:
        f1 < f2 < f3.
    Its boundary is the element of C1(B) with coefficients:
        1 for facet pair (f1, f3), -1 for facet pairs (f1, f2) and (f2, f3), 0 else
    These coefficients lined up with respect to an ordering of the basis of C1(B) comprises a row vector.

    After forming the matrix consisting of the rows as above, one needs a basis for the column space.
    ...
    I tried fpylll, but it was not able to do the CVP as far as I can tell. It seems not to support floating point lattice coordinates or other-vector coordinates.
    ...

    The dimension of C1(B) is equal to the number of pairs (f1, f2), f1 < f2 (excluding cases f1=top and f2=top).
    This dimension can also be calculated as follows:
    Each basis element of C1(B) labels a 1-cell that belongs to a unique open facet of the original standard cube complex.
    The number of such 1-cells in a given original facet is equal to the number of vertices of that facet, i.e. 2^dim(facet).
    Thus the number of 1-cells of the barycentric subdivision is equal to:
        sum over m: 2^m * #{ facets of standard cube which have dimension m }
    According to wikipedia, that # is:
        2^(n-m) (n choose m)
    Thus the formula for the number of basis elements of C1(B) is:
        sum over m: 2^n * (n choose m)
          = 2^n * sum over m: (n choose m)
          = 2^n * (2^n - 1)
          = 2^(2n) - 2^n
    (some error here)
    '''
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
        self.calculate_minimal_chain() # Unimplemented

    def set_feature_matrix(self, feature_matrix):
        self.print('Calculating ' + yellow + 'raw data 1-chain' + resetcode + ' (bipartite graph) ...')
        self.c_raw = DescriptionChain.get_raw_data_chain(feature_matrix, self.cube)
        self.print(self.c_raw)
        self.print('')

    def calculate_projection(self):
        '''
        Calculates Q-orthogonal projection of the raw data 1-chain onto the Q-orthogonal
        complement of the real space of the 1-boundaries subspace.
        '''
        with GeometricWeightedCoordinateChanger(self):
            self.print(yellow + 'Q' + resetcode + ' = quadratic form on C\u2081B induced by geometric length of 1-cells')
            self.print('')
            self.print('Calculating ' + yellow + 'projector' + resetcode + ' ...')
            ob = linalg.orth(self.cube.get_d())
            m = np.matmul(ob, ob.transpose())
            self.projector = np.identity(ob.shape[0]) - np.matmul(ob, ob.transpose())
            self.print(self.projector)
            self.print('')
            self.print('Calculating ' + yellow + 'Q-orthogonal projection' + resetcode + ' of raw data 1-chain into Q-complement of ' + yellow + 'Image(' + roundd + ')' + resetcode + ' ...')
            v = self.c_raw.to_vector()
            w = np.matmul(self.projector, v)
            self.projected_c_raw = DescriptionChain.from_vector(w, self.cube)
            self.print(self.projected_c_raw)
        return self.projected_c_raw

    def change_C1_coordinates(self, system):
        '''
        The geometric norm Q on C1(B) is conformal to the standard Euclidean norm with respect to the
        standard basis, by means of scalars for each basis element.
        This function converts the expressions of self.d (the boundary operator C2(B)->C1(B)) and cr
        into the corresponding expressions in Q-orthonormal coordinates, and back.
        It is supposed to be used to allow the application of a standardized closest-vector algorithm.
        '''
        if self.cube.current_coordinate_system == system:
            print(red + 'Warning' + resetcode + ': Trying to convert back to coordinates - ' + system + '?')
            return
        if system == 'standard Euclidean' and self.cube.current_coordinate_system == 'geometric weighted':
            self.cube.current_coordinate_system = system
            for pair, coefficient in self.c_raw.coefficients.items():
                norm = sqrt(DescriptionChain.get_basis_norm(pair))
                self.c_raw.coefficients[pair] = coefficient / norm
            for pair, coefficient in self.c_raw.coefficients.items():
                norm = sqrt(DescriptionChain.get_basis_norm(pair))
                if hasattr(self, 'projected_c_raw'):
                    self.projected_c_raw.coefficients[pair] = coefficient / norm
            for index, pair in self.cube.get_C1_basis().items():
                norm = sqrt(DescriptionChain.get_basis_norm(pair))
                v = np.multiply(self.cube.d[index,:], 1/norm)
                self.cube.get_d()[index,:] = v
            return
        if system == 'geometric weighted' and self.cube.current_coordinate_system == 'standard Euclidean':
            self.cube.current_coordinate_system = system
            for pair, coefficient in self.c_raw.coefficients.items():
                norm = sqrt(DescriptionChain.get_basis_norm(pair))
                self.c_raw.coefficients[pair] = coefficient * norm
            for pair, coefficient in self.c_raw.coefficients.items():
                norm = sqrt(DescriptionChain.get_basis_norm(pair))
                if hasattr(self, 'projected_c_raw'):
                    self.projected_c_raw.coefficients[pair] = coefficient * norm
            for index, pair in self.cube.get_C1_basis().items():
                norm = sqrt(DescriptionChain.get_basis_norm(pair))
                v = np.multiply(self.cube.get_d()[index,:], norm)
                self.cube.get_d()[index,:] = v
            return

    def calculate_minimal_chain(self):
        '''
        Algorithm:
        1. Start from raw data relative 1-chain cr.
        2. Consider simple homologies, the boundaries of all possible 2-simplices which pertain to the support of cr.
        3. Apply the one which gives the best geometric norm improvement.
        4. Repeat using the current 1-chain until there is no improvement to be had with this method.

        Implementation tasks:
        1. Iterator over simple homologies pertaining to a given 1-chain. Should give a 1-chain (2-simplex itself not necessary, I guess).
        2. Application of a homology to a 1-chain; add to DescriptionChain perhaps. It's just adding, but adding is non-trivial in current sparse representation. Drop 0-coefficiented elements along the way, for example?
        3. Evaluate all and pick loop.
        '''
        pass

