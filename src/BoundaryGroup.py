__author__ = 'Naheed'

from boundaryoperator import Boundary
import numpy as np


class KthBoundaryGroup:
    def __init__(self, K):
        self.kth_boundary = []  # list of Objects in Bk.
        self.kplus1_simplices = []
        self.ksimplices = []

    def construct_from_simplex(self, list_kplus1_simplices):
        '''
        Given a list of K+1 simplices, Compute the boundary of all those simplices and store them
        :param list_kplus1_simplices:
        :return: None
        '''
        self.kplus1_simplices = list_kplus1_simplices
        self.lookuptable = {}  # lookup table for storing K dimensional id,simplex mapping
        self.column_lookuptable = {}  # lookup table for storing K+1 dimensional id , simplex mapping
        self.uniqueid = 0

        for id, kplus1_simplex in enumerate(list_kplus1_simplices):
            boundary = Boundary()  # deltak+1(k+1 simplex) => k simplex
            boundary.compute_boundary(kplus1_simplex)
            self.column_lookuptable[tuple(kplus1_simplex.kvertices)] = id
            for sign, simplex in boundary.boundary:
                if self.lookuptable.get(tuple(simplex.kvertices), -1) == -1:
                    self.lookuptable[tuple(simplex.kvertices)] = self.uniqueid
                    # print tuple(simplex.kvertices)
                    self.uniqueid += 1
            self.kth_boundary.append(boundary)  # the boundary list is inside delta_kplus1
            kplus1_simplex.id = id

    def get_kth_boundarygroup(self):
        return self.kth_boundary

    def __str__(self):
        str_repr = ''
        for boundaryobj in self.kth_boundary:
            str_repr += str(boundaryobj) + '\n'
        return str_repr

    def build_transformation_matrix(self):
        '''
        Builds the transformation matrix from k+1 simplex to k-simplex.
        Objects in the row: sorted(self.lookuptable.items(),key = lambda x: x[1])
        Objects in the column: self.kplus1_simplices
        :return:
        '''
        # self.rows = np.array(len(self.kplus1_simplices))
        # self.cols = np.array(len(self.ksimplices))
        self.transformation_matrix = np.zeros((self.uniqueid, len(self.kplus1_simplices)), dtype=np.int32)

        # print self.lookuptable
        # print self.column_lookuptable
        for col, simplices in enumerate(self.kplus1_simplices):
            # print col, ' ', str(simplices)
            signs = np.zeros(self.uniqueid, dtype=np.int32)

            id_simplex_in_boundary = []
            for sign, simplex_in_boundary in self.kth_boundary[col].boundary:
                # print sign, ' ', str(simplex_in_boundary), ' ', self.lookuptable[tuple(simplex_in_boundary.kvertices)]
                id_of_simplex = self.lookuptable[tuple(simplex_in_boundary.kvertices)]
                id_simplex_in_boundary.append(id_of_simplex)
                signs[id_of_simplex] = sign

            # print signs

            self.transformation_matrix[:, col] = signs
        # print self.transformation_matrix

        return self.transformation_matrix

    def print_rowobjects(self):
        print [object for object, idx in
               sorted(self.lookuptable.items(), key=lambda x: x[1])]  # sorted list of objects along the row

    def print_columnobjects(self):
        print [str(simplex) for simplex in self.kplus1_simplices]