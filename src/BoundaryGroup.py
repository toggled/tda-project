__author__ = 'Naheed'

from boundaryoperator import Boundary


class KthBoundaryGroup:
    def __init__(self, K):
        self.kth_boundary = []  # list of Objects in Bk.

    def construct_from_simplex(self, list_kplus1_simplices):
        '''
        Given a list of K+1 simplices, Compute the boundary of all those simplices and store them
        :param list_kplus1_simplices:
        :return: None
        '''
        for kplus1_simplex in list_kplus1_simplices:
            boundary = Boundary()  # deltak+1(k+1 simplex) => k simplex
            boundary.compute_boundary(kplus1_simplex)
            self.kth_boundary.append(boundary)  # the boundary list is inside delta_kplus1

    def get_kth_boundarygroup(self):
        return self.kth_boundary

    def __str__(self):
        str_repr = ''
        for boundaryobj in self.kth_boundary:
            str_repr += str(boundaryobj) + '\n'
        return str_repr
