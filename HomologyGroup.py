__author__ = 'Naheed'
import numpy as np
from copy import deepcopy
from sympy import Matrix


class homology:
    def __init__(self, delta_k, delta_kplus1):
        self.delta_k = delta_k
        self.delta_kplus1 = delta_kplus1
        # self.rowops_delta_k = []*len(self.delta_k.shape[1]) # delta_k is used to compute the kernel dim. i.e num of dependent columns
        # self.rowops_delta_kplus1 = []*len(self.delta_kplus1.shape[0]) # delta_kplus1 is used to compute the image dim, i.e the num of independent rows

    def compute_homology(self):
        pass

    def compute_betti_number(self):
        # kernelDim = Matrix(self.delta_k.T).rref()[0].shape[0] - len(Matrix(self.delta_k.T).rref()[1])
        # imageDim = len(Matrix(self.delta_kplus1).rref()[1])
        # print Matrix(d_kplus1).rref()
        # return kernelDim - imageDim

        reduced_mat, operationsets = self.get_reduced_reform(self.delta_kplus1)

        print self.delta_k
        print self.delta_kplus1
        self.printop_list(operationsets)
        print reduced_mat

    def printop_list(self, listof_ops):
        print '\n'
        for i in listof_ops:
            print i, ' '
        print '\n'

    def get_reduced_reform(self, input_matrix):
        """
        Performs Elementary Row reduction- converting input_matrix into reduced-row-echelon-form(rref)
        :param input_matrix:
        :return: Reduced form of the input_matrix, list of operation sets performed on each rows
        """

        operations_list = [[i] for i in range(input_matrix.shape[0])]  # [] * number of rows

        matrix = deepcopy(input_matrix)
        numrows, numcols = matrix.shape

        i, j = 0, 0
        while True:
            if i >= numrows or j >= numcols:
                break

            if matrix[i, j] == 0:
                nonzeroRow = i
                while nonzeroRow < numrows and matrix[nonzeroRow, j] == 0:
                    nonzeroRow += 1

                if nonzeroRow == numrows:
                    j += 1
                    continue

                self.rowSwap(matrix, operations_list, i, nonzeroRow)

            pivot = matrix[i, j]
            self.scaleRow(matrix, operations_list, i, 1.0 / pivot)

            for otherRow in range(0, numrows):
                if otherRow == i:
                    continue
                if matrix[otherRow, j] != 0:
                    scaleAmt = -matrix[otherRow, j]
                    self.rowCombine(matrix, operations_list, otherRow, i, scaleAmt)

            i += 1
            j += 1

        return matrix, operations_list

    def rowSwap(self, A, rowops_list, i, j):
        temp, tempop = np.copy(A[i, :]), rowops_list[i]
        A[i, :], rowops_list[i] = A[j, :], rowops_list[j]
        A[j, :], rowops_list[j] = temp, rowops_list[i]

    def scaleRow(self, A, rowops_list, i, c):
        A[i, :] *= c * np.ones(A.shape[1])
        rowops_list[i] = (rowops_list[i], c)

    def rowCombine(self, A, rowops_list, addTo, scaleRow, scaleAmt):
        A[addTo, :] += scaleAmt * A[scaleRow, :]
        rowops_list[addTo] = (rowops_list[addTo], (rowops_list[scaleRow], scaleAmt), '+')

    def numPivotCols(A):
        z = np.zeros(A.shape[0])
        return [np.all(A[:, j] == z) for j in range(A.shape[1])].count(False)

    def numPivotRows(A):
        z = np.zeros(A.shape[1])
        return [np.all(A[i, :] == z) for i in range(A.shape[0])].count(False)
