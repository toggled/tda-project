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

    def compute_homology_groups(self):
        pass

    def compute_betti_number(self):

        reduced_mat, operationsets = self.get_reduced_reform(self.delta_kplus1)

        print 'print bk+1: \n', self.delta_kplus1
        # self.printop_list(operationsets)
        print 'reduced form of bk+1: \n', reduced_mat

        res = self.interpret_rref(reduced_mat, operationsets)
        print res

        reduced_mat2, operationsets2 = self.get_reduced_reform(self.delta_k.T)
        # self.printop_list(operationsets2)
        print 'print bk: \n', self.delta_k.T, '\n'
        print 'reduced form of bk: \n', reduced_mat2
        print self.interpret_rref(reduced_mat, operationsets)

        print "betti: ", (self.delta_k.shape[1] - self.numPivotRows(reduced_mat2)) - self.numPivotRows(reduced_mat)

    def printop_list(self, listof_ops):
        """
        Print the operations matrix
        :param listof_ops: operations matrix
        :return: None
        """
        print 'printing operations list: \n'
        for i in listof_ops:
            print i, ' '
        print '\n'

    def interpret_rref(self, reduced_matrix, listof_ops):
        """
        Compute the string representation of the operations matrix
        :param reduced_matrix: Reduced Row echelon form (rref) of the matrix
        :param listof_ops: List of operations performed on the operation matrix to make it rref
        :return: row representation of the rows in the operation matrix
        """
        string_repr = [""] * len(listof_ops)
        string_rowrepr = [[""]] * len(listof_ops)

        def process_next_op(operations_list):
            if len(operations_list) == 3:  # process row combine
                idxleft, left = process_next_op(operations_list[0])
                idxright, right = process_next_op(operations_list[1])
                # stack.append(left + '=' +left+'+'+right)
                # return left  # concatenate left string and right string
                string_repr[idxleft] = string_repr[idxleft] + "+" + string_repr[idxright]
                return idxleft, string_repr[idxleft]

            if len(operations_list) == 2:  # process row scaling
                if operations_list[1] == -1:
                    idx, op = process_next_op(operations_list[0])
                    # stack.append(op+'='+'-'+op)
                    # return op
                    string_repr[idx] = str(["-" + str([op])])
                    return idx, string_repr[idx]

                return process_next_op(operations_list[0])
            if len(operations_list) == 1:  # base case
                if type(operations_list) is list:
                    string_repr[operations_list[0]] = str(operations_list[0])
                    return operations_list[0], string_repr[operations_list[0]]

                    # if type(operations_list) is tuple:
                    #     idx, op = process_next_op(operations_list[0])
                    #     return idx,string_repr[idx]

        for idx, operations in enumerate(listof_ops):
            process_next_op(operations)
            string_rowrepr[idx] = string_repr[idx]
            # print string_repr

        return string_rowrepr

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
            self.scaleRow(matrix, operations_list, i, 1 * pivot)

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
        temp, tempop = np.copy(A[i, :]), deepcopy(rowops_list[i])
        A[i, :], rowops_list[i] = A[j, :], rowops_list[j]
        A[j, :], rowops_list[j] = temp, tempop

    def scaleRow(self, A, rowops_list, i, c):
        A[i, :] *= c * np.ones(A.shape[1], dtype=np.int)
        rowops_list[i] = (rowops_list[i], c)

    def rowCombine(self, A, rowops_list, addTo, scaleRow, scaleAmt):
        A[addTo, :] += scaleAmt * A[scaleRow, :]
        rowops_list[addTo] = (rowops_list[addTo], (rowops_list[scaleRow], scaleAmt), '+')

    def numPivotRows(self, A):
        z = np.zeros(A.shape[1])
        return [np.all(A[i, :] == z) for i in range(A.shape[0])].count(False)
