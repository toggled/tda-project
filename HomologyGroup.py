__author__ = 'Naheed'

from copy import deepcopy
from sympy import Matrix, Symbol
from src.BoundaryGroup import KthBoundaryGroup
import numpy as np


class Homology:
    def __init__(self, boundary_group_k, boundary_group_kplus1):
        """
        :type boundary_group_k: KthBoundaryGroup
        :type boundary_group_kplus1: KthBoundaryGroup
        """
        assert isinstance(boundary_group_k, KthBoundaryGroup)
        assert isinstance(boundary_group_kplus1, KthBoundaryGroup)
        self.delta_k = boundary_group_k.get_transformation_matrix()

        self.delta_kplus1 = boundary_group_kplus1.get_transformation_matrix()
        print self.delta_k
        print self.delta_kplus1
        self.kernel_basis_delta_k = boundary_group_k.get_columnobjects()
        self.image_basis_delta_kplus1 = boundary_group_kplus1.get_rowobjects()

        # self.rowops_delta_k = []*len(self.delta_k.shape[1]) # delta_k is used to compute the kernel dim. i.e num of dependent columns
        # self.rowops_delta_kplus1 = []*len(self.delta_kplus1.shape[0]) # delta_kplus1 is used to compute the image dim, i.e the num of independent rows

    # def compute_kth_homology_groups(self):
    #     """
    #     Compute the K-th Homology groups (Hk)
    #     :return: Tuple (Basis for kernel space(Zk) and Basis for Image space (Bk)
    #     TODO:   BUG HERE.
    #     """
    #
    #     # print self.kernel_basis_delta_k
    #     # print self.image_basis_delta_kplus1
    #
    #     def get_kernel_space(reduced_matrix, z_k):
    #         zeros = np.zeros(reduced_matrix.shape[1])
    #         l = [i for i in range(reduced_matrix.shape[0]) if np.all(reduced_matrix[i, :] == zeros) == True]
    #         return z_k[l]
    #
    #     reduced_mat2, operationsets2 = self.get_reduced_reform(self.delta_k.T)
    #     kernel_spaces = self.interpret_rref(self.kernel_basis_delta_k, operationsets2)
    #     kernel_generators = get_kernel_space(reduced_mat2, np.array(kernel_spaces))
    #     print 'x ',operationsets2
    #     print 'hi ',reduced_mat2
    #     for x in kernel_generators:
    #         print x
    #     def get_basis_imagespace(reduced_matrix): # Find the pivot columns
    #         #zeros = np.zeros(reduced_matrix.shape[1])
    #         #l = [i for i in range(reduced_matrix.shape[0]) if np.all(reduced_matrix[i, :] == zeros) == False]
    #         #return b_k[l]
    #         print self.image_basis_delta_kplus1
    #         basis_list = []
    #         for i in range(reduced_mat.shape[0]):
    #             for k in range(i,reduced_mat.shape[1]):
    #                 if reduced_mat[i][k] == 1:
    #                     basis_list.append(self.image_basis_delta_kplus1[k])
    #                     break
    #         return basis_list
    #
    #     reduced_mat, operationsets = self.get_reduced_reform(self.delta_kplus1)  # Find the image space (span of Bk)
    #     print self.delta_kplus1
    #     print reduced_mat
    #     print operationsets
    #
    #     image_spaces = self.interpret_rref(self.image_basis_delta_kplus1, operationsets)
    #
    #     generators_imagespace = get_basis_imagespace(reduced_mat) # To find the image space we don't need the basis for rows.
    #                                                              # we need to find pivot columns i.e minimum columns that are enough to represent image space
    #     import sys
    #     sys.exit(1)
    #     # print image_spaces
    #
    #     return (kernel_generators, generators_imagespace)

    def compute_kth_homology_groups(self):

        def get_kernel_space(reduced_matrix, z_k):
            zeros = np.zeros(reduced_matrix.shape[0], dtype=np.int)
            l = [i for i in range(reduced_matrix.shape[1]) if np.all(reduced_matrix[:, i] == zeros) == True]
            return list(z_k[l])

        def get_basis_imagespace(reduced_mat, image_basis_space):  # Find the pivot columns

            basis_list = []
            for i in range(reduced_mat.shape[0]):
                for k in range(i, reduced_mat.shape[1]):
                    if reduced_mat[i][k] == 1:
                        basis_list.append(image_basis_space[i])
                        break
            return basis_list

        rowops, colops, dk, dk_1 = self.simultaneousReduce(self.delta_k, self.delta_kplus1)
        # print dk
        dk_1, rowops = self.get_reduced_reform(dk_1, rowops)
        # print dk_1

        kernel_spaces = self.interpret_rref(self.kernel_basis_delta_k, colops)

        kernel_generators = get_kernel_space(dk, np.array(kernel_spaces))
        print kernel_generators
        self.play_with_kernel_space(kernel_generators)

        image_spaces = self.interpret_rref(self.image_basis_delta_kplus1, rowops)
        image_generators = get_basis_imagespace(dk_1, np.array(image_spaces))
        # print image_generators
        return (kernel_generators, image_generators)

    def play_with_kernel_space(self, kernel_generators):
        for i in range(len(kernel_generators)):
            for j in range(i + 1, len(kernel_generators)):
                if i == 2 or j == 2:
                    continue
                diff = kernel_generators[i] - kernel_generators[j]
                print diff, (kernel_generators[i], kernel_generators[j])



    def compute_betti_number(self):

        reduced_mat, operationsets = self.get_reduced_reform(self.delta_kplus1)

        #print 'print bk+1: \n', self.delta_kplus1
        # self.printop_list(operationsets)
        #print 'reduced form of bk+1: \n', reduced_mat



        reduced_mat2, operationsets2 = self.get_reduced_reform(self.delta_k.T)
        # self.printop_list(operationsets2)
        # print 'print bk: \n', self.delta_k.T, '\n'
        #print 'reduced form of bk: \n', reduced_mat2



        # print "betti: ", (self.delta_k.shape[1] - self.numPivotRows(reduced_mat2)) - self.numPivotRows(reduced_mat)
        return (self.delta_k.shape[1] - self.numPivotRows(reduced_mat2)) - self.numPivotRows(reduced_mat)

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

    def interpret_rref(self, repr_row_simplices, listof_ops):
        """
        Compute the string representation of the operations matrix
        :param repr_row_simplices: Row representation along the rows in the transformation matrix
        :param listof_ops: List of operations performed on the operation matrix to make it rref
        :return: row representation of the rows in the operation matrix
        """
        # string_repr = [""] * len(listof_ops)
        string_repr = None
        string_rowrepr = [0] * len(repr_row_simplices)

        def process_next_op(operations_list):
            if len(operations_list) == 3:  # process row combine
                #print type(operations_list[0]),'hi'
                idxleft, left = process_next_op(operations_list[0])
                idxright, right = process_next_op(operations_list[1])
                # stack.append(left + '=' +left+'+'+right)
                # return left  # concatenate left string and right string
                string_repr[idxleft] = left + right
                return idxleft, string_repr[idxleft]

            if len(operations_list) == 2:  # process row scaling
                if operations_list[1] == -1:
                    idx, op = process_next_op(operations_list[0])
                    # stack.append(op+'='+'-'+op)
                    # return op
                    string_repr[idx] = -op
                    return idx, string_repr[idx]

                return process_next_op(operations_list[0])
            if len(operations_list) == 1:  # base case
                # if type(operations_list) is list:
                #     string_repr[operations_list[0]] = str(operations_list[0])
                #     return operations_list[0], string_repr[operations_list[0]]
                return operations_list[0], Symbol(repr_row_simplices[operations_list[
                    0]])  # operation list is operation on row in transformation matrix
                # We need to return the corresponding symbolic representation
                # of the simplex at that row.
                    # if type(operations_list) is tuple:
                    #     idx, op = process_next_op(operations_list[0])
                    #     return idx,string_repr[idx]

        #print 'represnetation: ', len(string_rowrepr)
        for idx, operations in enumerate(listof_ops):
            # print operations
            string_repr = [Symbol(row) for row in repr_row_simplices]
            process_next_op(operations)
            string_rowrepr[idx] = string_repr[idx]
            #print string_repr
        return string_rowrepr

    def simultaneousReduce(self, A, B):
        # print 'inside simultaneous reduce\n'

        row_operations_list = [[i] for i in range(B.shape[0])]  # [] * number of rows
        col_operations_list = [[i] for i in range(A.shape[1])]  # [] * number of cols

        if A.shape[1] != B.shape[0]:
            raise Exception("Matrices have the wrong shape.")

        numRows, numCols = A.shape

        i, j = 0, 0
        operation = 1
        while True:
            if i >= numRows or j >= numCols:
                break

            if A[i, j] == 0:
                nonzeroCol = j
                while nonzeroCol < numCols and A[i, nonzeroCol] == 0:
                    nonzeroCol += 1

                if nonzeroCol == numCols:
                    i += 1
                    continue

                self.colSwap(A, col_operations_list, j, nonzeroCol)
                self.rowSwap(B, row_operations_list, j, nonzeroCol)

            pivot = A[i, j]
            self.scaleCol(A, col_operations_list, j, 1 * pivot)
            self.scaleRow(B, row_operations_list, j, 1 * pivot)

            for otherCol in range(0, numCols):
                if otherCol == j:
                    continue
                if A[i, otherCol] != 0:
                    scaleAmt = -A[i, otherCol]
                    self.colCombine(A, col_operations_list, otherCol, j, scaleAmt)
                    self.rowCombine(B, row_operations_list, j, otherCol, -scaleAmt)

            i += 1
            j += 1

        # print 'exiting simultaneous reduce\n'
        return (row_operations_list, col_operations_list, A, B)

    def get_reduced_reform(self, input_matrix, operations_list=[]):
        """
        Performs Elementary Row reduction- converting input_matrix into reduced-row-echelon-form(rref)
        :param input_matrix:
        :return: Reduced form of the input_matrix, list of operation sets performed on each rows. i.e the operation list
        """
        if operations_list == []:
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
        temp, tempop = np.copy(A[i, :]), list(deepcopy(rowops_list[i]))
        A[i, :], rowops_list[i] = A[j, :], rowops_list[j]
        assert isinstance(tempop, list)
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

    def colSwap(self, A, colops_list, i, j):
        temp, tempop = np.copy(A[:, i]), list(deepcopy(colops_list[i]))
        A[:, i], colops_list[i] = A[:, j], colops_list[j]
        A[:, j], colops_list[j] = temp, tempop

    def scaleCol(self, A, colops_list, i, c):
        A[:, i] *= c * np.ones(A.shape[0], dtype=np.int)
        colops_list[i] = (colops_list[i], c)

    def colCombine(self, A, colops_list, addTo, scaleCol, scaleAmt):
        A[:, addTo] += scaleAmt * A[:, scaleCol]
        colops_list[addTo] = (colops_list[addTo], (colops_list[scaleCol], scaleAmt), '+')
