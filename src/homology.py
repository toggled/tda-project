import numpy
import numpy.linalg
import sys
from sympy import Matrix


def rowSwap(A, i, j):
    temp = numpy.copy(A[i, :])
    A[i, :] = A[j, :]
    A[j, :] = temp


def colSwap(A, i, j):
    temp = numpy.copy(A[:, i])
    A[:, i] = A[:, j]
    A[:, j] = temp


def scaleCol(A, i, c):
    # print (c * numpy.ones(A.shape[0], dtype = int)).dtype
    # print A.dtype
    A[:, i] *= c * numpy.ones(A.shape[0], dtype=int)


def scaleRow(A, i, c):
    A[i, :] *= c * numpy.ones(A.shape[1], dtype=int)


def colCombine(A, addTo, scaleCol, scaleAmt):
    A[:, addTo] += scaleAmt * A[:, scaleCol]


def rowCombine(A, addTo, scaleRow, scaleAmt):
    A[addTo, :] += scaleAmt * A[scaleRow, :]


def simultaneousReduce(A, B):
    print 'inside simultaneous reduce\n'
    # print A.shape
    # print B.shape
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

            colSwap(A, j, nonzeroCol)
            rowSwap(B, j, nonzeroCol)
            print 'operation: ', str(operation), '\n'
            print A
            print B
        pivot = A[i, j]
        scaleCol(A, j, 1 * pivot)
        scaleRow(B, j, 1 * pivot)
        print A
        print B
        for otherCol in range(0, numCols):
            if otherCol == j:
                continue
            if A[i, otherCol] != 0:
                scaleAmt = -A[i, otherCol]
                colCombine(A, otherCol, j, scaleAmt)
                rowCombine(B, j, otherCol, -scaleAmt)
                print A
                print B
        i += 1
        j += 1
    print 'exiting simultaneous reduce\n'
    return A, B


# Doing Elementary Row reduction. converting B into reduced-row-echelon-form(rref)
def finishRowReducing(B):
    numRows, numCols = B.shape

    i, j = 0, 0
    while True:
        if i >= numRows or j >= numCols:
            break

        if B[i, j] == 0:
            nonzeroRow = i
            while nonzeroRow < numRows and B[nonzeroRow, j] == 0:
                nonzeroRow += 1

            if nonzeroRow == numRows:
                j += 1
                continue

            rowSwap(B, i, nonzeroRow)

        pivot = B[i, j]
        scaleRow(B, i, 1 * pivot)

        for otherRow in range(0, numRows):
            if otherRow == i:
                continue
            if B[otherRow, j] != 0:
                scaleAmt = -B[otherRow, j]
                rowCombine(B, otherRow, i, scaleAmt)

        i += 1;
        j += 1

    return B


def numPivotCols(A):
    z = numpy.zeros(A.shape[0])
    return [numpy.all(A[:, j] == z) for j in range(A.shape[1])].count(False)


def numPivotRows(A):
    z = numpy.zeros(A.shape[1])
    return [numpy.all(A[i, :] == z) for i in range(A.shape[0])].count(False)


def bettiNumber(d_k, d_kplus1):
    A, B = numpy.copy(d_k), numpy.copy(d_kplus1)

    simultaneousReduce(A, B)
    print A
    print B
    finishRowReducing(B)

    print "After finish row reducing:\n", B

    dimKChains = A.shape[1]
    # print(dimKChains)
    kernelDim = dimKChains - numPivotCols(A)
    # print(kernelDim)
    imageDim = numPivotRows(B)
    # print(imageDim)

    return kernelDim - imageDim


def another_bettinum(d_k, d_kplus1):
    kernelDim = Matrix(d_k.T).rref()[0].shape[0] - len(Matrix(d_k.T).rref()[1])
    imageDim = len(Matrix(d_kplus1).rref()[1])
    print 'image space: \n', Matrix(d_kplus1).rref()
    print 'kernel space: \n', Matrix(d_k.T).rref()
    return kernelDim - imageDim


if __name__ == '__main__':
    '''
    bd0 = numpy.array([[0,0,0,0,0]])
    bd1 = numpy.array([[-1,-1,-1,-1,0,0,0,0], [1,0,0,0,-1,-1,0,0],    [0,1,0,0,1,0,-1,-1], [0,0,1,0,0,1,1,0], [0,0,0,1,0,0,0,1]])
    bd2 = numpy.array([[1,1,0,0],[-1,0,1,0],[0,-1,-1,0],
         [0,0,0,0],[1,0,0,1],[0,1,0,-1],
         [0,0,1,1],[0,0,0,0]])
    bd3 = numpy.array([[-1],[1],[-1],[1]])
    '''
    bd0 = numpy.array([[0, 0, 0, 0]])
    # bd1 = numpy.array([[-1,0,0,-1,-1], [1,-1,0,0,0], [0,1,-1,0,1], [0,0,1,1,0]])
    # bd1 = numpy.array([[-1,0,0,-1,-1,0], [1,-1,0,0,0,-1], [0,1,-1,0,1,0], [0,0,1,1,0,1]])

    # bd1 = numpy.array([[1,-1,0,0,0,-1], [-1,0,0,1,-1,0], [0,1,-1,0,1,0], [0,0,1,-1,0,1]])

    # bd2 = numpy.array([[0,0,1, -1, 1,0],[1,1,0,0,-1,0]]).T

    bd2 = numpy.array([[0, 0, 1, -1, 1, 0]]).T

    # bd0 = numpy.array([[0,0,0,0,0]])
    # bd1 = numpy.array([[-1,0,0,-1,1], [1,-1,0,0,0], [0,1,-1,0,-1], [0,0,1,1,0], [0,0,0,0,0]])
    # bd1 = numpy.array(
    #     [[-1, 0, 1, 1, 0, 1, 0, 0], [1, -1, 0, 0, 0, 0, -1, 0], [0, 1, -1, 0, 0, 0, 0, -1], [0, 0, 0, -1, -1, 0, 1, 1],
    #      [0, 0, 0, 0, 1, -1, 0, 0]])
    # bd2 = numpy.array(
    #     [[1, 0, 0, 1], [1, 0, 1, 0], [1, -1, 0, 0], [0, 1, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, -1, 1],
    #      [0, 1, 1, 0]])
    # bd 1 = numpy.array([[-1,0,0,-1,-1], [1,-1,0,0,0], [0,1,-1,0,1], [0,0,1,1,0]])

    bd1 = numpy.array([[-1, 0, 0, 0, -1, 0, -1, -1],
                       [1, -1, 0, 0, 0, -1, 0, 0],
                       [0, 1, -1, 0, 0, 0, 1, 0],
                       [0, 0, 1, -1, 0, 1, 0, 1],
                       [0, 0, 0, 1, 1, 0, 0, 0]])
    bd2 = numpy.array([[1],
                       [1],
                       [0],
                       [0],
                       [0],
                       [0],
                       [-1],
                       [0]])

    print("Example complex from post")
    # print("0th homology: %d" % bettiNumber(bd0,bd1))
    # print "another way to compute betti number: ",another_bettinum(bd1,bd2)
    # print Matrix(bd1.T).rref()

    print("1st homology: %d" % bettiNumber(bd1, bd2))
    #print "another way to compute betti number: ", another_bettinum(bd1, bd2)

    # print("2nd homology: %d" % bettiNumber(bd2, bd3))
    # print "another way to compute betti number: ", another_bettinum(bd2, bd3)

    import sys

    sys.exit(1)
    mobiusD1 = numpy.array([
        [-1, -1, -1, -1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, -1, -1, -1, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, -1, -1, 0],
        [0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
    ])

    mobiusD2 = numpy.array([
        [1, 0, 0, 0, 1],
        [0, 0, 0, 1, 0],
        [-1, 0, 0, 0, 0],
        [0, 0, 0, -1, -1],
        [0, 1, 0, 0, 0],
        [1, -1, 0, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 1, 1, 0, 0],
        [0, 0, -1, 1, 0],
        [0, 0, 1, 0, 0],
    ])

    print("Mobius Band")
    print("1st homology: %d" % bettiNumber(mobiusD1, mobiusD2))
    print "another way to compute betti number: ", another_bettinum(mobiusD1, mobiusD2)
