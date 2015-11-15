__author__ = 'Naheed'
from src.simplex import KSimplex, SimplicialComplex
from src.BoundaryGroup import KthBoundaryGroup
from src.homology import another_bettinum

if __name__ == '__main__':
    sigma = SimplicialComplex()
    sigma.add_simplex_fromfile('test_simplexfromfile.txt')

    # for k in range(sigma.maxK):  # not maxK + 1
    for k in range(-1, sigma.maxK - 1):  # not maxK + 1
        print str(k + 1) + '-th Transformation Matrix:'
        Bk = KthBoundaryGroup(k)
        Bk.construct_from_simplex(sigma.get_allkth_simplices(k + 1))  # Send all the k+1 simplices to it
        print Bk.get_transformation_matrix()

        m = k + 1
        print str(m + 1) + '-th Transformation Matrix:'
        Bk1 = KthBoundaryGroup(m)
        Bk1.construct_from_simplex(sigma.get_allkth_simplices(m + 1))  # Send all the k+1 simplices to it
        print Bk1.get_transformation_matrix()

        # Bk.print_columnobjects()
        # Bk.print_rowobjects()
        print 'H', str(k + 1), ' = ', another_bettinum(Bk.get_transformation_matrix(), Bk1.get_transformation_matrix())
