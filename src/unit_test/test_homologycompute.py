__author__ = 'Naheed'
from src.simplex import KSimplex, SimplicialComplex
from src.BoundaryGroup import KthBoundaryGroup
from src.homology import another_bettinum, bettiNumber
from HomologyGroup import homology

if __name__ == '__main__':
    # FilePath = "../../data/simplices.out"
    FilePath = "../../data/test_simplexfromfile.txt"

    sigma = SimplicialComplex()
    sigma.add_simplex_fromfile(FilePath)

    # for k in range(sigma.maxK):  # not maxK + 1
    for k in range(0, 1):  # not maxK
        # print str(k + 1) + '-th Transformation Matrix:'
        Bk = KthBoundaryGroup(k)
        Bk.construct_from_simplex(sigma.get_allkth_simplices(k + 1))  # Send all the k+1 simplices to it
        print Bk.get_transformation_matrix()

        m = k + 1
        # print str(m + 1) + '-th Transformation Matrix:'
        Bk1 = KthBoundaryGroup(m)
        Bk1.construct_from_simplex(sigma.get_allkth_simplices(m + 1))  # Send all the k+1 simplices to it
        #print Bk1.get_transformation_matrix()

        # Bk.print_columnobjects()
        # Bk.print_rowobjects()
        homology_obj = homology(Bk.get_transformation_matrix(), Bk1.get_transformation_matrix())
        print 'H', str(k + 1), ' = ', str(homology_obj.compute_betti_number())
        # print 'H', str(k + 1), ' = ', bettiNumber(Bk.get_transformation_matrix(), Bk1.get_transformation_matrix())
        # print 'H', str(k + 1), ' = ', another_bettinum(Bk.get_transformation_matrix(), Bk1.get_transformation_matrix())
        # if result == 0:  # if K-th Betti number is 0. the rest will be so
        #break
