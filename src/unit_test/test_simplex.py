from __future__ import absolute_import
from src.simplex import KSimplex, SimplicialComplex
from src.boundaryoperator import Boundary
from src.BoundaryGroup import KthBoundaryGroup

def test_manual_simplex():
    sim = KSimplex([1, 2, 3])
    C = SimplicialComplex()
    C.add_simplex(sim)
    for i in range(2):
        sim = KSimplex(range(i + 1))
        C.add_simplex(sim)
    C.add_simplex(KSimplex([2, 3]))
    print C


def test_file_simplex():
    sigma = SimplicialComplex()
    sigma.add_simplex_fromfile('test_simplexfromfile.txt')
    print sigma


def test_boundary_op():
    delta = Boundary()
    delta.compute_boundary(KSimplex([1, 2, 3]))
    print delta


def test_nested_boundary():
    delta_k = Boundary()
    delta_k.compute_boundary(KSimplex([1, 2, 3]))
    print delta_k
    for sign, kmin1_simpl in delta_k.get_boundary():
        # print kmin1_simpl.k
        delta_kmin1 = Boundary()
        delta_kmin1.compute_boundary(kmin1_simpl)
        print delta_kmin1


def test_nested_boundary_simplicialcomplex():
    sigma = SimplicialComplex()
    sigma.add_simplex_fromfile('test_simplexfromfile.txt')

    for k in range(sigma.maxK + 1):
        print str(k) + '-th Chain group:'
        for k_sim in sigma.get_allkth_simplices(k):
            delta_k = Boundary()
            delta_k.compute_boundary(k_sim)
            print 'Boundary of: ', str(k_sim), ': ', str(delta_k)


def test_kth_boundary_group():
    sigma = SimplicialComplex()
    sigma.add_simplex_fromfile('test_simplexfromfile.txt')

    for k in range(sigma.maxK):  # not maxK + 1
        print str(k) + '-th Boundary group:'
        Bk = KthBoundaryGroup(k)
        Bk.construct_from_simplex(sigma.get_allkth_simplices(k + 1))  # Send all the k+1 simplices to it
        print Bk

if __name__ == "__main__":
    # test_file_simplex()
    # test_boundary_op()
    # test_nested_boundary()
# test_nested_boundary_simplicialcomplex()
# test_kth_boundary_group()
