r"""
Schur Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.data_structures.blas_dict import convert_remove_zeroes
import sage.libs.lrcalc.lrcalc as lrcalc
from sage.matrix.constructor import matrix
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.partition import Partitions
from sage.categories.tensor import tensor

class SupersymFunctionAlgebra_schur(super_sfa.SuperSymAlgebra_generic):
    def __init__(self, SuperSym):
        super().__init__(SuperSym=SuperSym, graded=False, prefix='s', basis_name='Schur')

    def product_on_basis(self, left, right):
        return self.element_class(self, convert_remove_zeroes(lrcalc.mult(left, right),
                                                                  self.base_ring()))
    
    def coproduct_on_basis(self, mu):
        T = self.tensor_square()
        return T.element_class(T, convert_remove_zeroes(lrcalc.coprod(mu, all=1),
                                                        self.base_ring()))

    def lift_on_basis(self, part, basis_name='homogeneous'):
        req_matrix = matrix.identity(l)
        if basis_name == 'homogeneous':
            h = self.realization_of().h()
            l = len(part)
            for i in range(l):
                for j in range(l):
                    req_matrix[i][j] = h[part[i] - i + j]
            return req_matrix.det()
        elif basis_name == 'elementary':
            e = self.realization_of().e()
            part = part.conjugate()
            l = len(part)
            for i in range(l):
                for j in range(l):
                    req_matrix[i][j] = e[part[i] - i + j]
            return req_matrix.det()

    def supersym_to_sym(self, part):
        R = self.base_ring()
        s = SymmetricFunctions(R).s()
        T = s.tensor_square()
        req_sum = R.zero()
        for mu in Partitions(sum(part), ending=part):
            for nu in Partitions(sum(part), ending=part):
                nu = nu.conjugate()
                req_sum += lrcalc.lrcoef(mu, nu, part) * T.sum_of_monomials((mu, nu))
        return req_sum
