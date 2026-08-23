r"""
Schur Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.data_structures.blas_dict import convert_remove_zeroes
import sage.libs.lrcalc.lrcalc as lrcalc
from sage.matrix.constructor import matrix

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
        
