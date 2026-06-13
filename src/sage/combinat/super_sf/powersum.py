r"""
Powersum basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.partition import Partition
from sage.categories.tensor import tensor
from sage.combinat.sf.sf import SymmetricFunctions

class SupersymFunctionAlgebra_powersum(super_sfa.SuperSymAlgebra_multiplicative):
    def __init__(self, Supersym):
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=False)

    def coproduct_on_generators(self, n):
            part = Partition([n])
            return tensor([self[part], self[part]])

    def antipode_on_basis(self, n):
            part = Partition([n])
            return -self[part]

    class Element(super_sfa.SuperSymAlgebra_multiplicative.Element):
        def expand(self, part=[], alphabet_x='x', alphabet_y='y'):
            part.sort(reverse=True)
            part = Partition(part)
            R = self.base_ring()
            p_sym = SymmetricFunctions(R).p()
            self_parts = self.monomial_coefficients()
            sum = R.zero()
            for k in self_parts:
                prod = R.one()
                for p in part:
                    prod *= p_sym[k].expand(p, alphabet=alphabet_x) - p_sym[k].expand(p, alphabet=alphabet_y)
                sum += self_parts[k] * prod
            return sum
       
        def plethysm(self, part):
            R = self.base_ring()
            phi = R.one()
            p_sym = SymmetricFunctions(R).p()
            for p in part:
                a = tensor([p_sym[p], R.one()])
                b = (-1) ** (p) * tensor([R.one(), p_sym[p]])
                phi *= a + b
            return phi.section()