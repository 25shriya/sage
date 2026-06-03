r"""
Powersum basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.partition import Partition

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
            part = Partition(part.sort(reverse=True))
            prod = self.base_ring.one()
            for p in part:
                x_gens = [alphabet_x+str(i).format(i) for i in range(1, p+1)]
                y_gens = [alphabet_y+str(i).format(i) for i in range(1, p+1)]
                variables = x_gens + y_gens
                R = PolynomialRing(self.base_ring, variables)
                l = min(len(self.x_gens), len(self.y_gens))
                req_sum = R.zero()
                for j in range(l):
                    req_sum += self.x_gens[j] ** p - self.y_gens[j] ** p
                prod *= req_sum
            return prod

        def plethysm(self, i):
            T = self.tensor_square()
            return
            # return T.sum_of_monomials((Partition([i]), Partition([])), (Partition()))
        
        # write embedding function - It's upper triangular
