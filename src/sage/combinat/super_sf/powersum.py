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
            prod = self.parent().base_ring().one()
            for p in part:
                x_gens = [alphabet_x+str(i).format(i) for i in range(1, p+1)]
                y_gens = [alphabet_y+str(i).format(i) for i in range(1, p+1)]
                variables = x_gens + y_gens
                R = PolynomialRing(self.parent().base_ring(), variables)
                R_gens = R.gens_dict()
                x_gens = [R_gens[gen] for gen in x_gens]
                y_gens = [R_gens[gen] for gen in y_gens]
                l = min(len(x_gens), len(y_gens))
                req_sum = R.zero()
                for j in range(l):
                    req_sum += x_gens[j] ** p - y_gens[j] ** p
                prod *= req_sum
            return prod

        def plethysm(self, part):
            R = self.base_ring()
            phi = R.one()
            p_sym = SymmetricFunctions(R).p()
            for p in part:
                a = tensor([p_sym[p], R.one()])
                b = (-1) ** (p) * tensor([R.one(), p_sym[p]])
                phi *= a + b
            return phi.section()
            
