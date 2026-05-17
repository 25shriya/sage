from . import super_sfa
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class SupersymFunctionAlgebra_powersum(super_sfa.SuperSymAlgebra_multiplicative):
    def __init__(self, base_ring):
        super.__init__(self, base_ring)

    class Element(super_sfa.SuperSymAlgebra_multiplicative.Element):
        def expand(self, i):
            variables = self.alpha_gens + self.x_gens + self.y_gens
            R = PolynomialRing(self.base_ring, variables)
            l = min(len(self.x_gens), len(self.y_gens))
            req_sum = R.zero()
            for j in range(l):
                req_sum += self.x_gens[j] ** i - self.y_gens[j] ** i
            return req_sum