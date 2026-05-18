from . import super_sfa
from itertools import combinations_with_replacement, combinations
from sage.combinat.partition import Partitions
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class SupersymFunctionAlgebra_hom_el(super_sfa.SuperSymAlgebra_multiplicative):
    def __init__(self, base_ring, basis='homogeneous'):
        self.basis =  basis
        self.element = self.Element(self.basis)
        super.__init__(self, base_ring)
    
    class Element(super_sfa.SuperSymAlgebra_multiplicative):
        def __init__(self, basis):
            self.basis = basis

        def expand(self, part=[], alphabet_x='x', alphabet_y='y', alpha='alpha'):
            res = self.base_ring.one()
            for p in part:
                x_gens = [alphabet_x+str(i).format(i) for i in range(1, p+1)]
                y_gens = [alphabet_y+str(i).format(i) for i in range(1, p+1)]
                alpha_gens = [alpha+str(i).format(i) for i in range(1, p+1)]
                variables = x_gens + y_gens + alpha_gens
                R = PolynomialRing(self.base_ring, variables)
                req_sum = R.zero()
                parts = Partitions(p, length=2)
                for k in parts:
                    for seq1 in combinations(p, k[1]):
                        for seq2 in combinations_with_replacement(p, k[0]):
                            new_seq1 = seq1.sort()
                            new_seq2 = seq2.sort(reverse=True)
                            prod = R.one()
                            for i in range(len(new_seq1)):
                                prod *= y_gens[i]
                            for j in range(len(new_seq2)):
                                prod *= x_gens[j]
                                req_sum += prod
                    res *= req_sum
            return res

# How do I add alpha gens given the negative indexing?
# Isn't the above double supersym? Should that be a different class?
# Above case is just for homogeneous supersym. Make it double + get elementary into the picture
                
                            

                
            
