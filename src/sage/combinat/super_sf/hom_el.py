from . import super_sfa
from itertools import combinations_with_replacement, combinations
from sage.combinat.partition import Partitions, Partition
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class SupersymFunctionAlgebra_hom_el(super_sfa.SuperSymAlgebra_multiplicative):
    def __init__(self, Supersym, basis_name='homogeneous'):
        self.basis_name =  basis_name
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=False)
    
    class Element(super_sfa.SuperSymAlgebra_multiplicative):
        def expand(self, part=[], alphabet_x='x', alphabet_y='y'):
            basis_name = self.parent().basis_name
            res = self.base_ring.one()
            for p in part:
                x_gens = [alphabet_x+str(i).format(i) for i in range(1, p+1)]
                y_gens = [alphabet_y+str(i).format(i) for i in range(1, p+1)]
                variables = x_gens + y_gens
                R = PolynomialRing(self.base_ring, variables)
                req_sum = R.zero()
                parts = Partitions(p, length=2)
                for k in parts:
                    for seq1 in combinations(p, k[1]):
                        for seq2 in combinations_with_replacement(p, k[0]):
                            if basis_name == 'homogeneous':
                                new_seq1 = seq1.sort()
                                new_seq2 = seq2.sort(reverse=True)
                                prod = R.one()
                                for i in range(len(new_seq1)):
                                    prod *= y_gens[i]
                                for j in range(len(new_seq2)):
                                    prod *= x_gens[j]
                                    req_sum += prod
                            elif basis_name == 'elementary':
                                new_seq1 = seq1.sort(reverse=True)
                                new_seq2 = seq2.sort()
                                prod = R.one()
                                for i in range(len(new_seq2)):
                                    prod *= y_gens[i]
                                for j in range(len(new_seq1)):
                                    prod *= x_gens[j]
                    res *= req_sum
            return res

        def coproduct_on_basis(self, n):
            T = self.tensor_square()
            return T.sum_of_monomials((Partition([i]), Partition([n-i])) for i in range(1, n+1))

# What are super_categories? Figure out the math
# Current problems - creating an element seems impossible. How does the base ring get set? How can you define elements for each partition?
# The above two questions will be answered through symmetric functions.

                
                            

                
            
