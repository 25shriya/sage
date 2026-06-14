from . import super_sfa
from itertools import combinations_with_replacement, combinations
from sage.combinat.partition import Partitions, Partition
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.functions.other import factorial

class SupersymFunctionAlgebra_hom_el(super_sfa.SuperSymAlgebra_multiplicative):
    def __init__(self, Supersym, basis_name='homogeneous'):
        self.basis_name =  basis_name
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=False)

    def coproduct_on_generators(self, n):
            T = self.tensor_square()
            return T.sum_of_monomials((Partition([i]), Partition([n-i])) for i in range(1, n+1))

    class Element(super_sfa.SuperSymAlgebra_multiplicative.Element):
        def expand(self, part=[], alphabet_x='x', alphabet_y='y'):
            basis_name = self.parent().basis_name
            res = self.base_ring().one()
            monomial_coeff = self.monomial_coefficients()
            for k in monomial_coeff:
                for ki in k:
                    for p in part:
                        x_gens = [alphabet_x+str(i).format(i) for i in range(1, ki+1)]
                        y_gens = [alphabet_y+str(i).format(i) for i in range(1, ki+1)]
                        variables = x_gens + y_gens
                        R = PolynomialRing(self.base_ring(), variables)
                        R_gens = R.gens_dict()
                        x_gens = [R_gens[gen] for gen in x_gens]
                        y_gens = [R_gens[gen] for gen in y_gens]
                        req_sum = R.zero()
                        parts = Partitions(p, length=2)
                        for k in parts:
                            for seq1 in combinations(range(p), k[1]):
                                for seq2 in combinations_with_replacement(range(p), k[0]):
                                    if basis_name == 'homogeneous':
                                        new_seq1 = sorted(seq1)
                                        new_seq2 = sorted(seq2, reverse=True)
                                        prod = R.one()
                                        for i in range(len(new_seq1)):
                                            prod *= y_gens[i]
                                        for j in range(len(new_seq2)):
                                            prod *= x_gens[j]
                                        req_sum += prod
                                    elif basis_name == 'elementary':
                                        new_seq1 = sorted(seq1, reverse=True)
                                        new_seq2 = sorted(seq2)
                                        prod = R.one()
                                        for i in range(len(new_seq2)):
                                            prod *= y_gens[i]
                                        for j in range(len(new_seq1)):
                                            prod *= x_gens[j]
                                        req_sum += prod
                            res *= req_sum
                    return res

        def change_of_basis(self, n):
            basis_name = self.parent().basis_name
            parts = Partitions(n)
            Supersym = self.parent().SuperSym
            R = self.base_ring()
            sum = R.zero()

            def epsilon(part):
                return (-1) ** (n - len(part))

            def z(part):
                prod = R.one()
                for p in part:
                    m = part.to_exp()[p-1]
                    prod *= (p ** m) * factorial(p)
                return prod

            if basis_name == 'homogeneous':
                for part in parts:
                    sum += (1 / z(part)) * Supersym.p()[part]
            elif basis_name == 'elementary':
                for part in parts:
                    sum += (1 / z(part)) * (epsilon(part)) * Supersym.p()[part]
            return sum



                
                            

                
            
