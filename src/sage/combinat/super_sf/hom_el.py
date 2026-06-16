r"""
Homogeneous and Elementary Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from itertools import combinations_with_replacement, combinations
from sage.combinat.partition import Partitions, Partition
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.functions.other import factorial

class SupersymFunctionAlgebra_hom_el(super_sfa.SuperSymAlgebra_multiplicative):
    r"""
    Homogeneous and elementary supersymmetric functions.

    The homogeneous supersymmetric function defined on variables `\mathbb{x}` and
    `\mathbb{y}`, `h_k(\mathbb{x} \mid \mathbb{y})`, is given as

    .. MATH::

        h_k(\mathbb{x} \mid \mathbb{y}) = \sum_{a+b=k}\sum_{i_1 \geq \ldots \geq i_b \\ j_1 < \ldots < j_a}
        y_{j_1} \ldots y_{j_a} x_{i_1} \ldots x_{i_b}.

    The elementary supersymmetric function defined on variables `\mathbb{x}` and
    `\mathbb{y}`, `e_k(\mathbb{x} \mid \mathbb{y})`, is given as

    .. MATH::

        e_k(\mathbb{x} \mid \mathbb{y}) = \sum_{a+b=k}\sum_{i_1 > \ldots > i_b \\ j_1 \leq \ldots \leq j_a}
        y_{j_1} \ldots y_{j_a} x_{i_1} \ldots x_{i_b}.
    

    These form a multiplicative non-graded basis for the ring of supersymmetric
    functions.

    REFERENCES:

    - [BHS25]_

    INPUT:

    - ``Supersym`` -- the ring of supersymmetric functions
    - ``basis_name`` -- string (default: ``'homogeneous'``); one of the following

        * ``'homogeneous'`` - homogeneous basis of supersymmetric functions
        * ``'elementary'`` - elementary basis of supersymmetric functions
    """
    def __init__(self, Supersym, basis_name='homogeneous'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: TestSuite(h).run()
        """
        self.basis_name =  basis_name
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: e = s.e()
            sage: e._repr_()
            Supersymmetric functions over Rational Field in the elementary basis
        """
        return "%s in the %s basis" % (self.realization_of(), self.basis_name)

    def coproduct_on_generators(self, n):
        r"""
        Return the coproduct on `h_i`.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: h.coproduct_on_generators(5)
            B[[1]] # B[[4]] + B[[2]] # B[[3]] + B[[3]] # B[[2]] + B[[4]] # B[[1]] + B[[5]] # B[[]]
        """
        T = self.tensor_square()
        return T.sum_of_monomials((Partition([i]), Partition([n-i])) for i in range(1, n+1))

    def lift_on_gens(self, n): # Debug
        r"""
        Return the homogeneous or elementary basis in terms of powersum
        basis for a given ``n``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: h.change_of_basis(3)
        """
        # basis_name = self.parent().basis_name
        # parts = Partitions(n)
        # Supersym = self.parent().realization_of()
        # R = self.base_ring()
        # sum = Supersym.zero()

        # def epsilon(part):
        #     return (-1) ** (n - len(part))

        # def z(part):
        #     prod = Supersym.one()
        #     for p in part:
        #         m = part.to_exp()[p-1]
        #         prod *= (p ** m) * factorial(p)
        #     return prod

        # if basis_name == 'homogeneous':
        #     for part in parts:
        #         sum += (1 / z(part)) * Supersym.p()[part]
        # elif basis_name == 'elementary':
        #     for part in parts:
        #         sum += (1 / z(part)) * (epsilon(part)) * Supersym.p()[part]
        # return sum
        return

    class Element(super_sfa.SuperSymAlgebra_multiplicative.Element):
        def expand(self, n, alphabet_x='x', alphabet_y='y'):
            r"""
            Expand the supersymmetric function ``self`` as a supersymmetric
            polynomial in ``n`` variables.

            INPUT:

            - ``n`` -- nonnegative integer
            - ``alphabet_x`` -- (default: ``'x'``) a variable for the expansion `x`
            - ``alphabet_y`` -- (default: ``'y'``) a variable for the expansion `y`

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ)
                sage: h = s.h()
                sage: h[3,2].expand(3)
                324*x1^2*x2^2*y1^2 + 162*x1^2*x2*y1^2
            """
            basis_name = self.parent().basis_name
            res = self.base_ring().one()
            monomial_coeff = self.monomial_coefficients()
            x_gens = [alphabet_x+str(i).format(i) for i in range(1, n+1)]
            y_gens = [alphabet_y+str(i).format(i) for i in range(1, n+1)]
            variables = x_gens + y_gens
            R = PolynomialRing(self.base_ring(), variables)
            R_gens = R.gens_dict()
            x_gens = [R_gens[gen] for gen in x_gens]
            y_gens = [R_gens[gen] for gen in y_gens]
            req_sum = R.zero()
            fin_res = R.zero()
            for k in monomial_coeff:
                for ki in k:
                    for p in range(1, ki+1):
                        for seq1 in combinations(range(n), (ki - p)):
                            for seq2 in combinations_with_replacement(range(n), p):
                                if basis_name == 'homogeneous':
                                    prod = R.one()
                                    for i in range(len(seq1)):
                                        prod *= y_gens[i]
                                    for j in range(len(seq2)):
                                        prod *= x_gens[j]
                                    req_sum += prod
                                elif basis_name == 'elementary':
                                    prod = R.one()
                                    for i in range(len(seq2)):
                                        prod *= y_gens[i]
                                    for j in range(len(seq1)):
                                        prod *= x_gens[j]
                                    req_sum += prod
                    res *= req_sum
                fin_res += monomial_coeff[k] * res
            return fin_res
