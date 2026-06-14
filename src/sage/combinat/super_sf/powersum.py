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
        r"""
        Class for methods associated with powersum supersymmetric functions.

        INPUT:

        - ``Supersym`` -- the ring of supersymmetric functions

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: TestSuite(p).run()
        """
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p._repr_()
            Supersymmetric functions over Rational Field in the powersum basis
        """
        return "%s in the powersum basis" % (self.realization_of())

    def coproduct_on_generators(self, n):
        r"""
        Return coproduct on generators for power sums `p_i`
        (for `i > 0`).

        INPUT:

        - ``n`` -- a positive integer

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p.coproduct_on_generators(3)
            B[[3]] # B[[3]]
        """
        part = Partition([n])
        return tensor([self[part], self[part]]) # Shouldn't this take in a partition?

    def antipode_on_basis(self, part):
        r"""
        Return the antipode of ``self[partition]``.

        INPUT:

        - ``partition`` -- a list or partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: from sage.combinat.partition import Partition
            sage: part = Partition([4,4,2])
            sage: p.antipode_on_basis(part)
            -B[[4, 4, 2]]
        """
        if not isinstance(part, Partition):
            part.sort(reverse=True)
            part = Partition(part)
        return -self[part]

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
                sage: p = s.p()
                sage: p[[4]].expand(3)
                x1^4 + x2^4 + x3^4 - y1^4 - y2^4 - y3^4
            """
            self_parts = self.monomial_coefficients()
            x_gens = [alphabet_x+str(i).format(i) for i in range(1, n+1)]
            y_gens = [alphabet_y+str(i).format(i) for i in range(1, n+1)]
            variables = x_gens + y_gens
            R = PolynomialRing(self.base_ring(), variables)
            R_gens = R.gens_dict()
            x_gens1 = [R_gens[gen] for gen in x_gens]
            y_gens1 = [R_gens[gen] for gen in y_gens]
            res = R.zero()
            for part in self_parts:
                prod = R.one()
                for p in part:
                    sum = R.zero()
                    for i in range(0, n):
                        sum += x_gens1[i] ** p - y_gens1[i] ** p
                    prod *= sum                    
                res += self_parts[part] * prod
            return res
       
        def plethysm(self, part):
            r"""
                Return the plethysm of ``self`` with ``part``.

                INPUT:

                - ``part`` -- a partition
            """
            R = self.base_ring()
            phi = R.one()
            p_sym = SymmetricFunctions(R).p()
            for p in part:
                a = tensor([p_sym[p], R.one()])
                b = (-1) ** (p) * tensor([R.one(), p_sym[p]])
                phi *= a + b
            return phi.section() # Debug