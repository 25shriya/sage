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
from sage.misc.misc_c import prod

class SupersymFunctionAlgebra_powersum(super_sfa.SuperSymAlgebra_multiplicative):
    r"""
    Powersum supersymmetric functions.

    The powersum supersymmetric function defined on variables `\mathbb{x}` and
    `\mathbb{y}`, `p_i(\mathbb{x} \mid \mathbb{y})`, is given as

    .. MATH::

        p_i(\mathbb{x} \mid \mathbb{y}) = \sum_{k=1}^\infty x_k^i - y_k^i

    These form a non-graded multiplicative basis for the ring of supersymmetric
    functions.

    REFERENCES:

    - [BHS25]_

    INPUT:

        - ``Supersym`` -- the ring of supersymmetric functions
    """
    def __init__(self, Supersym):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: TestSuite(p).run()
        """
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=False,
                                                   prefix='p', basis_name='powersum')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p._repr_()
            'Supersymmetric functions over Rational Field in the powersum basis'
        """
        return "%s in the %s basis" % (self.realization_of(), self.basis_name)

    def coproduct_on_generators(self, part):
        r"""
        Return coproduct on generators for power sums `p_\lambda`
        (for partition `\lambda`).

        INPUT:

        - ``part`` -- a partition or list

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p.coproduct_on_generators([3, 2])
            p[3, 2] # p[3, 2]
        """
        if not isinstance(part, Partition):
            part = Partition(part)
        return tensor([self[part], self[part]])

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
            -p[4, 4, 2]
        """
        if not isinstance(part, Partition):
            part.sort(reverse=True)
            part = Partition(part)
        return -self[part]

    def plethysm(self, part):
        r"""
        Return the plethysm of ``self`` with ``part``.

        INPUT:

        - ``part`` -- a partition
        """
        R = self.base_ring()
        phi = R.one()
        p_sym = SymmetricFunctions(R).p()
        phi = prod([tensor([p_sym[p], p_sym.one()]) + (-1 ** p) * tensor([p_sym.one(), p_sym[p]]) for p in part])
        return phi #.section() - Debug

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
            res = sum([self_parts[part] *
                       prod([sum([x_gens1[i] ** p - y_gens1[i] ** p
                                  for i in range(n)])
                                  for p in part])
                                  for part in self_parts])
            return res

# Debug .section()
# I added self._base to sym. It's messing with the ring.
# TestSuite() cases?
