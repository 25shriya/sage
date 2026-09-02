r"""
Schur Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.data_structures.blas_dict import convert_remove_zeroes
from sage.libs.lrcalc import lcalc
from sage.matrix.constructor import matrix
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.partition import Partitions
from sage.categories.tensor import tensor

class SupersymFunctionAlgebra_schur(super_sfa.SuperSymAlgebra_generic):
    r"""
    Schur supersymmetric functions.
    
    The Schur supersymmetric function defined on variables `\mathbf{x}` and
    `\mathbf{y}`, `s_\lambda(\mathbf{x} \mid \mathbf{y})`, is given as

    .. MATH::
    
        s_\lambda(\mathbf{x} \mid \mathbf{y}) = \text{det}([h_{\lambda_i - i + j}(\mathbf{x} \mid \mathbf{y})]_{i,j=1}^{l})

    where, `\lambda = (\lambda_1, \lambda_2, \ldots, \lambda_l)` is a partition, and `l = \mid \lambda \mid`.
    
    These form a non-graded multiplicative basis for the ring of supersymmetric
    functions.
    
    REFERENCES:
    
    - [BHS25]_
    
    INPUT:

    - ``Supersym`` -- ring of supersymmetric functions
    
    EXAMPLES::
    
        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: p = s.s()
        sage: p
        Supersymmetric functions over Rational Field in the Schur basis
    """
    def __init__(self, SuperSym):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.s()
            sage: TestSuite(s).run()
        """
        super().__init__(SuperSym=SuperSym, graded=False, prefix='s', basis_name='Schur')

    def product_on_basis(self, left, right):
        r"""
        Return the product of ``left`` and ``right``.
        
        INPUT:

        - ``left``, ``right`` -- partitions

        TESTS::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ).s()
            sage: a = s([2,1]) + 1; a
            s[] + s[2, 1]
            sage: a^2   # indirect doctest
            s[] + 2*s[2, 1] + s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1]
                + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]
        """
        return self.element_class(self, convert_remove_zeroes(lrcalc.mult(left, right),
                                                                  self.base_ring()))
    
    def coproduct_on_basis(self, mu):
        r"""
        Return the coproduct of ``self(mu)``.

        Here ``self`` is the basis of Schur supersymmetric functions.

        INPUT:

        - ``mu`` -- a partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: s.coproduct_on_basis([2])
            s[] # s[2] + s[1] # s[1] + s[2] # s[]
        """
        T = self.tensor_square()
        return T.element_class(T, convert_remove_zeroes(lrcalc.coprod(mu, all=1),
                                                        self.base_ring()))

    def lift_on_basis(self, part, basis_name='homogeneous'):
        r"""
        Return the Schur basis in terms of homogeneous or elementary basis
        for a given partition, ``part``.

        INPUT:

        - ``part`` -- a partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: from sage.combinat.partition import _Partitions
            sage: part = _Partitions([4,4,2])
            sage: s.lift_on_basis(part)
            0
        """
        l = len(part)
        req_matrix = [[0]*l]*l
        if basis_name == 'homogeneous':
            h = self.realization_of().h()
            for i in range(l):
                for j in range(l):
                    req_matrix[i][j] = h[part[i] - i + j]
            req_matrix = matrix(req_matrix)
            return req_matrix.det()
        elif basis_name == 'elementary':
            e = self.realization_of().e()
            part = part.conjugate()
            l = len(part)
            for i in range(l):
                for j in range(l):
                    req_matrix[i][j] = e[part[i] - i + j]
            req_matrix = matrix(req_matrix)
            return req_matrix.det()

    def supersym_to_sym(self, part):
        r"""
        Return the Schur basis in terms of Schur symmetric basis
        for a given partition, ``part``.

        INPUT:

        - ``part`` -- a partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: from sage.combinat.partition import _Partitions
            sage: part = _Partitions([6,5])
            sage: s.supersym_to_sym(part)
            0
        """
        R = self.base_ring()
        s = SymmetricFunctions(R).s()
        T = s.tensor_square()
        req_sum = R.zero()
        for mu in Partitions(sum(part), ending=part):
            for nu in Partitions(sum(part), ending=part):
                nu = nu.conjugate()
                req_sum += lrcalc.lrcoef(mu, nu, part) * T.sum_of_monomials((mu, nu))
        return req_sum

    def murnaghan_nakayama(self, part):
        r"""
        Return the Schur basis in terms of Schur symmetric basis
        for a given partition, ``part``.

        INPUT:

        - ``part`` -- a partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: from sage.combinat.partition import Partition
            sage: part = Partition([6,5])
            sage: s.murnaghan_nakayama(part)
            1/302400*p[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] +
            1/15120*p[2, 1, 1, 1, 1, 1, 1, 1, 1, 1] +
            1/2016*p[2, 2, 1, 1, 1, 1, 1, 1, 1] +
            1/720*p[2, 2, 2, 1, 1, 1, 1, 1] +
            1/576*p[2, 2, 2, 2, 1, 1, 1] +
            1/20160*p[3, 1, 1, 1, 1, 1, 1, 1, 1] +
            1/720*p[3, 2, 1, 1, 1, 1, 1, 1] +
            1/288*p[3, 2, 2, 1, 1, 1, 1] + 1/144*p[3, 2, 2, 2, 1, 1] -
            1/576*p[3, 2, 2, 2, 2] + 1/360*p[3, 3, 1, 1, 1, 1, 1] +
            1/72*p[3, 3, 2, 2, 1] - 1/108*p[3, 3, 3, 1, 1] +
            1/108*p[3, 3, 3, 2] - 1/2520*p[4, 1, 1, 1, 1, 1, 1, 1] +
            1/72*p[4, 3, 1, 1, 1, 1] - 1/36*p[4, 3, 3, 1] +
            1/48*p[4, 4, 1, 1, 1] - 1/48*p[4, 4, 3] -
            1/450*p[5, 1, 1, 1, 1, 1, 1] - 1/120*p[5, 2, 1, 1, 1, 1] -
            1/120*p[5, 2, 2, 2] + 1/90*p[5, 3, 1, 1, 1] +
            1/30*p[5, 3, 2, 1] + 1/90*p[5, 3, 3] + 1/20*p[5, 4, 1, 1] +
            1/25*p[5, 5, 1] - 1/180*p[6, 1, 1, 1, 1, 1] -
            1/36*p[6, 2, 1, 1, 1] - 1/36*p[6, 3, 1, 1] +
            1/36*p[6, 3, 2] + 1/30*p[6, 5] - 1/168*p[7, 1, 1, 1, 1] -
            1/28*p[7, 2, 1, 1] - 1/56*p[7, 2, 2] - 1/21*p[7, 3, 1] -
            1/28*p[7, 4]
        """
        R = self.base_ring()
        Sym = SymmetricFunctions(R)
        p_sym = Sym.p()
        s_sym = Sym.s()
        res = p_sym(s_sym[part])
        p = self.realization_of().p()
        return p._from_dict(res.monomial_coefficients())
