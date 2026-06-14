r"""
Bases of Supersymmetric functions
"""
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.combinat.partition import Partition
from sage.rings.integer_ring import ZZ

class SuperSymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of bases of supersymmetric functions.

    EXAMPLES::

        sage: SuperSymmetricFunctionsBases()    
        Category of bases of supersymmetric functions
    """
    def _repr_(self):
        r"""
        Return the string representation of this category.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: s1 = SuperSymmetricFunctionsBases(s)
            sage: s1._repr_()
            Category of bases of Supersymmetric functions over Rational Field
        """
        return "Category of bases of %s" % self.base()

    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: s1 = SuperSymmetricFunctionsBases(s)
            sage: s1
            Category of bases of Supersymmetric functions over Rational Field
            sage: s1.super_categories()
            [Category of realizations of Supersymmetric functions over Rational Field,
             Category of commutative Hopf algebras with basis over Rational Field,
             Join of Category of realizations of Hopf algebras over Rational Field and Category of graded algebras over Rational Field and Category of graded coalgebras over Rational Field,
             Category of unique factorization domains]
        """
        R = self.base().base_ring()
        cat = HopfAlgebras(R)
        categories = [self.base().Realizations(),
                      cat.Commutative().WithBasis(),
                      cat.Graded().Realizations()]
        if R in PrincipalIdealDomains:
            categories.append(UniqueFactorizationDomains())
        return categories

class FilteredSuperSymBases(Category_realization_of_parent):
    r"""
    The category of filtered bases of supersymmetric functions.

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: from sage.combinat.super_sf.super_sfa import FilteredSuperSymBases
        sage: f = FilteredSuperSymBases(s)
        sage: f
        Category of filtered bases of Supersymmetric functions over Rational Field
    """
    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import FilteredSuperSymBases
            sage: f = FilteredSuperSymBases(s)
            sage: f._repr_()
            Category of filtered bases of Supersymmetric functions over Rational Field
        """
        return "Category of filtered bases of %s" % self.base()

    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import FilteredSuperSymBases
            sage: f = FilteredSuperSymBases(s)
            sage: f.super_categories()
            [Category of bases of Supersymmetric functions over Rational Field,
             Category of commutative filtered Hopf algebras with basis over Rational Field]
        """
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Filtered()
        return [SuperSymmetricFunctionsBases(self.base()), cat]

class GradedSuperSymBases(Category_realization_of_parent):
    r"""
    The category of graded bases of supersymmetric functions.

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: from sage.combinat.super_sf.super_sfa import GradedSuperSymBases
        sage: f = GradedSuperSymBases(s)
        sage: f
        Category of graded bases of Supersymmetric functions over Rational Field
    """
    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import FilteredSuperSymBases
            sage: f = GradedSuperSymBases(s)
            sage: f._repr_()
            Category of graded bases of Supersymmetric functions over Rational Field
        """
        return "Category of graded bases of %s" % self.base()
    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import GradedSuperSymBases
            sage: g = GradedSuperSymBases(s)
            sage: g.super_categories()
            [Category of filtered bases of Supersymmetric functions over Rational Field,
             Category of commutative graded Hopf algebras with basis over Rational Field]
        """
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Graded()
        return [FilteredSuperSymBases(self.base()), cat]

class SuperSymAlgebra_generic(CombinatorialFreeModule):
    r"""
    Abstract class for Supersymmetric Function Algebras.
    """
    def __init__(self, SuperSym=None, graded=True):
        r"""
        Initialize a supersymmetric function algebra.

        INPUT:

        - ``SuperSym`` -- the ring of supersymmetric functions (default: ``None``)
        - ``graded`` -- boolean (default: ``True``); if ``True``, then the basis is
          considered to be graded, otherwise the basis is filtered

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import SuperSymAlgebra_generic
            sage: p = s.p()
            sage: isinstance(p, SuperSymAlgebra_generic)
            True
        """
        self.SuperSym = SuperSym
        R = SuperSym.base_ring()
        from sage.categories.commutative_rings import CommutativeRings
        if R not in CommutativeRings():
            raise TypeError("argument R must be a commutative ring")
        if graded:
            cat = GradedSuperSymBases(SuperSym)
        else:  # Right now, there are no non-graded bases. Do we have filtered bases in Supersym case?
            cat = FilteredSuperSymBases(SuperSym)
        CombinatorialFreeModule.__init__(self, R, category=cat)

    def __getitem__(self, c):
        r"""
        Return the monomial corresponding to the given list/partition/integer.

        INPUT:

        - ``c`` -- list, integer or partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p[5,4,3]
            B[[5, 4, 3]]
            sage: from sage.combinat.partition import Partition
            sage: part = Partition([3,2])
            sage: p[part]
            B[[3, 2]]
        """
        if not isinstance(c, Partition):
            if c not in ZZ:
                list(c).sort(reverse=True)
            else:
                c = [c]
            c = Partition(mu=c)
        return self.monomial(c)

class SuperSymAlgebra_multiplicative(SuperSymAlgebra_generic):
    r"""
    Abstract class of multiplicative supersymmetric function algebras.

    A realization `h` of the ring of supersymmetric functions is multiplicative
    if for a partition `\lambda = (\lambda_1,\lambda_2,\ldots)` we have 
    `h_\lambda = h_{\lambda_1} h_{\lambda_2} \cdots`.
    """
    def product_on_basis(self, left, right):
        r"""
        Return the product of ``left`` and ``right``.

        INPUT:

        - ``left``, ``right`` -- partitions

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: e = s.h_e('elementary')
            sage: e.product_on_basis([3,2], [6,5])
            B[[6, 5, 3, 2]]
        """
        m = list(left) + list(right)
        m.sort(reverse=True)
        return self.monomial(Partition(m))

    def coproduct_on_basis(self, mu):
        r"""
        Return the coproduct on a basis element for multiplicative bases.

        INPUT:

        - ``mu`` -- a partition
        """
        # Debug this function
        # T = self.tensor_square()
        # return T.prod(self.coproduct_on_generators(p) for p in mu)
        return