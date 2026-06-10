r"""
Bases of Supersymmetric functions
"""
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
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
        """
        return "Category of bases of %s" % self.base()

    def super_categories(self):
        R = self.base().base_ring()
        cat = HopfAlgebras(R)
        categories = [self.base().Realizations(),
                      cat.Commutative().WithBasis(),
                      cat.Graded().Realizations()]
        if R in PrincipalIdealDomains:
            categories.append(UniqueFactorizationDomains())
        return categories

class FilteredSuperSymBases(Category_realization_of_parent):
    def super_categories(self):
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Filtered()
        return [SuperSymmetricFunctionsBases(self.base()), cat]

class GradedSuperSymBases(Category_realization_of_parent):
    def super_categories(self):
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Graded()
        return [FilteredSuperSymBases(self.base()), cat]

class SuperSymAlgebra_generic(CombinatorialFreeModule):
    def __init__(self, SuperSym=None, graded=True):
        R = SuperSym.base_ring()
        from sage.categories.commutative_rings import CommutativeRings
        if R not in CommutativeRings():
            raise TypeError("argument R must be a commutative ring")
        if graded:
            cat = GradedSuperSymBases(SuperSym)
        else:  # Right now, there are no non-filtered bases. Do we have filtered bases in Supersym case?
            cat = FilteredSuperSymBases(SuperSym)
        CombinatorialFreeModule.__init__(self, R, category=cat)

    def __getitem__(self, c):
        return self.monomial(c)

class SuperSymAlgebra_multiplicative(SuperSymAlgebra_generic):
    def __init__(self, base):
        CombinatorialFreeModule.__init__(self, base, category=SuperSymmetricFunctionsBases(base))