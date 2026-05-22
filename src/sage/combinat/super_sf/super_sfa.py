r"""
Bases of Supersymmetric functions
"""
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains

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

class SuperSymAlgebra_generic(CombinatorialFreeModule):
    def __init__(self, base):
        CombinatorialFreeModule.__init__(self, base, category=SuperSymmetricFunctionsBases(base))

    class Element(CombinatorialFreeModule.Element):
        """
            Generic bases of supersymmetric functions.
        """
        pass

class SuperSymAlgebra_multiplicative(SuperSymAlgebra_generic):
    def __init__(self, base):
        CombinatorialFreeModule.__init__(self, base, category=SuperSymmetricFunctionsBases(base))
    
    class Element(SuperSymAlgebra_generic.Element):
        """
            Multiplicative bases of supersymmetric functions. 
        """
        pass