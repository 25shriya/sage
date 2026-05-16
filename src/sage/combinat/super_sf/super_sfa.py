r"""
Bases of Supersymmetric functions
"""
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule

class SuperSymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of bases of supersymmetric functions.
    sage: SuperSymmetricFunctionsBases()    
    Category of bases of supersymmetric functions
    """
    def _repr_(self):
        r"""
        Return the string representation of this category.
        """
        return "Category of bases of %s" % self.base()

class SuperSymAlgebra_generic(CombinatorialFreeModule):
    class Element(CombinatorialFreeModule.Element):
        def __repr__(self):
            return "%s(%s)" % (self.parent().name(), self)

class SuperSymAlgebra_multiplicative(SuperSymAlgebra_generic):
    def __init__(self, base):
        CombinatorialFreeModule.__init__(self, base, base().index_set(), prefix="s", category=SuperSymmetricFunctionsBases(base))

    def _repr_(self):
        return "Multiplicative basis of %s" % self.base()
    
    class Element(SuperSymAlgebra_generic.Element):
        def __repr__(self):
            return "%s(%s)" % (self.parent().name(), self)