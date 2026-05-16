r"""
Bases of Supersymmetric functions
"""
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partition

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

class SuperSymAlgebra_generic(CombinatorialFreeModule):
    def __init__(self, base):
        CombinatorialFreeModule.__init__(self, base, base().index_set(), prefix="s", category=SuperSymmetricFunctionsBases(base))

    class Element(CombinatorialFreeModule.Element):
        def __repr__(self):
            return "%s(%s)" % (self.parent().name(), self)

class SuperSymAlgebra_multiplicative(SuperSymAlgebra_generic):
    def __init__(self, base):
        CombinatorialFreeModule.__init__(self, base, base().index_set(), prefix="s", category=SuperSymmetricFunctionsBases(base))

    def product_on_basis(self, left, right):
        m = list(left) + list(right)
        m.sort(reverse=True)
        return self.monomial(Partition(m))

    def coproduct_on_basis(self, mu):
        T = self.tensor_square()
        return T.prod(self.coproduct_on_generators(p) for p in mu)

    def _repr_(self):
        return "Multiplicative basis of %s" % self.base()
    
    class Element(SuperSymAlgebra_generic.Element):
        def __init__(self, x_gens, y_gens, alpha_gens):
            self.x_gens = x_gens
            self.y_gens = y_gens
            self.alpha_gens = alpha_gens

        def __repr__(self):
            return "%s(%s)" % (self.parent().name(), self)