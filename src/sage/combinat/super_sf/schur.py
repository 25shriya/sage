r"""
Schur Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.partition import _Partitions, Partition
from sage.categories.tensor import tensor
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.misc_c import prod
from sage.misc.lazy_attribute import lazy_attribute

class SupersymFunctionAlgebra_schur(super_sfa.SuperSymAlgebra_generic):
    def __init__(self, SuperSym):
        super().__init__(SuperSym=SuperSym, graded=False, prefix='s', basis_name='schur')

    def 