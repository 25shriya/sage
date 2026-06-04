r"""
Supersymmetric functions, with their realizations

AUTHORS:

- Shriya M
"""
# ****************************************************************************
#       Copyright (C) 2026 Shriya M <25shriya at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.rational_field import QQ
from sage.categories.fields import Fields
from sage.categories.rings import Rings
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.combinat.super_sf.powersum import SupersymFunctionAlgebra_powersum
from sage.combinat.super_sf.hom_el import SupersymFunctionAlgebra_hom_el

class SuperSymmetricFunctions(UniqueRepresentation, Parent):
    def __init__(self, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: R = SuperSymmetricFunctions(QQ)
            sage: TestSuite(R).run()
        """
        self._base = R
        assert R in Fields() or R in Rings()
        cat = GradedHopfAlgebras(R).Commutative().Cocommutative()
        if R in PrincipalIdealDomains():
            cat &= UniqueFactorizationDomains()
        Parent.__init__(self, category=cat.WithRealizations())

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: R = SuperSymmetricFunctions(QQ)
            sage: R._repr_()
            'Supersymmetric functions over Rational Field'
        """
        return "Supersymmetric functions over %s" % self._base

    def powersum(self):
        return SupersymFunctionAlgebra_powersum(self)

    p = powersum

    def homogeneous_elementary(self):
        return SupersymFunctionAlgebra_hom_el(self)

    h_e = homogeneous_elementary