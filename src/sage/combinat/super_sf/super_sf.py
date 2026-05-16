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
from unittest import TestSuite
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.principal_ideal_domains import PrincipalIdealDomains

class SuperSymmetricFunctions(Parent, UniqueRepresentation):
    def __init__(self, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: R = SuperSymmetricFunctions(QQ)
            sage: TestSuite(R).run()
        """
        self._base_ring = R
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
        return "Supersymmetric functions over %s" % self._base_ring