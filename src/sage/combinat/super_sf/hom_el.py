from . import super_sfa

class SupersymFunctionAlgebra_hom_el(super_sfa.SuperSymAlgebra_multiplicative):
    def __init__(self, base_ring, basis='homogeneous'):
        self.basis =  basis
        super.__init__(self, base_ring)
    
    class Element(super_sfa.SuperSymAlgebra_multiplicative):
        def expand(self, i):
            # yet to finish
            return i
