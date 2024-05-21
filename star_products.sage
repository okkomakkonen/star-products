from sage.combinat.q_analogues import q_binomial
from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator

from rich.progress import track

def is_mds(C):
    """Return True if and only if C is MDS"""
    
    return C.minimum_distance() == C.length() - C.dimension() + 1

def star_product(C, D):
    """Return the star product of two codes"""

    G = matrix([c.pairwise_product(d) for c in C.generator_matrix() for d in D.generator_matrix()])

    return LinearCode(G)

def sum_dimensions(C, l):
    
    G = C.generator_matrix()
    n = C.length()
    F = C.base_field()
    
    total = 0

    # instead of running through linear codes, we run through all reduced echelon matrices, since this is equivalent
    for D in track(reduced_echelon_matrix_iterator(F, l, n), total=q_binomial(n, l, F.order())):
        
        M = matrix([c.pairwise_product(d) for c in G for d in D])
        
        total += M.rank()
        
    return total