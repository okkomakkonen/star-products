"""
Implement some basic methods for studying star products of codes
"""

from itertools import permutations, product
from collections import defaultdict

from sage.combinat.q_analogues import q_binomial
from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator
from sage.coding.grs_code import GeneralizedReedSolomonCode
from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel

from rich.progress import track

def is_mds(C):
    """Return True if and only if C is MDS"""
    
    return C.minimum_distance() == C.length() - C.dimension() + 1

def star_product(C, D):
    """Return the star product of two codes"""

    G = matrix([c.pairwise_product(d) for c in C.generator_matrix() for d in D.generator_matrix()])

    return LinearCode(G)

def support(M):

    return sum(m != 0 for m in M.transpose())

def sum_dimensions(C, l, verbose=True):
    
    G = C.generator_matrix()
    n = C.length()
    F = C.base_field()

    if verbose:
        print(f"Testing with l={l}:\n{G}")

    total = 0

    # instead of running through linear codes, we run through all reduced echelon matrices, since this is equivalent
    for D in track(reduced_echelon_matrix_iterator(F, l, n), total=q_binomial(n, l, F.order())):

        M = matrix([c.pairwise_product(d) for c in G for d in D])
        
        total += M.rank()

    if verbose:
        print("Sum of dimensions:", total)
        
    return total

def dimension_distribution(C, l, verbose=True):
    
    G = C.generator_matrix()
    n = C.length()
    F = C.base_field()

    if verbose:
        print(f"Testing with l={l}:\n{G}")

    distribution = defaultdict(int)
    total = 0

    # instead of running through linear codes, we run through all reduced echelon matrices, since this is equivalent
    for D in track(reduced_echelon_matrix_iterator(F, l, n), total=q_binomial(n, l, F.order())):

        M = matrix([c.pairwise_product(d) for c in G for d in D])
        
        distribution[M.rank()] += 1
        total += M.rank()

    distribution = dict(distribution)

    if verbose:
        print("Dimension_distribution:", distribution)
        print("Sum of dimensions:", total)
        
    return distribution

def sum_supports(F, n, l):
    
    total = 0

    # instead of running through linear codes, we run through all reduced echelon matrices, since this is equivalent
    for D in track(reduced_echelon_matrix_iterator(F, l, n), total=q_binomial(n, l, F.order())):
        
        total += support(D)
        
    return total

def sum_lower_bound(F, n, l, k):

    total = 0

    # instead of running through linear codes, we run through all reduced echelon matrices, since this is equivalent
    for D in track(reduced_echelon_matrix_iterator(F, l, n), total=q_binomial(n, l, F.order())):

        supp = support(D)

        total += min(supp, min(supp, k) + l - 1)
        
    return total

def closed_sum_support(F, n, l):

    q = F.order()

    return sum(s * binomial(n, s) * sum((-1)^i * q_binomial(s - i, l, q)*binomial(s, s - i) for i in range(s - l + 1)) for s in range(l, n + 1))

def random_mds_code(F, n, k, num_tries=1000):

    for _ in range(num_tries):

        G = random_matrix(F, k, n)
        C = LinearCode(G)

        if is_mds(C):
            return C

    raise Exception("number of tries exceeded")

def is_rs(C):

    if not is_mds(C):
        return False

    F = C.base_field()
    n = C.length()

    for alpha in permutations(F.list(), n):

        if C == GeneralizedReedSolomonCode(alpha, dimension=C.dimension()):
            return True
        
    return False

def is_grs(C):

    if not is_mds(C):
        return False

    F = C.base_field()
    n = C.length()

    assert n >= 2

    elements = F.list()[1:-1]

    for evals in permutations(elements, int(n - 2)):

        alpha = vector(F, [F.zero(), F.one()] + list(evals))

        for cols in product(F.list()[1:], repeat=int(n-1)):

            cols = vector(F, [F.one()] + list(cols))

            if C == GeneralizedReedSolomonCode(alpha, dimension=C.dimension(), column_multipliers=cols):
                return True
        
    return False

def are_equivalent(C, D):

    return LinearCodeAutGroupCanLabel(C).get_canonical_form() == LinearCodeAutGroupCanLabel(D).get_canonical_form()

def Perm(C):

    n = C.length()
    CF = LinearCodeAutGroupCanLabel(C)
    Sn = SymmetricGroup(n)
    G = Sn.subgroup([s.get_perm() for s in CF.get_autom_gens()])

    return G

def ExtendedReedSolomonCode(F, dimension):

    q = F.order()
    k = dimension

    G = matrix(F, [[F.list()[i]^j if i != q else 0 if j != k - 1 else 1 for i in range(q + 1)] for j in range(k)])

    return LinearCode(G)

def FullLengthReedSolomonCode(F, dimension):

    return GeneralizedReedSolomonCode(F.list(), dimension=dimension)

def punctured_code(C, subset):

    subset = [s - 1 for s in subset]

    return LinearCode(C.generator_matrix()[:, subset])