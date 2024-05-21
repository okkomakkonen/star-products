load("star_products.sage")

from sage.coding.grs_code import *

alpha = vector(GF(7), [0, 1, 2, 3, 4, 5])
beta  = vector(GF(7), [1, 2, 3, 4, 5, 6])

A = GeneralizedReedSolomonCode(alpha, dimension=3)
B = GeneralizedReedSolomonCode(beta, dimension=3)

assert is_mds(A) and is_mds(B)

l = 2

sum_A = sum_dimensions(A, l)

print(sum_A)  # prints 39462150

sum_B = sum_dimensions(B, l)

print(sum_B)  # prints 39462150

# These are equal