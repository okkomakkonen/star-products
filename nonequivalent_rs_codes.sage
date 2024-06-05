load("star_products.sage")

A = LinearCode(matrix(GF(11), [[1, 0, 0, 6, 4, 6], [0, 1, 0, 10, 10, 3], [0, 0, 1, 7, 9, 3]]))
B = LinearCode(matrix(GF(11), [[1, 0, 0, 2, 6, 1], [0, 1, 0, 2, 7, 4], [0, 0, 1, 8, 10, 7]]))
C = LinearCode(matrix(GF(11), [[1, 0, 0, 10, 8, 8], [0, 1, 0, 6, 8, 2], [0, 0, 1, 7, 7, 2]]))
D = LinearCode(matrix(GF(11), [[1, 0, 0, 6, 1, 8], [0, 1, 0, 7, 4, 7], [0, 0, 1, 10, 7, 8]]))

assert is_mds(A) and is_mds(B) and is_mds(C) and is_mds(D)

print("A pair of these codes are equivalent:", are_equivalent(A, B) or are_equivalent(A, C) or are_equivalent(A, D) or are_equivalent(B, C) or are_equivalent(B, D) or are_equivalent(C, D))

l = 2

sum_dimensions(A, l)  # prints 1393999830

sum_dimensions(B, l)  # prints 1393999830

sum_dimensions(C, l)  # prints 1393999830

sum_dimensions(D, l)  # prints 1393999830

# These are all the same even though these codes are pairwise nonequivalent