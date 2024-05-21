load("star_products.sage")

A = LinearCode(matrix(GF(11), [[1, 0, 0, 6, 4, 6], [0, 1, 0, 10, 10, 3], [0, 0, 1, 7, 9, 3]]))
B = LinearCode(matrix(GF(11), [[1, 0, 0, 2, 6, 1], [0, 1, 0, 2, 7, 4], [0, 0, 1, 8, 10, 7]]))

assert is_mds(A) and is_mds(B)

l = 2

sum_A = sum_dimensions(A, l)

print(sum_A)

sum_B = sum_dimensions(B, l)

print(sum_B)
