load("star_products.sage")

A = LinearCode(matrix(GF(7), [[1, 0, 0, 4, 5, 2], [0, 1, 0, 6, 1, 1], [0, 0, 1, 5, 6, 5]]))
B = LinearCode(matrix(GF(7), [[1, 0, 0, 1, 1, 6], [0, 1, 0, 4, 1, 4], [0, 0, 1, 6, 2, 4]]))

assert is_mds(A) and is_mds(B)

l = 2

sum_A = sum_dimensions(A, l)

print(sum_A)  # prints 39415494

sum_B = sum_dimensions(B, l)

print(sum_B)  # prints 39462150

# These will be different even though both codes are [6, 3] MDS codes