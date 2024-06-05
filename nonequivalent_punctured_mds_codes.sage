load("star_products.sage")

F.<a> = GF(8)

C = LinearCode(matrix(F, [[a, a^2, 1, a^2 + a, a^2 + a + 1, 1, a^2 + a], [0, a^2 + 1, a^2 + 1, a + 1, a^2 + a + 1, 1, a], [a + 1, 0, a^2 + 1, 0, a^2 + a, 1, a^2 + 1]]))

# These two puncturings of C happen to not be equivalent, these are [6, 3]_8 MDS codes
A = C.punctured([1])
B = C.punctured([4])

print("A and B are equivalent:", are_equivalent(A, B))

assert is_mds(C) and is_mds(A) and is_mds(B)

l = 2

sum_A = sum_dimensions(A, l)  # prints 112664981

sum_B = sum_dimensions(B, l)  # prints 112782630

# These are different even though both come from the same large MDS code