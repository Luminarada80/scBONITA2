# Rules for three incoming nodes
def A_AND_B_AND_C(A, B, C, not_a, not_b, not_c):
    A, B, C = (not A if not_a else A), (not B if not_b else B), (not C if not_c else C)
    return A and B and C

def A_AND_B_OR_C(A, B, C, not_a, not_b, not_c):
    A, B, C = (not A if not_a else A), (not B if not_b else B), (not C if not_c else C)
    return A and B or C

def A_OR_B_OR_C(A, B, C, not_a, not_b, not_c):
    A, B, C = (not A if not_a else A), (not B if not_b else B), (not C if not_c else C)
    return A or B or C

# Two incoming nodes
def A_AND_B(A, B, not_a, not_b):
    A, B = (not A if not_a else A), (not B if not_b else B)
    return A and B

def A_OR_B(A, B, not_a, not_b):
    A, B = (not A if not_a else A), (not B if not_b else B)
    return A or B

def A(A, not_a):
    A = not A if not_a else A
    return A