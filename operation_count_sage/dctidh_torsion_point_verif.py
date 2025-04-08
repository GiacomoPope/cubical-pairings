from sage.all import ZZ, EllipticCurve, GF, prod, floor, log, randrange, proof, factor
from cubical_arithmetic import tate_pairing_fp, ladder, ladder_proj
from utilities import batched_inversion, lucas_sequences, batch_cof_trace_total, lucas_squaring, seeded_elligator
from jacobian_pairing import miller_modified_odd_naf_trace_fp

# Precomputed public keys for dCTIDH
from dctidh_public_keys import public_keys

proof.all(False)

##################################
# SCRIPT CONSTANTS
##################################

NUM_MONT_CURVES = len(public_keys)

# M = multiplications in Fp2
# S = squarings in Fp2
# A = additions in Fp2
# Mm = mixed multiplications Fp2 - Fp
# m = multiplications in Fp
# s = squarings in Fp
# a = additions in Fp2
# i = inversions in Fp
# sqrt = square roots in Fp


cost_model = {
    "M": 3,
    "S": 2,
    "A": 0.00,
    "Mm": 2,
    "m": 1,
    "s": 1,
    "a": 0.00,
    "i": 57,
    "sqrt": 2150,
}
print()
print("cost model (relative to one Fp multiplication)".center(70))
print("-" * 70)
print(" | ".join((key.center(5)) for key in cost_model))
print(" | ".join((str(cost).center(5)) for cost in cost_model.values()))
print()


def reset_cost():
    return {key: 0 for key in cost_model}


def compute_fp_cost(cost_count, cost_model):
    fp_cost = {key: 0 for key in cost_model}
    fp_cost["m"] = sum([cost_count[key] * cost_model[key] for key in {"M", "S", "Mm", "i", "sqrt", "m"}])
    fp_cost["s"] = cost_count["s"]
    fp_cost["a"] = cost_count["a"] + 2*cost_count["A"]
    return fp_cost

def compute_cost(cost_count, cost_model):
    return int(round(sum([cost_count[key] * cost_model[key] for key in cost_count])))


# Assumption: a global cost_count dictionary is used to keep track of arithmetic costs
# we initialize it as
cost_count = reset_cost()

##################################
# Parameter generation
##################################

n = 194
ells = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187
]

ellsg = [
    3, 5, 7*7, 11, 13, 17, 19, 23, 29, 31, 37, 41*41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187
]

assert len(ells) == n
L = prod(ells)
f = 64*6 + 3    # 2-torsion
g = 7*41        # cofactor, important to keep track of

# Compute fields
p = 2**f * g * L - 1 #dCTIDH 2048-194 prime
Fp = GF(p)
Fp2 = GF(p**2, "i", modulus=[1, 0, 1])


##################################
# Torsion point verification
# Algorithm 4 in eprint 2023/858
# also presented as algorithm eprint 2024/575
##################################

# We use the same algorithm with the following changes:
# - we use x-only Montgomery arithmetic
# - we use biextension pairings


def test_torsion_basis(A, u):
    # A is the A-coefficient of a Montgomery curve
    A = Fp(A)
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])
    [P, Q] = seeded_elligator(A, u, E, Fp, Fp2, cost_count)

    # 1. clear cofactor
    Px = [P.x(), 1]
    nPx = ladder_proj(Px, A24, 2**f, cost_count, extdegs=[1, 1])

    # renormalize Px and recover its y coordinate
    nP = E.lift_x(nPx[0] / nPx[1])
    # in 2024/575 the ladder also returns (X,Z) coords of [n+1]P (for free),
    # then recovers y([n]P) from (X,Y,Z)(P), (X,Z)([n]P), (X,Z)([n+1]P)
    # the following cost count refers to this algorithm
    cost_count["m"] += 11
    cost_count["s"] += 1
    cost_count["a"] += 8

    Px = [nP.x(), 1]
    Qx = [Q.x(), 1]
    
    # derive right PQ coordinate by jacobian addition
    PQx = [(nP - Q).x(), 1]
    cost_count["m"] += 2
    cost_count["s"] += 3
    cost_count["A"] += 1
    cost_count["a"] += 7

    # batch compute inverses
    xPQconj = PQx[0].conjugate()
    ixP, ixQ, ixQnorm = batched_inversion(
        [Px[0], Qx[0], PQx[0].norm()], cost_count, extdeg=1
    )
    ixPQ = xPQconj * ixQnorm
    cost_count["Mm"] += 1
    cost_count["s"] += 2


    # 2. compute order-N pairing
    try:
        e = tate_pairing_fp(Px, Qx, PQx, ixP, ixQ, ixPQ, A24, L*g, cost_count)
    except ValueError:  # P is not N-torsion ->
        # the original point is not of (p+1)-torsion, hence E is not supersingular
        return False

    # 3. compute (the trace of) the final exponentiation

    # first get rid of the Fp^* part
    e = e.conjugate()/e
    # computed as e.conjugate^2 / norm(e)
    # but we only need the real part anyway, so we compute
    # e = (a^2 - b^2) / (a^2 + b^2)
    cost_count["s"] += 2
    cost_count["a"] += 2
    cost_count["i"] += 1
    cost_count["m"] += 1

    e = e + e.conjugate()
    # this really just means to get 2*Re(e) so its free
    # then get rid of the remaining exponentiation where we can now use lucas seqs
    # to compute e^(2^f) where 2^f == (p+1) // L*g
    e = lucas_squaring(e, f, cost_count)        

    # 4. check order of e.
    # This batch order verification is the same used in 2024/575
    # for more info, compare with https://github.com/LearningToSQI/SQISign-SageMath/blob/1217b8fd68dc750364ab5096cb14bce6e2df1d1d/utilities.py#L267C1-L322C16
    accs = [1, 0]
    ind = [i for i in range(len(ellsg))]
    batch_cof_trace_total(accs, e, ind, len(ellsg), ellsg, p, cost_count)
    return accs[1]


# Same algorithm, but now compute pairing using Miller's algorithm as in 2024/575
def test_torsion_basis_miller(A, u):
    # A is the A-coefficient of a Montgomery curve
    A = Fp(A)
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])

    [P, Q] = seeded_elligator(A, u, E, Fp, Fp2, cost_count)

    # 1. clear cofactor
    Px = [P.x(), 1]
    nPx = ladder_proj(Px, A24, 2**f, cost_count, extdegs=[1, 1])

    # renormalize Px and recover its y coordinate
    nP = E.lift_x(nPx[0] / nPx[1])
    # in 2024/575 the ladder also returns (X,Z) coords of [n+1]P (for free),
    # then recovers y([n]P) from (X,Y,Z)(P), (X,Z)([n]P), (X,Z)([n+1]P)
    # the following cost count refers to this algorithm
    cost_count["m"] += 11
    cost_count["s"] += 1
    cost_count["a"] += 8

    # compute an isomorphic Montgomery model via a translation x -> x + A/3
    AA = Fp2(1 - A**2 / 3)
    BB = Fp2(A**3 / 27 * 2 - A / 3)
    EE = EllipticCurve(Fp2, [AA, BB])

    Pp = EE([nP.x() + A / 3, nP.y()])
    Qp = EE([Q.x() + A / 3, Q.y()])
    cost_count["a"] += 2

    XP, YP = Pp.xy()
    XQ, YQ = Qp.xy()

    # 2. compute the order-N pairing
    [e, label] = miller_modified_odd_naf_trace_fp(
        XQ, YQ, XP, YP, 1, AA, AA, L*g, L*g, p, cost_count
    )

    if label == 0:
        return False

    # 3. compute (the trace of) the final exponentiation

    # first get rid of the Fp^* part
    e = e.conjugate()/e
    # computed as e.conjugate^2 / norm(e)
    # but we only need the real part anyway, so we compute
    # e = (a^2 - b^2) / (a^2 + b^2)
    cost_count["s"] += 2
    cost_count["a"] += 2
    cost_count["i"] += 1
    cost_count["m"] += 1

    e = e + e.conjugate()
    # this really just means to get 2*Re(e) so its free
    # then get rid of the remaining exponentiation where we can now use lucas seqs
    # to compute e^(2^f) where 2^f == (p+1) // L*g
    e = lucas_squaring(e, f, cost_count) 


    # 4. check order of e.
    # This batch order verification uses the same approach as 2024/575
    # for more info, compare with https://github.com/LearningToSQI/SQISign-SageMath/blob/1217b8fd68dc750364ab5096cb14bce6e2df1d1d/utilities.py#L267C1-L322C16
    accs = [1, 0]
    ind = [i for i in range(len(ellsg))]
    batch_cof_trace_total(accs, e, ind, len(ellsg), ellsg, p, cost_count)
    return accs[1]


# ==================== #
#       Testing        #
# ==================== #

function_info = [
    ("Rei23 algorithm 4", test_torsion_basis),
    ("Rei23 algorithm 4-CLZ", test_torsion_basis_miller),
]

for name, h in function_info:
    for A, u in public_keys:
        if not h(A, u):
            print(f"{name} has a false negative for {A = }")
        if h(A,u - 1):
            print(f"{name} has a false positive for {u - 1 = }")
        if h(A + 1,u):
            print(f"{name} has a false positive for {A + 1 = }")

# ==================== #
#     Benchmarking     #
# ==================== #

cost_count = reset_cost()

for name, h in function_info:
    total_costs = []
    for A, u in public_keys:
        cost_count = reset_cost()
        h(A, u)
        total_costs.append(compute_cost(cost_count, cost_model))

    print(cost_count)
    print(compute_cost(cost_count, cost_model))
    print(compute_fp_cost(cost_count, cost_model))
    print(f"Cost of {name} in Fp-multiplications m")
    print(f"Min cost: \033[1m {min(total_costs)} \033[0m")
    print(f"Max cost: \033[1m {max(total_costs)} \033[0m")
    print(f"Avg cost: \033[1m {sum(total_costs) / NUM_MONT_CURVES} \033[0m")
