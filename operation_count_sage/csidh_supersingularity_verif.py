from sage.all import ZZ, EllipticCurve, GF, prod, floor, log, randrange, proof
from cubical_arithmetic import tate_pairing_fp, ladder, ladder_proj
from utilities import batched_inversion, elligator, lucas_sequences, batch_cof_trace_rec
from jacobian_pairing import miller_modified_odd_naf_trace_fp


proof.all(False)

##################################
# SCRIPT CONSTANTS
##################################

NUM_MONT_CURVES = 10

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
    "A": 0.16,
    "Mm": 2,
    "m": 1,
    "s": 0.67,
    "a": 0.08,
    "i": 64,
    "sqrt": 459,
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

# fmt: off
ells = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 
    163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 
    241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 
    337, 347, 349, 353, 359, 367, 373, 587
]
# fmt: on

# Compute fields
p = ZZ(4 * prod(ells) - 1)  # CSIDH-512 prime
Fp = GF(p)
Fp2 = GF(p**2, "i", modulus=[1, 0, 1])

# Compute cofactor values
ellsN = ells[41:74]
N = prod(ellsN)  # N >= 4 sqrt(p)
assert N >= 4 * p.sqrt()
cofactor = (p + 1) // N


# Toy implementation of csidh to generate random public keys
def csidh(pub, priv):
    """
    Toy CSIDH implementation

    Given a public key (Montgomery coefficient) and a private key
    (an array of integers), output the corresponding curve
    """
    # Compute Montgomery coefficient from public key
    E = EllipticCurve(Fp, [0, pub, 0, 1, 0])

    # Ensure that the curve is supersingular
    assert (p + 1) * E.random_point() == E(0)

    for es in ([max(0, +e) for e in priv], [max(0, -e) for e in priv]):
        while any(es):
            js = [j for j, e in enumerate(es) if e]

            k = prod(ells[j] for j in js)
            P = E.random_point()
            P *= (p + 1) // k
            if not P:
                continue

            for j in js[:]:
                k //= ells[j]
                Q = k * P

                if not Q:
                    continue

                Q.set_order(ells[j])
                phi = E.isogeny(Q, model="montgomery")
                E, P = phi.codomain(), phi(P)
                es[j] -= 1
        E = E.quadratic_twist()
        E = E.montgomery_model()

    A = E.a_invariants()[1]
    assert (0, A, 0, 1, 0) == E.a_invariants(), f"{E.a_invariants() = }"
    return A


# Precompute some "random" A-invariants
print(f"Generating {NUM_MONT_CURVES} supersingular curves...")
supersingular_As = [Fp(6)]
while len(supersingular_As) != NUM_MONT_CURVES:
    pub = supersingular_As[-1]
    priv = [randrange(-1, 1) for _ in ells]
    A = csidh(pub, priv)
    supersingular_As.append(A)


##################################
# Supersingularity verification
# algorithm eprint 2024/575
# also presented as Algorithm 5 in eprint 2023/858
##################################

# We use the same algorithm with the following changes:
# - we use x-only Montgomery arithmetic
# - we use biextension pairings


def test_supersingular_alg5(A):
    # A is the A-coefficient of a Montgomery curve
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])
    [P, Q] = elligator(A, E, Fp, Fp2, cost_count)

    # the Elligator used above gives the (x, y) coordinates of P and Q. 
    Px = [P.x(), 1]
    Qx = [Q.x(), 1]

    # by jacobian addition  2M + 3S + 8A
    PQx = [(P - Q).x(), 1]
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

    # 1. clear cofactor
    # cof = (p+1) / N, N > 4*sqrt(p)
    Px = ladder(Px, ixP, A24, cofactor, cost_count, extdegs=[1, 1])

    # renormalize Px setting zP = 1
    xP = Px[0]
    zP = Px[1]
    xzP = xP * zP
    ixzP = 1 / xzP
    Px = [xP * xP * ixzP, 1]
    ixP = zP * zP * ixzP
    cost_count["m"] += 3
    cost_count["s"] += 2
    cost_count["i"] += 1

    # 2. compute order-N pairing
    try:
        e = tate_pairing_fp(Px, Qx, PQx, ixP, ixQ, ixPQ, A24, N, cost_count)
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
    e = lucas_sequences(e, (p + 1) // N, cost_count)

    # 4. check order of e.
    # This batch order verification is the same used in 2024/575
    # for more info, compare with https://github.com/LearningToSQI/SQISign-SageMath/blob/1217b8fd68dc750364ab5096cb14bce6e2df1d1d/utilities.py#L267C1-L322C16
    accs = [1, 0]
    ind = [i for i in range(len(ellsN))]
    batch_cof_trace_rec(accs, e, ind, len(ellsN), ellsN, p, cost_count)
    if accs[1]:
        return True
    else:   
        # there was too little torsion, we need to run it again
        # this can be done more efficiently, we could reuse the torsion
        # already computed in accs[0] with more sophisticated techniques
        # but the gain is minimal, as the failure rate is low
        return test_supersingular_alg5(A)


# Same algorithm, but now compute pairing using Miller's algorithm as in 2024/575
def test_supersingular_alg5_miller(A):
    # A is the A-coefficient of a Montgomery curve
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])

    [P, Q] = elligator(A, E, Fp, Fp2, cost_count)

    # clear cofactor
    Px = [P.x(), 1]
    nPx = ladder_proj(Px, A24, cofactor, cost_count, extdegs=[1, 1])

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

    # compute the order-N pairing
    [e, label] = miller_modified_odd_naf_trace_fp(
        XQ, YQ, XP, YP, 1, AA, AA, N, N, p, cost_count
    )

    if label == 0:
        return False

    # get rid of the Fp^* part
    e = e.conjugate()/e
    # computed as e.conjugate^2 / norm(e)
    # but we only need the real part anyway, so we compute
    # e = (a^2 - b^2) / (a^2 + b^2)
    cost_count["s"] += 2
    cost_count["a"] += 2
    cost_count["i"] += 1
    cost_count["m"] += 1

    e = e + e.conjugate()
    e = lucas_sequences(e, (p + 1) // N, cost_count)

    # check order of e.
    # This batch order verification uses the same approach as 2024/575
    # for more info, compare with https://github.com/LearningToSQI/SQISign-SageMath/blob/1217b8fd68dc750364ab5096cb14bce6e2df1d1d/utilities.py#L267C1-L322C16
    accs = [1, 0]
    ind = [i for i in range(len(ellsN))]
    batch_cof_trace_rec(accs, e, ind, len(ellsN), ellsN, p, cost_count)
    if accs[1]:
        return True
    else:   
        # there was too little torsion, we need to run it again
        # this can be done more efficiently, we could reuse the torsion
        # already computed in accs[0] with more sophisticated techniques
        # but the gain is minimal, as the failure rate is low
        return test_supersingular_alg5(A)

##################################
# Supersingularity verification
# using Approach 1 - Algorithm 4 in eprint 2023/858
##################################

# instead of computing the order-N pairing between [cof]P and Q,
# we now compute the order-(p+1) pairing between P and Q
# and check that the result has order exactly p+1
# using the same batch order verification as above

# since p+1 is even, and we want to check it has full order,
# we need to use the algorithm for even pairing.
# otherwise we get e(P,Q)^2: we'd lose 1 bit of information about its order


def test_supersingular_alg4(A):
    # A is the A-coefficient of a Montgomery curve
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])

    [P, Q] = elligator(A, E, Fp, Fp2, cost_count)
    Px = [P.x(), 1]
    Qx = [Q.x(), 1]

    # derive right PQ coordinate by jacobian addition
    PQx = [(P - Q).x(), 1]
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

    # 2. compute order-(p+1) pairing
    try:
        e = tate_pairing_fp(Px, Qx, PQx, ixP, ixQ, ixPQ, A24, p + 1, cost_count)
    except ValueError:  # P is not (p+1)-torsion -> E is not supersingular
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

    # 4. check order of e.
    # This batch order verification uses the same approach as 2024/575
    # for more info, compare with https://github.com/LearningToSQI/SQISign-SageMath/blob/1217b8fd68dc750364ab5096cb14bce6e2df1d1d/utilities.py#L267C1-L322C16
    accs = [1, 0]
    ind = [i for i in range(len(ells))]
    batch_cof_trace_rec(accs, e, ind, len(ells), ells, p, cost_count)
    return accs[1]


##################################
# Supersingularity verification
# using Doliskani's test
##################################


# Doliskani's test reformulated as:
# E is (most probably) supersingular if, for a random P in E(Fp2), the reduced Tate pairing t(P,P) is 1
# i.e. the non-reduced self pairing e(P,P) belongs to Fp


def test_supersingular_doliskani_bgs(A):
    # A is the A-coefficient of a Montgomery curve
    A24 = (A + 2) / 4

    xP = Fp2.random_element()
    Px = (xP, 1)
    pP = ladder_proj(Px, A24, p, cost_count, extdegs=[2, 1])

    cost_count["M"] += 1
    if pP[0] != xP * pP[1]:
        return False

    q = Fp(4) * xP.conjugate()
    cost_count["A"] += 2
    for _ in range(floor(log(p, 2)) + 1):
        q = q**2
    cost_count["S"] += floor(log(p, 2) + 1)

    cost_count["M"] += 1
    cost_count["A"] += 2
    if Fp(4) * xP * pP[1] == q:
        return True
    return False


def test_supersingular_doliskani_pp1(A):
    # this is the same as the previous test, but we check the pairing with P+1
    # we also compute the ladder forgetting about factors 1/4 in diff_add
    # this ladder gives a non-reduced Tate pairing.
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])

    P = E.random_point()
    xP = P.x()
    Px = (xP, 1)
    ixP = 1 / xP
    # for an Fp2-inversion, compute xP.norm(), invert it in Fp,
    # then compute xP.conjugate() * inverse(xP.norm())
    cost_count["i"] += 1
    cost_count["s"] += 2
    cost_count["Mm"] += 1

    p1P = ladder(Px, ixP, A24, p + 1, cost_count, extdegs=[2, 2])
    if p1P[1] != 0:
        # P is not a (p+1)-torsion point
        return False
    if p1P[0].list()[1] == 0:  # X([p+1]P) lies in Fp
        # X([p+1]P) = e_{p +- 1}(P,P)^2 is the square of the NON-reduced Tate pairing.
        # this lies in Fp iff the square of the reduced Tate pairing is 1
        return True
    return False


def test_supersingular_doliskani_p(A):
    # this is the test proposed in Rob24
    A24 = (A + 2) / 4
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])

    P = E.random_point()
    xP = P.x()
    Px = (xP, 1)
    ixP = 1 / xP
    # for an Fp2-inversion, compute xP.norm() (cost=2s), invert it in Fp,
    # then compute xP.conjugate() * inverse(xP.norm()) (cost=1Mm)
    cost_count["i"] += 1
    cost_count["s"] += 2
    cost_count["Mm"] += 1

    pP = ladder(Px, ixP, A24, p, cost_count, extdegs=[2, 2], correct=True)
    if (pP[0] == xP) and (pP[1] == 1):
        # this is up to a (p-1)-th root of unity, i.e. up to any Fp element...
        return True
    return False


# Simple ladder test


def test_simple_ladder(A):
    A24 = (A + 2) / 4

    while True:
        cost_count["m"] += 2
        cost_count["s"] += 2
        cost_count["a"] += 2
        r = Fp.random_element()
        if r.is_zero():
            r = Fp.one()
        r = -r * r
        rhs = r**3 + A * r**2 + r
        if rhs.is_square():
            break

    P = [r, 1]
    Pn = ladder_proj(P, A24, p, cost_count, extdegs=[1, 1])

    cost_count["m"] += 2
    return P[0] * Pn[1] == Pn[0] * P[1]


# ==================== #
#       Testing        #
# ==================== #

function_info = [
    ("Rei23 algorithm 4", test_supersingular_alg4),
    ("Rei23 algorithm 5", test_supersingular_alg5),
    ("Rei23 algorithm 5-CLZ", test_supersingular_alg5_miller),
    ("Doliskani-BGS", test_supersingular_doliskani_bgs),
    ("Doliskani-biext-p", test_supersingular_doliskani_p),
    ("Doliskani-biext-(p+1)", test_supersingular_doliskani_pp1),
    ("Simple Ladder by (p+1)", test_simple_ladder),
]

for name, f in function_info:
    for A in supersingular_As:
        if not f(A):
            print(f"{name} has a false negative for {A = }")
        if f(A + 1):
            print(f"{name} has a false positive for {A + 1 = }")

# ==================== #
#     Benchmarking     #
# ==================== #

cost_count = reset_cost()

for name, f in function_info:
    total_costs = []
    for A in supersingular_As:
        cost_count = reset_cost()
        f(A)
        total_costs.append(compute_cost(cost_count, cost_model))

    # print(cost_count)
    print(f"\nCost of {name}")
    fp_cost = compute_fp_cost(cost_count, cost_model)
    for key in ["m", "s", "a"]:
        print(f" {round(fp_cost[key]*1.0)}", end=" &")
    print(f"\033[1m {sum(total_costs) / NUM_MONT_CURVES} \033[0m \\\\")    
    print(f"Avg cost: \033[1m {sum(total_costs) / NUM_MONT_CURVES} \033[0m")
