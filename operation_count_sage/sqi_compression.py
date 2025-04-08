from sage.all import ZZ, GF, EllipticCurve, randint, proof
from utilities import (
    batched_inversion,
    random_supersingular_montgomery,
    torsion_basis,
    dlog_2pt_tate,
    jac_ADD
)
from jacobian_pairing import dlog_2pt_tate_clz
from cubical_arithmetic import dlog_2pt_tate_biext

proof.all(False)

ITER_COUNT = 10

##################################
# SCRIPT CONSTANTS
##################################

# M = multiplications in Fp2
# S = squarings in Fp2
# A = additions in Fp2
# I = inversions in Fp2
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
    "M": 1,
    "S": 0.8,
    "A": 0.15,
    "I": 47,
    "sqrt": 222,
    "Mm": 0,
    "a": 0,
    "m": 0,
    "s": 0,
    "i": 0,  # all these should be ignored, everything is over Fp2
}
print()
print("cost model (relative to one Fp multiplication)".center(70))
print("-" * 70)
print(" | ".join((key.center(5)) for key in cost_model))
print(" | ".join((str(cost).center(5)) for cost in cost_model.values()))
print()


def reset_cost():
    return {key: 0 for key in cost_model}


def average_cost(cost_count, runs):
    return {key: cost_count[key]*1.0/runs for key in cost_count}

def compute_fp2_cost(cost_count, cost_model):
    fp2_cost = {key: 0 for key in cost_model}
    fp2_cost["M"] = sum([cost_count[key] * cost_model[key] for key in {"M", "I", "sqrt", "Mm", "a", "m", "s", "i"}])
    for key in {"S", "A"}:
        fp2_cost[key] = cost_count[key]
    return fp2_cost

def compute_cost(cost_count, cost_model):
    assert cost_count["Mm"] == 0
    assert cost_count["m"] == 0
    assert cost_count["s"] == 0
    assert cost_count["a"] == 0
    assert cost_count["i"] == 0
    return int(round(sum([cost_count[key] * cost_model[key] for key in cost_count])))


# Assumption: a global cost_count dictionary is used to keep track of arithmetic costs
# we initialize it as
cost_count = reset_cost()

##################################
# Parameter generation: SQIsign
##################################

# level 1: p = 5 * 2**248 - 1
# level 3:  p = 65 * 2**376 - 1
# level 5:  p = 27 * 2**500 - 1
SQISIGN_PARAMETERS = {
    1: (248, 5),
    3: (376, 65),
    5: (500, 27),
}


def gen_params(level):
    """
    Generate the field and parameters for SQIsign given a security
    level 1, 3 or 5
    """
    try:
        e, cofactor = SQISIGN_PARAMETERS[level]
    except KeyError:
        raise ValueError("The security level must be in [1, 3, 5]")

    p = ZZ(cofactor * 2**e - 1)
    assert p.is_prime()
    assert p % 4 == 3

    Fp2 = GF(p**2, "i", modulus=[1, 0, 1])
    return Fp2, e, cofactor


def gen_curve_and_bases(Fp2, e, cofactor):
    """
    Generate a random supersingular curve E with
    torsion basis E[2^e]
    """
    E = random_supersingular_montgomery(Fp2)

    # generate torsion basis
    P, Q = torsion_basis(E, 2**e, cofactor)
    PQ = P - Q

    return E, P, Q, PQ


##################################
# SQIsign point compression
##################################

# Setting:
# We're given
# - a Montgomery curve E,
# - a 2^e-torsion basis in x-only coordinates P=[XP,ZP], Q=[XQ,ZQ], PQ=P-Q=[XPQ,ZPQ],
# - a 2^f-torsion basis in x-only coordinates R=[XR,ZR], S=[XS,ZS], RS=R-S=[XRS,ZRS].
# Usually these bases are computed deterministically and come equipped with the difference point.
##
# Goal:
# Compute coefficients r1, r2, s1, s2 such that
# R = 2^(e-f) * (r1*P + r2*Q), S = 2^(e-f) * (s1*P + s2*Q)
# using pairings

for level in [1, 3, 5]:
    # generate parameters
    Fp2, e, cofactor = gen_params(level)
    E, P, Q, PQ = gen_curve_and_bases(Fp2, e, cofactor)
    A = E.a_invariants()[1]
    A24 = (A + 2) / 4

    # Init the cost counting dictionary
    cost_count = reset_cost()

    for iter in range(ITER_COUNT):
        # Pick a random f for some partial E[2^f] torsion to check
        f = randint(
            e // 2 + (iter - 1) * (e // (2 * ITER_COUNT)),
            e // 2 + iter * (e // (2 * ITER_COUNT)),
        )

        # for benchmarking, we fix this to f = e/2
        f = e // 2

        # Generate a random E[2^f] = <R, S>
        R, S = torsion_basis(E, 2**f, cofactor * 2 ** (e - f))

        # We work with x-only coordinates
        XP, ZP = P.x(), Fp2.one()
        XQ, ZQ = Q.x(), Fp2.one()
        XPQ, ZPQ = PQ.x(), Fp2.one()
        XR, ZR = R.x(), Fp2.one()
        XS, ZS = S.x(), Fp2.one()

        # 1. normalize points (we don't assume P, Q are normalized at the beginning)
        invs = batched_inversion([XP, ZP, XQ, ZPQ, ZR, ZS], cost_count, extdeg=2)
        # 15M 1I

        # Compute normalised x-coordinates
        xP = XP * invs[1]
        ixP = ZP * invs[0]
        ixQ = ZQ * invs[2]
        xPQ = XPQ * invs[3]
        xR = XR * invs[4]
        xS = XS * invs[5]
        cost_count["M"] += 6

        # Package points as (X : Z) projective x-only
        Px = (xP, Fp2.one())
        PQx = (xPQ, Fp2.one())
        Rx = (xR, Fp2.one())
        Sx = (xS, Fp2.one())

        # Package up basis information
        PQ_basis = (Px, PQx, ixP, ixQ, e)
        RS_basis = (Rx, Sx, f)

        # 2. compute pairwise differences
        RP = P - R
        RQ = R - Q
        SP = P - S
        SQ = S - Q

        # Package points as (X : Z) projective x-only
        RPx = (RP.x(), Fp2.one())
        RQx = (RQ.x(), Fp2.one())
        SPx = (SP.x(), Fp2.one())
        SQx = (SQ.x(), Fp2.one())

        diff_points = (RPx, RQx, SPx, SQx)

        # in SQIsign points are given in x-only form,
        # but we want P, Q in modified Jacobian coordinates:
        # we need to lift bases (P,Q) and (R,S)
        # cost per basis:
        # # recover one y(P): 1S 2M 2A 1sqrt for sqrt(x^3 + Ax^2 + 1)
        # # recover y(Q) from (x(P), x(Q), x(P-Q), y(P)): 1S 5M 8A
        # # # see https://github.com/LinKaizhan/Pairingoptimizations/blob/main/verify_sup.sage with all Z coordinates = 1
        # lift_basis: 2S 7M 10A 1sqrt
        ##
        # Jacobian addition with Z = 1:
        # cost per addition 2M + 3S + 8A
        ##
        # total cost: 2*lift_basis + 4*add = 46M + 12S + 2sqrt + 72A
        cost_count["M"] += 22
        cost_count["S"] += 16
        cost_count["A"] += 52
        cost_count["sqrt"] += 2

        # 3. compute pairings and dlogs
        r1, r2, s1, s2 = dlog_2pt_tate_biext(
            PQ_basis, RS_basis, diff_points, A24, e, f, cofactor, cost_count
        )
        assert 2 ** (e - f) * (r1 * P + r2 * Q) == R
        assert 2 ** (e - f) * (s1 * P + s2 * Q) == S

        expected_coeffs = dlog_2pt_tate(P, Q, R, S, e, f)
        assert expected_coeffs == (r1, r2, s1, s2)

    cost_count = average_cost(cost_count,ITER_COUNT)
    print(f"Level {level}, {f = }: {compute_fp2_cost(cost_count, cost_model) = }")
    print(
        f"Cost in Fp2-multiplications:  \033[1m {compute_cost(cost_count, cost_model)} \033[0m"
    )
    print("All tests passed!")

print("----CLZ--------------------------")

for level in [1, 3, 5]:

    # generate parameters
    Fp2, e, cofactor = gen_params(level)
    E, P, Q, PQ = gen_curve_and_bases(Fp2, e, cofactor)
    A = E.a_invariants()[1]

    cost_count = reset_cost()

    # for benchmarking, we fix this to f = e/2
    f = e // 2
    R, S = torsion_basis(E, 2**f, cofactor * 2 ** (e - f))
    RS = R - S

    # in SQIsign points are given in x-only form,
    # but we want P, Q in modified Jacobian coordinates:
    # we need to lift bases (P,Q) and (R,S)
    # cost per basis:
    # # recover one y(P): 1S 2M 2A 1sqrt for sqrt(x^3 + Ax^2 + 1)
    # # recover y(Q) from (x(P), x(Q), x(P-Q), y(P)): 1S 5M 8A
    # # # see https://github.com/LinKaizhan/Pairingoptimizations/blob/main/verify_sup.sage with all Z coordinates = 1
    # total cost: 14M + 4S + 10A + 2sqrt
    cost_count["M"] += 14
    cost_count["S"] += 4
    cost_count["A"] += 20
    cost_count["sqrt"] += 2

    ## assume A/3 is precomputed
    AA = Fp2(1 - A**2 / 3)
    BB = Fp2(A**3 / 27 * 2 - A / 3)
    EE = EllipticCurve(Fp2, [AA, BB])
    Pp = EE([P.x() + A / 3, P.y()])
    Qp = EE([Q.x() + A / 3, Q.y()])
    Rp = EE([R.x() + A / 3, R.y()])
    Sp = EE([S.x() + A / 3, S.y()])
    cost_count["A"] += 6

    PQ_basis = Pp, Qp, e
    RS_basis = Rp, Sp, f

    # 3. compute pairings and dlogs
    r1, r2, s1, s2 = dlog_2pt_tate_clz(
        PQ_basis, RS_basis, AA, e, f, cofactor, cost_count
    )
    assert 2 ** (e - f) * (r1 * P + r2 * Q) == R
    assert 2 ** (e - f) * (s1 * P + s2 * Q) == S

    expected_coeffs = dlog_2pt_tate(P, Q, R, S, e, f)
    assert expected_coeffs == (r1, r2, s1, s2)

    print(f"Level {level}, {f = }: {compute_fp2_cost(cost_count, cost_model) = }")
    print(
        f"Cost in Fp2-multiplications: \033[1m {compute_cost(cost_count, cost_model)} \033[0m"
    )
