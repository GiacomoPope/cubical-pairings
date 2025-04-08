from sage.all import EllipticCurve, GF, random_prime, randint, gcd

from cubical_pairings import cubical_weil_pairing, cubical_tate_pairing
from cubical_arithmetic import cubical_dbl, cubical_diff_add, cubical_ladder

TEST_CURVES = 10
TEST_ITER = 1000


def random_prime_3_mod_4():
    """
    Helper function to find a small prime such that p = 3 mod 4
    """
    p = random_prime(200)
    while p % 4 != 3 or p == 3:
        p = random_prime(200)
    return p


def random_supersingular_curve(F):
    """
    Helper function which takes a random step on the isogeny graph
    assumes that the torsion has only small factors and that p = 3 mod 4
    """
    # Create a "random" supersingular elliptic curve E / Fq in Montgomery form
    E = EllipticCurve(F, [0, 6, 0, 1, 0])
    E_prime = E.isogeny(E.random_point(), algorithm="factored").codomain()

    # Find a curve with a valid Montgomery model in the base field
    while True:
        try:
            return E_prime.montgomery_model()
        except ValueError:
            E_prime = E.isogeny(E.random_point(), algorithm="factored").codomain()


def random_fp_curve():
    """
    Helper function which computes a random supersingular elliptic curve
    in Montgomery form over a random field E / Fp
    """
    # Create a finite field Fp
    p = random_prime_3_mod_4()
    F = GF(p)

    # Create a "random" supersingular elliptic curve E / Fp in Montgomery form
    return random_supersingular_curve(F)


def random_fp2_curve():
    """
    Helper function which computes a random supersingular elliptic curve
    in Montgomery form over a random field E / Fp^2
    """
    # Create a finite field Fp^2
    p = random_prime_3_mod_4()
    F = GF((p, 2), names="i", modulus=[1, 0, 1])

    # Create a "random" supersingular elliptic curve E / Fp^2 in Montgomery form
    return random_supersingular_curve(F)


def random_curve():
    """
    Helper function which computes a random supersingular elliptic curve
    in Montgomery form over a random field E / Fp^k for k in {1, 2}
    """
    # Test E/Fp and E/Fp^2
    if randint(0, 1):
        return random_fp_curve()
    return random_fp2_curve()


def random_distinct_points(E):
    """
    Helper function which finds P, Q which are distinct and not zero
    """
    # Generate random distinct points which are not the identity
    # and not points of order two (to avoid edge cases)
    P = E.random_point()
    Q = E.random_point()

    # Ensure that P, Q, P + Q and P - Q are not the identity or the
    # point (0 : 1) on the Kummer line
    while P.is_zero() or Q.is_zero() or P == Q:
        P = E.random_point()
        Q = E.random_point()

    return P, Q


def random_torsion_points(E):
    """
    Helper function to compute two points in E[ell]. We require that ell
    divides p^2 - 1 for the reduced tate pairing.
    """
    while True:
        P = E.random_point()
        Q = E.random_point()
        ell = gcd(P.order(), Q.order())
        if ell > 1:
            break

    # All curves are supersingular, so the embedding degree should be 2
    p = E.base_ring().characteristic()
    assert ell.divides(p**2 - 1)

    # Clear the cofactor to ensure P, Q in E[ell]
    P = (P.order() / ell) * P
    Q = (Q.order() / ell) * Q

    return P, Q, ell


def test_point_equality(P, X, Z):
    """
    Helper function which checks if a x-only point (X : Z) is equal
    to a Sage type projective point (X : Y : Z)
    """
    if P.is_zero():
        return Z.is_zero()
    return P.x() == X / Z


def test_cubical_dbl(E):
    """
    Compare output of the cubical doubling against sage output
    """
    # Generate random points and their sum
    A24 = (E.a2() + 2) / 4
    P = E.random_point()
    while P.is_zero():
        P = E.random_point()
    PP = P + P

    # Convert to x-only representation
    P_tilde = (P.x(), 1)

    # Use cubical arithmetic to compute (XPQ : ZPQ)
    X2, Z2 = cubical_dbl(P_tilde, A24)

    # Ensure that the double is correct
    assert test_point_equality(PP, X2, Z2)


def test_cubical_diff_add(E):
    """
    Compare output of the cubical diff add against sage output
    """
    # Generate random distinct points where x(P - Q) is invertible
    P, Q = random_distinct_points(E)
    while (P - Q).x().is_zero():
        P, Q = random_distinct_points(E)
    PQ = P + Q

    # Convert to x-only representation
    P_tilde = (P.x(), 1)
    Q_tilde = (Q.x(), 1)
    ixPQ = 1 / (P - Q).x()

    # Use cubical arithmetic to compute (XPQ : ZPQ)
    XPQ, ZPQ = cubical_diff_add(P_tilde, Q_tilde, ixPQ)

    # Ensure that the sum is correct
    assert test_point_equality(PQ, XPQ, ZPQ)


def test_cubical_ladder(E):
    """
    Compare output of the cubical ladder against sage output
    """
    # Generate random distinct points with invertible x(P), x(Q), x(P - Q)
    # which is required for the ladder
    P, Q = random_distinct_points(E)
    while any(X.is_zero() or X.x().is_zero() for X in [P, Q, P - Q]):
        P, Q = random_distinct_points(E)

    # Generate a random scalar smaller than point order
    n_max = min(P.order(), Q.order())
    n = randint(1, n_max - 1)
    nP = n * P
    QnP = nP + Q

    # Compute the cubical ladder to obtain [n]P and [n]P + [n]
    (X1, Z1), (X2, Z2) = cubical_ladder(P, Q, n)

    assert test_point_equality(nP, X1, Z1)
    assert test_point_equality(QnP, X2, Z2)


def test_self_cubical_ladder(E):
    """
    Compare output of the cubical ladder against sage output
    """
    # Generate random distinct points with invertible x(P), x(Q), x(P - Q)
    # which is required for the ladder
    P = E.random_point()
    while P.is_zero() or P.x().is_zero():
        P = E.random_point()

    # Generate a random scalar smaller than point order
    n_max = P.order()
    n = randint(1, n_max - 1)
    nP = n * P
    nPP = nP + P

    # Compute the cubical ladder to obtain [n]P and [n + 1]P
    (X1, Z1), (X2, Z2) = cubical_ladder(P, P, n)

    assert test_point_equality(nP, X1, Z1)
    assert test_point_equality(nPP, X2, Z2)


def test_weil_pairing(E):
    """
    Compare output of the cubical weil pairing against sage output
    """
    # Compute a basis <P, Q> of E[ell] compatible with the algorithm
    P, Q, ell = random_torsion_points(E)

    # Compute the pairing with cubical arithmetic
    ePQ = cubical_weil_pairing(P, Q, ell, check=True)

    # For the odd case we compute the squared pairing
    if ell % 2:
        assert ePQ == P.weil_pairing(Q, ell) ** 2
    else:
        assert ePQ == P.weil_pairing(Q, ell)


def test_tate_pairing(E):
    """
    Compare output of the cubical tate pairing against sage output
    """
    # For the Tate pairing we need P in E[ell] and Q any other point
    P, _, ell = random_torsion_points(E)
    Q = E.random_point()

    # Compute the pairing with cubical arithmetic
    ePQ = cubical_tate_pairing(P, Q, ell, k=2, check=True)

    # For the odd case we compute the squared pairing
    # We use supersingular elliptic curves, so the embedding
    # degree is always 2
    if ell % 2:
        assert ePQ == P.tate_pairing(Q, ell, 2) ** 2
    else:
        assert ePQ == P.tate_pairing(Q, ell, 2)


def test_self_tate_pairing(E):
    """
    Compare output of the cubical tate self-pairing against sage output
    """
    # For the Tate pairing we need P in E[ell] and Q any other point
    P, _, ell = random_torsion_points(E)

    # Compute the pairing with cubical arithmetic
    ePQ = cubical_tate_pairing(P, P, ell, k=2, check=True)

    # For the odd case we compute the squared pairing
    # We use supersingular elliptic curves, so the embedding
    # degree is always 2
    if ell % 2:
        assert ePQ == P.tate_pairing(P, ell, 2) ** 2
    else:
        assert ePQ == P.tate_pairing(P, ell, 2)


def test_tate_pairing_T(E):
    # For the Tate pairing we need P in E[ell] and Q any other point
    P, _, ell = random_torsion_points(E)
    T = E((0, 0, 1))

    # Compute the pairing with cubical arithmetic
    ePT = cubical_tate_pairing(P, T, ell, k=2, check=True)
    if ell % 2:
        eTP = cubical_tate_pairing(T, P, 2 * ell, k=2, check=True)
    else:
        eTP = cubical_tate_pairing(T, P, ell, k=2, check=True)

    # For the odd case we compute the squared pairing
    # We use supersingular elliptic curves, so the embedding
    # degree is always 2
    if ell % 2:
        assert ePT == P.tate_pairing(T, ell, 2) ** 2
        assert eTP == T.tate_pairing(P, 2 * ell, 2) ** 2
        assert ePT == eTP == 1
    else:
        assert ePT == P.tate_pairing(T, ell, 2)
        assert eTP == T.tate_pairing(P, ell, 2)
        assert ePT**2 == ePT**2 == 1


if __name__ == "__main__":
    for _ in range(TEST_CURVES):
        E = random_curve()
        print(f"Testing curve... {E = }")
        for _ in range(TEST_ITER):
            # Arithmetic tests
            test_cubical_dbl(E)
            test_cubical_diff_add(E)
            test_cubical_ladder(E)
            test_self_cubical_ladder(E)

            # Pairing tests
            test_tate_pairing_T(E)
            test_self_tate_pairing(E)
            test_weil_pairing(E)
            test_tate_pairing(E)
    print("All tests passed")
