from sage.all import ZZ, EllipticCurve_from_j, supersingular_j, floor, prod, factor, discrete_log
from sage.groups.generic import order_from_multiple


##################################
# Utilities for testing
##################################

def jac_ADD(P, Q, A, cost_count, extdegs=[2, 2, 2]):
    """
    Algorithm 8.13 of the SQISign spec.

    extdeg indicates Fp or Fp2 for [P, Q, A]

    Given distinct points P, Q and the Montgomery coefficient A compute u, v, w
    such that x(P + Q) = (u - v : w) and x(P - Q) = (u + v : w)
    """
    # Unpack coordinates
    XP, YP, ZP = P
    XQ, YQ, ZQ = Q

    # Compute the slope dy / dx for the case P != Q
    t0 = ZP * ZP
    t1 = ZQ * ZQ
    t2 = t1 * XP
    t3 = t0 * XQ
    t4 = YP * ZQ
    t4 = t4 * t1
    t5 = ZP * YQ
    t5 = t5 * t0
    t0 = t0 * t1
    t6 = t4 * t5
    v = t6 + t6
    t4 = t4 * t4
    t5 = t5 * t5
    t4 = t4 + t5
    t5 = t2 + t3
    t6 = t3 + t3
    t6 = t5 - t6
    t6 = t6 * t6
    t1 = A * t0
    t1 = t5 + t1
    t1 = t1 * t6
    u = t4 - t1
    w = t6 * t0

    # Update the costs
    if extdegs == [1, 1, 1]:
        cost_count["m"] += 11
        cost_count["s"] += 5
        cost_count["a"] += 7
    elif (extdegs == [1, 2, 1]) or (extdegs == [2, 1, 1]):
        cost_count["M"] += 6
        cost_count["Mm"] += 5
        cost_count["S"] += 4
        cost_count["s"] += 1
        cost_count["A"] += 7
    elif (extdegs == [1, 2, 2]) or (extdegs == [2, 1, 2]):
        cost_count["M"] += 7
        cost_count["Mm"] += 4
        cost_count["S"] += 5
        cost_count["s"] += 1
        cost_count["A"] += 7
    elif extdegs == [2, 2, 2]:
        cost_count["M"] += 11
        cost_count["S"] += 5
        cost_count["A"] += 7

    return u - v, u + v, w

def torsion_basis(E, n, cofactor=1):
    """
    Compute a basis for the n-torsion of E.
    where E is a supersingular elliptic curve
    with rational n*cofactor torsion.
    """

    def point_order_n(E, n, cofactor):
        """
        Compute a random point P of order n
        """
        while True:
            P = cofactor * E.random_point()
            P.set_order(multiple=n)
            if P.order() == n:
                return P

    # We can pick any point P of order n and then find some
    # Q such that e(P, Q) has full order
    P = point_order_n(E, n, cofactor)
    while True:
        Q = point_order_n(E, n, cofactor)
        ePQ = P.weil_pairing(Q, n)
        o = order_from_multiple(ePQ, n, operation="*")
        if o == n:
            return P, Q


def random_supersingular_montgomery(F):
    """
    Generate a "random" supersingular elliptic curve
    E / F for a given field Fq in Montgomery form.

    Assumes that p+1 is smooth
    """
    j = supersingular_j(F)
    E0 = EllipticCurve_from_j(j)
    ker = E0.random_point()
    return E0.isogeny(ker, algorithm="factored", model="montgomery").codomain()


def batched_inversion(input, cost_count, extdeg=2):
    """
    Given k values, compute the inverse using
    3(k-1) multiplications and one inversion

    The operations are assumed to be in F_{p^extdeg}
    """
    # for values a0, a1, ..., an compute
    # a0, a0a1, a0a1a2, ... a0...an using
    # (k-1) multiplications
    multiples = [input[0]]
    for ai in input[1:]:
        multiples.append(ai * multiples[-1])

    # Compute 1 / (a0a1a2...an)
    last_multiple = multiples[-1]
    inverses_multiples = [1 / last_multiple]

    # Compute (a0a1a2...an)^-1, (a0a1a2...a(n-1)) ... a0^-1
    # using k-1 multiplications
    for ai in reversed(input[1:]):
        inverses_multiples.append(ai * inverses_multiples[-1])

    # Reverse for easy ordering below
    # inverses_multiples = a0^-1, (a0a1)^-1 ...(a0a1...an)^-1
    inverses_multiples = inverses_multiples[::-1]
    inverses = [inverses_multiples[0]]

    # Compute the inverse of each element from multiples and
    # their inverses using k-1 multiplications
    k = len(input)
    for i in range(1, k):
        inverses.append(inverses_multiples[i] * multiples[i - 1])

    # Update costs
    if extdeg == 1:
        cost_count["m"] += 3 * (k - 1)
        cost_count["i"] += 1
    elif extdeg == 2:
        cost_count["M"] += 3 * (k - 1)
        cost_count["I"] += 1

    return inverses


# source: https://github.com/LinKaizhan/Pairingoptimizations/blob/main/verify_sup.sage
def elligator(A, E, Fp, Fp2, cost_count):
    # sample P in E(Fp) and Q in the other eigenspace of Frobenius

    NQR = Fp(2)
    while NQR.is_square():
        NQR += 1

    u = Fp2(NQR)
    sqrt_u = u.sqrt()

    flag_2 = 0
    while flag_2 == 0:
        r = Fp.random_element()
        tmp = Fp(u) * r**2
        cost_count["m"] += 1
        cost_count["s"] += 1
        cond1 = 1 + tmp
        cond2l = A**2 * tmp
        cost_count["m"] += 1
        cost_count["s"] += 1
        cond2r = cond1**2
        cost_count["s"] += 1
        if cond1 != 0 and cond2l != cond2r:
            v = -A / cond1
            cost_count["i"] += 1
            v = Fp(v)
            e = Fp(v * (v**2 + A * v + Fp(1))).is_square()
            # In fact we do not need to compute e as y1 should be computed
            x1 = Fp2(v)
            x2 = -v - A
            y1 = (x1 * (x1**2 + A * x1 + 1)).sqrt()
            cost_count["m"] += 2
            cost_count["s"] += 1
            cost_count["a"] += 3
            cost_count["sqrt"] += 1
            y2 = y1 * sqrt_u * r
            # We can select a small r to save one multiplication
            cost_count["m"] += 1
            if e:
                return E([x1, y1]), E([x2, y2])
            return E([x2, y2]), E([x1, y1])

def seeded_elligator(A, seed, E, Fp, Fp2, cost_count):
    # sample P in E(Fp) and Q in the other eigenspace of Frobenius

    # normally, these values are precomputed given the parameter set
    # however, we use several different primes, hence compute these on the spot
    # we do not count these operations, as they would not affect performance
    # in any real-world implementation
    NQR = Fp(2)
    while NQR.is_square():
        NQR += 1

    u = Fp2(NQR)
    sqrt_u = u.sqrt()

    # start of elligator
    flag_2 = 0
    while flag_2 == 0:
        r = Fp(seed)
        tmp = Fp(u) * r**2
        cost_count["m"] += 1
        cost_count["s"] += 1
        cond1 = 1 + tmp
        cond2l = A**2 * tmp
        cost_count["m"] += 1
        cost_count["s"] += 1
        cond2r = cond1**2
        cost_count["s"] += 1
        if cond1 != 0 and cond2l != cond2r:
            v = -A / cond1
            # in a seeded variant, -1/cond1 can be precomputed
            # so it just takes a single mul
            cost_count["m"] += 1        
            v = Fp(v)
            e = Fp(v * (v**2 + A * v + Fp(1))).is_square()
            # In fact we do not need to compute e as y1 should be computed
            x1 = Fp2(v)
            x2 = -v - A
            y1 = (x1 * (x1**2 + A * x1 + 1)).sqrt()
            cost_count["m"] += 2
            cost_count["s"] += 1
            cost_count["a"] += 3
            cost_count["sqrt"] += 1
            y2 = y1 * sqrt_u * r
            # We can select a small r to save one multiplication
            cost_count["m"] += 1
            if e:
                return E([x1, y1]), E([x2, y2])
            return E([x2, y2]), E([x1, y1])


def lucas_squaring(ft, n, cost_count):
    # ft = f + f.conjugate() for some f in Fp2
    # lucas_sequences(ft, n) should compute f^(2^n) + (f^(2^n)).conjugate()
    v0 = ft    
    for j in range(n):
        v0 = v0 * v0 - 2
        cost_count["s"] += 1
        cost_count["a"] += 1
    return v0

def lucas_sequences(ft, ell, cost_count):
    # ft = f + f.conjugate() for some f in Fp2
    # lucas_sequences(ft, ell) should compute f^ell + (f^ell).conjugate()
    v0 = 2
    v1 = ft
    string = ZZ(ell).str(base=2)
    Len = len(string)
    for j in range(Len):
        if string[j] == "1":
            v0 = v0 * v1 - ft
            cost_count["m"] += 1
            cost_count["a"] += 1
            v1 = v1 * v1 - 2
            cost_count["s"] += 1
            cost_count["a"] += 1
        else:
            v1 = v0 * v1 - ft
            cost_count["m"] += 1
            cost_count["a"] += 1
            v0 = v0 * v0 - 2
            cost_count["s"] += 1
            cost_count["a"] += 1
    return v0


def batch_cof_trace_rec(accs, f, ind, n, L, p, cost_count):
    """
    Given f = tr(e) = e + e.conjugate(),
    compute whether e has multiplicative order larger than 4*sqrt(p) via a product tree

    L is a list of primes, ord(e) divides prod(L)

    accs = [acc_order, acc_enough_order] is a pair of accumulators, initialized as [1,0].
    acc_order: for every prime l in L, the base case (n=1) tests whether e^((p+1)/l)!=1 (i.e. whether f = trace(e) != 2).
    Iff this happens, l divides ord(e), so l is multiplied into acc_order
    acc_enough_order: if acc_order reaches the threshold 4*sqrt(p), set accs[1]:=1 and cut the computation beforehand

    ind, n manage the recursion. first call: n = len(L), ind = [1,...,n]
    """
    if accs[1] == 1:
        return
    if n == 1:
        if f != 2:
            accs[0] *= L[ind[0]]
            if accs[0] > 4 * p.sqrt():  # threshold reached
                accs[1] = 1
        return

    m = floor(n / 2)
    u = 1
    v = 1

    for i in range(m):
        u = L[ind[i]] * u
    left = lucas_sequences(f, u, cost_count)
    batch_cof_trace_rec(
        accs, left, [ind[i] for i in range(m, n)], n - m, L, p, cost_count
    )
    if accs[1] == 1:
        return

    for j in range(m, n):
        v = L[ind[j]] * v
    right = lucas_sequences(f, v, cost_count)
    batch_cof_trace_rec(accs, right, [ind[i] for i in range(m)], m, L, p, cost_count)

def batch_cof_trace_total(accs, f, ind, n, L, p, cost_count):
    """
    Given f = tr(e) = e + e.conjugate(),
    compute whether e has multiplicative order prod(L) via a product tree

    L is a list of primes, ord(e) divides prod(L)

    accs = [acc_order, acc_enough_order] is a pair of accumulators, initialized as [1,0].
    acc_order: for every prime l in L, the base case (n=1) tests whether e^((p+1)/l)!=1 (i.e. whether f = trace(e) != 2).
    Iff this happens, l divides ord(e), so l is multiplied into acc_order
    acc_enough_order: if acc_order reaches the threshold 4*sqrt(p), set accs[1]:=1 and cut the computation beforehand

    ind, n manage the recursion. first call: n = len(L), ind = [1,...,n]
    """
    if accs[1] == 1:
        return
    if n == 1:
        if f != 2:
            accs[0] *= L[ind[0]]
            if accs[0] == prod(L):  # threshold reached
                accs[1] = 1
        return

    m = floor(n / 2)
    u = 1
    v = 1

    for i in range(m):
        u = L[ind[i]] * u
    left = lucas_sequences(f, u, cost_count)
    batch_cof_trace_total(
        accs, left, [ind[i] for i in range(m, n)], n - m, L, p, cost_count
    )
    if accs[1] == 1:
        return

    for j in range(m, n):
        v = L[ind[j]] * v
    right = lucas_sequences(f, v, cost_count)
    batch_cof_trace_total(accs, right, [ind[i] for i in range(m)], m, L, p, cost_count)

# Reference function that uses built-in pairings in Sage
# We're going to implement a version that uses the cubical arithmetic
def dlog_2pt_tate(P, Q, R, S, e, f, ell=2):
    """
    Input: torsion bases <P,Q> = E[ell^e], <R,S> = E[ell^f], where
    - E is an elliptic curve over a finite field Fq,
    - e, f are integers,
    - ell^e is the maximum power of ell dividing #E(Fq)
    """
    tPQ = P.tate_pairing(Q, ell**e, 2)

    D = ell**f
    tRP = R.tate_pairing(P, D, 2)
    tRQ = R.tate_pairing(Q, D, 2)
    tSP = S.tate_pairing(P, D, 2)
    tSQ = S.tate_pairing(Q, D, 2)

    tPQ = tPQ ** (ell ** (e - f))

    # r1 = tRQ.log(tPQ, order=D)
    # r2 = -tRP.log(tPQ, order=D)
    # s1 = tSQ.log(tPQ, order=D)
    # s2 = -tSP.log(tPQ, order=D)
    r1 = discrete_log(tRQ, tPQ, ell**f)
    r2 = -discrete_log(tRP, tPQ, ell**f)
    s1 = discrete_log(tSQ, tPQ, ell**f)
    s2 = -discrete_log(tSP, tPQ, ell**f)

    assert ell**(e - f) * (r1 * P + r2 * Q) == R
    assert ell**(e - f) * (s1 * P + s2 * Q) == S
    return r1, r2, s1, s2


# Translate a string to its naf form
def int_to_naf(N):
    N_str = bin(N)[2:]
    l = len(N_str)
    string_naf = []
    for j in range(l + 1):
        string_naf.append(0)
    j = 0
    while j < l:
        if N_str[l - j - 1] == "1":
            k = j
            while k != l - 1 and N_str[l - k - 2] == "1":
                k = k + 1
            if k != j:
                string_naf[l - j] = -1
                string_naf[l - k - 1] = 1
                j = k
            else:
                string_naf[l - j] = 1
        j = j + 1
    j = 0
    while j < l:
        if string_naf[l - j] == 1:
            k = j
            while k != l and string_naf[l - k - 1] == 1:
                k = k + 1
            if k != j:
                string_naf[l - j] = -1
                for jj in range(j, k):
                    string_naf[l - jj - 1] = 0
                string_naf[l - k - 1] = 1
                j = k
            if string_naf[l - j - 1] == -1:
                string_naf[l - j - 1] = 0
                string_naf[l - j] = -1
                j = j - 1
        j = j + 1
    if string_naf[0] == 0:
        for j in range(l):
            string_naf[j] = string_naf[j + 1]
        string_naf.pop()
    return string_naf
