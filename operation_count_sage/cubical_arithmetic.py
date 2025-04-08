from sage.all import ZZ, ceil, log, discrete_log
from utilities import batched_inversion



# ============================== #
#       Cubical arithmetic       #
# ============================== #


def translate(P, T, cost_count, extdeg=[2, 2]):
    """
    Translate a point P by a point T in E[2]

    extdeg[0]:  =1 --> P is Fp,  =2 --> P is Fp2
    extdeg[1]:  =1 --> T is Fp,  =2 --> T is Fp2
    """

    # Update the costs
    if extdeg == [1, 1]:
        cost_count["m"] += 2
        cost_count["a"] += 6
    elif extdeg == [1, 2] or extdeg == [2, 1]:
        cost_count["Mm"] += 2
        cost_count["A"] += 4
        cost_count["a"] += 2
    elif extdeg == [2, 2]:
        cost_count["M"] += 2
        cost_count["A"] += 6

    if T[1] == 0:
        return (P[0], P[1])

    if T[0] == 0:
        return (P[1], P[0])

    return (P[0] * T[0] - P[1] * T[1], P[0] * T[1] - P[1] * T[0])


def diff_add(P, Q, ixPmQ, cost_count, extdegs=[2, 2, 2], correct=False):
    """
    Compute x(P + Q) given P, Q and 1/x(P-Q) with x-only arithmetic

    ixPmQ = 1/x(P-Q)
    extdeg[0]: =1 --> P is Fp,  =2 --> P is Fp2
    extdeg[1]: =1 --> Q is Fp,  =2 --> Q is Fp2
    extdeg[2]: =1 --> P-Q is Fp,  =2 --> P-Q is Fp2

    correct: the complete formula for the cubical diffAdd requires to divide XPQ and ZPQ by 4.
    For the following pairing computations: Weil; Tate with (p+1)|final_exponentiation_exponent,
    this factor is killed, so we can omit it and save a few operations.
    """

    XP, ZP = P
    XQ, ZQ = Q
    a = XP + ZP
    b = XP - ZP
    c = XQ + ZQ
    d = XQ - ZQ
    da = d * a  # 1m
    cb = c * b  # 1m
    dapcb = da + cb
    damcb = da - cb
    XPQ = dapcb**2 * ixPmQ  # 1m 1s
    ZPQ = damcb**2  # 1s

    if correct:
        XPQ = XPQ / 4
        ZPQ = ZPQ / 4

    # Update the costs
    if extdegs == [1, 1, 1]:
        cost_count["m"] += 3
        cost_count["s"] += 2
        cost_count["a"] += 6
    elif (extdegs == [1, 2, 1]) or (extdegs == [2, 1, 1]):
        cost_count["Mm"] += 3
        cost_count["S"] += 2
        cost_count["a"] += 2
        cost_count["A"] += 4
    elif (extdegs == [1, 2, 2]) or (extdegs == [2, 1, 2]):
        cost_count["M"] += 1
        cost_count["Mm"] += 2
        cost_count["S"] += 2
        cost_count["a"] += 2
        cost_count["A"] += 4
    elif extdegs == [2, 2, 2]:
        cost_count["M"] += 3
        cost_count["S"] += 2
        cost_count["A"] += 6

    return (XPQ, ZPQ)


def diff_add_proj(P, Q, xPmQ, cost_count, extdegs=[2, 2, 2]):
    XP, ZP = P
    XQ, ZQ = Q
    a = XP + ZP
    b = XP - ZP
    c = XQ + ZQ
    d = XQ - ZQ
    da = d * a  # 1M
    cb = c * b  # 1M
    dapcb = da + cb
    damcb = da - cb
    XPQ = dapcb * dapcb  # 1S
    ZPQ = damcb * damcb * xPmQ  # 1M 1S

    # Update the costs
    if extdegs == [1, 1, 1]:
        cost_count["m"] += 3
        cost_count["s"] += 2
        cost_count["a"] += 6
    elif (extdegs == [1, 2, 1]) or (extdegs == [2, 1, 1]):
        cost_count["Mm"] += 3
        cost_count["S"] += 2
        cost_count["a"] += 2
        cost_count["A"] += 4
    elif (extdegs == [1, 2, 2]) or (extdegs == [2, 1, 2]):
        cost_count["M"] += 1
        cost_count["Mm"] += 2
        cost_count["S"] += 2
        cost_count["a"] += 2
        cost_count["A"] += 4
    elif extdegs == [2, 2, 2]:
        cost_count["M"] += 3
        cost_count["S"] += 2
        cost_count["A"] += 6

    return (XPQ, ZPQ)


def double(P, A24, cost_count, extdegs=[2, 2]):
    """
    Compute x([2]P) with x-only arithmetic

    A24 = (A+2)/4, where A is the montgomery coefficient of the curve where P lives

    extdeg [1,1] --> P is Fp and A is Fp
    extdeg [2,1] --> P is Fp2 and A is Fp
    extdeg [2,2] --> P is Fp2 and A is Fp2
    """

    X, Z = P

    XpZ = X + Z
    XmZ = X - Z
    a = XpZ**2  # 1s
    b = XmZ**2  # 1s
    c = a - b
    X2 = a * b  # 1m
    Z2 = c * (b + A24 * c)  # 2m

    # Update the costs
    if extdegs == [1, 1]:
        cost_count["m"] += 3
        cost_count["s"] += 2
        cost_count["a"] += 4
    elif extdegs == [2, 1]:
        cost_count["M"] += 2
        cost_count["Mm"] += 1
        cost_count["S"] += 2
        cost_count["A"] += 4
    elif extdegs == [1, 2]:
        cost_count["m"] += 1
        cost_count["Mm"] += 2
        cost_count["s"] += 2
        cost_count["A"] += 4
    elif extdegs == [2, 2]:
        cost_count["M"] += 3
        cost_count["S"] += 2
        cost_count["A"] += 4

    return (X2, Z2)

def triple(P, A24, A24m1, cost_count, extdegs=[2, 2]):
    """
    Compute x([3]P) with x-only arithmetic. Algorithm 6 in Appendix A of https://sike.org/files/SIDH-spec.pdf

    A24 = (A+2)/4, where A is the montgomery coefficient of the curve where P lives
    A24m1 = A24 - 1 = (A-2)/4

    extdeg [1,1] --> P is Fp and A is Fp
    extdeg [2,1] --> P is Fp2 and A is Fp
    extdeg [2,2] --> P is Fp2 and A is Fp2
    """

    X, Z = P
    
    a = (X - Z)**2
    b = (X + Z)**2
    apb = a + b
    t1 = (X + X)**2 - apb
    # t0 = 2*Z : not used

    # in the SIDH spec:
    # t5 = b * A24
    # t3 = t5 * b
    # t6 = a * A24m1
    # t2 = a * t6

    # t2 = t5 - t6 = b * A24 - a * (A24 - 1) = (b - a) * A24 + a
    # t3 = t2 - t3 = a**2 * A24m1 - b**2 * A24 = (a**2 - b**2) * A24 - a**2 =
    #    = - (b - a) * A24 * (a + b) - a**2 
    s1 = (b - a) * A24
    t2 = s1 + a
    a2 = a**2
    s2 = s1 * apb
    t3 = -s2 - a2

    t1 = t1 * t2

    t2 = (t3 + t1)**2
    X3 = t2 * X

    t1 = (t3 - t1)**2
    Z3 = t1 * Z

    # Update the costs
    if extdegs == [1, 1]:
        cost_count["a"] += 10
        cost_count["m"] += 5
        cost_count["s"] += 6
    elif extdegs == [2, 1]:
        cost_count["A"] += 10
        cost_count["M"] += 4
        cost_count["Mm"] += 1
        cost_count["S"] += 6
    elif extdegs == [1, 2]:
        cost_count["a"] += 5
        cost_count["A"] += 5
        cost_count["Mm"] += 5
        cost_count["s"] += 4
        cost_count["S"] += 2
    elif extdegs == [2, 2]:
        cost_count["A"] += 10
        cost_count["M"] += 5
        cost_count["S"] += 6
    return (X3, Z3)


def ladder_offset(
    P, Q, PQ, ixP, ixQ, ixPQ, A24, n, cost_count, extdegs=[2, 2, 2, 2], correct=False
):
    """
    Montgomery ladder (left to right) with offset point P
    Compute [n]P, [n]P + Q
    given n >= 2, the points P, Q, P-Q and the inverses of their x-coordinates.

    A24 = (A+2)/4, where A is the montgomery coefficient of the curve where P lives

    extdeg[0]: =1 --> P is Fp,  =2 --> P is Fp2
    extdeg[1]: =1 --> Q is Fp,  =2 --> Q is Fp2
    extdeg[2]: =1 --> P-Q is Fp,  =2 --> P-Q is Fp2
    extdeg[3]: =1 --> A24 is Fp,  =2 --> A24 is Fp2

    correct: see diff_add
    """

    return multiladder_offset(
        P, [Q], [PQ], ixP, [ixQ], [ixPQ], A24, n, cost_count, extdegs=extdegs, correct=correct
    )

def multiladder_offset(
    P, Qs, PQs, ixP, ixQs, ixPQs, A24, n, cost_count, extdegs=[2, 2, 2, 2], correct=False
):
    """
    Montgomery ladder (left to right) with offset point P
    Compute [n]P, [n]P + Q
    given n >= 2, the points P, Q, P-Q and the inverses of their x-coordinates.

    A24 = (A+2)/4, where A is the montgomery coefficient of the curve where P lives

    extdeg[0]: =1 --> P is Fp,  =2 --> P is Fp2
    extdeg[1]: =1 --> Q is Fp,  =2 --> Q is Fp2
    extdeg[2]: =1 --> P-Q is Fp,  =2 --> P-Q is Fp2
    extdeg[3]: =1 --> A24 is Fp,  =2 --> A24 is Fp2

    correct: see diff_add
    """

    nP = (1, 0)
    nPP = P
    nPQs = Qs

    if n == 0:
        return (nP, nPQs)

    # Montgomery-ladder
    for bit in bin(n)[2:]:
        R = diff_add(
            nPP,
            nP,
            ixP,
            cost_count,
            extdegs=[extdegs[0], extdegs[0], extdegs[0]],
            correct=correct,
        )
        if bit == "0":
            nPQs = [diff_add(
                nPQ,
                nP,
                ixQ,
                cost_count,
                extdegs=[extdegs[2], extdegs[0], extdegs[1]],
                correct=correct,
            ) for nPQ, ixQ in zip(nPQs, ixQs)]
            nP = double(nP, A24, cost_count, extdegs=[extdegs[0], extdegs[3]])
            nPP = R
        else:
            nPQs = [diff_add(
                nPQ,
                nPP,
                ixPQ,
                cost_count,
                extdegs=[extdegs[2], extdegs[0], extdegs[2]],
                correct=correct,
            ) for nPQ, ixPQ in zip(nPQs, ixPQs)]
            nPP = double(nPP, A24, cost_count, extdegs=[extdegs[0], extdegs[3]])
            nP = R
    return (nP, *nPQs)


def ladder(P, ixP, A24, n, cost_count, extdegs=[1, 1], correct=False):
    """
    Compute x([n]P) with an x-only Montgoomery ladder (left to right)
    given n >= 2, P and 1/x(P)

    A24 = (A+2)/4, where A is the montgomery coefficient of the curve where P lives

    extdeg[0]: =1 --> P is Fp,  =2 --> P is Fp2
    extdeg[1]: =1 --> A24 is Fp,  =2 --> A24 is Fp2

    correct: see diff_add
    """
    nP = (1, 0)
    nPP = P

    if n == 0:
        return nP

    # Montgomery-ladder
    for bit in bin(n)[2:]:
        R = diff_add(
            nPP,
            nP,
            ixP,
            cost_count,
            extdegs=[extdegs[0], extdegs[0], extdegs[0]],
            correct=correct,
        )
        if bit == "0":
            nP = double(nP, A24, cost_count, extdegs=extdegs)
            nPP = R
        else:
            nPP = double(nPP, A24, cost_count, extdegs=extdegs)
            nP = R

    return nP


def ladder_proj(P, A24, n, cost_count, extdegs=[2, 1]):
    """
    usual Montgomery ladder left to right (without offset)
    Montgomery-ladder
    returns [n]P

    extdegs are the extension degrees of the points P and coefficient A + 2 / 4
    """
    nP = (1, 0)
    nPP = P

    for bit in bin(n)[2:]:
        R = diff_add_proj(
            nPP, nP, P[0], cost_count, extdegs=[extdegs[0], extdegs[0], extdegs[0]]
        )
        if bit == "0":
            nP = double(nP, A24, cost_count, extdegs=extdegs)
            nPP = R
        else:
            nPP = double(nPP, A24, cost_count, extdegs=extdegs)
            nP = R
    return nP


def ladder_power2(P, PQ, ixQ, A24, m, cost_count, extdegs=[2, 2, 2, 2], correct=False):
    """
    Compute [2^m]P, [2^m]P - Q (notice the sign change wrt the usual ladder)
    given m >= 0, the points P, PQ = P-Q, and ixQ = 1/x(Q).
    A24 = (A+2)/4, where A is the montgomery coefficient of the curve where P lives

    extdeg[i] for i = 0 (resp. 1, 2, 3) is the degree of the extension field
    where x(P) (resp. x(Q), x(P-Q), A24) is defined

    correct: see diff_add
    """

    return multiladder_power2(
        P, [PQ], [ixQ], A24, m, cost_count, extdegs=extdegs, correct=correct
    )


def multiladder_power2(
    P, PQs, ixQs, A24, m, cost_count, extdegs=[2, 2, 2, 2], correct=False
):
    """
    Compute [2^m]P and [2^m]P - Q (notice the sign change wrt the usual ladder)
    for every Q in a list of points Qs
    
    if [Q1, ..., Qk] = Qs, then
    PQs = [P-Q1, ..., P-Qk] and ixQs = [1/x(Q1), ..., 1/x(Qk)]

    Assumes the extension degrees of the Qs are all the same (see ladder_power2 for more info)

    correct: see diff_add
    """
    nP = P
    nPQs = PQs
    if m == 0:
        return (nP, *PQs)

    # Doublings in the biextension
    for i in range(0, m):
        nPQs = [
            diff_add(
                nPQ,
                nP,
                ixQ,
                cost_count,
                extdegs=[extdegs[2], extdegs[0], extdegs[1]],
                correct=correct,
            )
            for nPQ, ixQ in zip(nPQs, ixQs)
        ]
        nP = double(nP, A24, cost_count, extdegs=[extdegs[0], extdegs[3]])

    return (nP, *nPQs)


def dlog_2pt_tate_biext(
    PQ_basis, RS_basis, diff_points, A24, e, f, cofactor, cost_count
):
    """
    Given bases <P, Q> of E[2^e] and <R, S> of E[2^f] with e >= f, compute
    scalars r1, r2, s1, s2 such that:
    - R = 2^(e-f) * (r1*P + r2*Q),
    - S = 2^(e-f) * (s1*P + s2*Q)

    Assumes the curve and all the points are defined over F_{p^2}.
    """

    # A24 = (A+2)/4, where A is the montgomery coefficient of the curve we're working on
    # PQ_basis, RS_basis and diff_points look like this:
    P, PQ, ixP, ixQ, e = PQ_basis
    R, S, f = RS_basis
    RP, RQ, SP, SQ = diff_points
    # for any A,B:
    # AB is the point A-B in x-only coordinates
    # ixAB = 1/x(A-B)

    # compute ladder for P, Q
    P2e, P2eQ = ladder_power2(P, PQ, ixQ, A24, e - 1, cost_count)

    # we're computing an even pairing:
    # compute multiladder [2**(f-1)]R, [2**(f-1)]R - P, [2**(f-1)]R - Q
    # and the same for S
    R2f, R2fP, R2fQ = multiladder_power2(
        R, [RP, RQ], [ixP, ixQ], A24, f - 1, cost_count
    )
    S2f, S2fP, S2fQ = multiladder_power2(
        S, [SP, SQ], [ixP, ixQ], A24, f - 1, cost_count
    )

    # translate points to get [2**f]R - P and analogous points
    P2eQ = translate(P2eQ, P2e, cost_count, extdeg=2)
    R2fP = translate(R2fP, R2f, cost_count, extdeg=2)
    S2fP = translate(S2fP, S2f, cost_count, extdeg=2)
    R2fQ = translate(R2fQ, R2f, cost_count, extdeg=2)
    S2fQ = translate(S2fQ, S2f, cost_count, extdeg=2)

    # and now points [2**e]P, [2**f]R, [2**f]S
    P2e = translate(P2e, P2e, cost_count, extdeg=2)
    R2f = translate(R2f, R2f, cost_count, extdeg=2)
    S2f = translate(S2f, S2f, cost_count, extdeg=2)

    # the nonreduced pairings are:
    # tPQ = e(P, -Q) = e(P, Q)^(-1) = z(P2eQ) / x(P2e)
    # tRP = e(R, -P) = e(R, P)^(-1) = z(R2fP) / x(R2f)
    # tRQ = e(R, -Q) = e(R, Q)^(-1) = z(R2fQ) / x(R2f)
    # tSP = e(S, -P) = e(S, P)^(-1) = z(S2fP) / x(S2f)
    # tSQ = e(S, -Q) = e(S, Q)^(-1) = z(S2fQ) / x(S2f)

    # they should be raised to the power cof * 2**(e-f) * (p-1)
    # in order to do the p-1 exponentiation, we compute (a/b)**(p-1) = a**(p-1) / b**(p-1) = Frob(a)b / Frob(b)a.
    # note Frob(x) = x^p = x.conjugate() in Fp2
    As = [P2eQ[1], R2fP[1], R2fQ[1], S2fP[1], S2fQ[1]]
    Bs = [P2e[0], R2f[0], R2f[0], S2f[0], S2f[0]]
    nums = [A.conjugate() * B for A, B in zip(As, Bs)]
    dens = [B.conjugate() * A for A, B in zip(As, Bs)]
    cost_count["M"] += 10
    idens = batched_inversion(dens, cost_count, extdeg=2)

    tPQ, tRP, tRQ, tSP, tSQ = [
        (num * iden) ** (cofactor * 2 ** (e - f)) for num, iden in zip(nums, idens)
    ]

    # cof-exponentiation costs (worst-case) log_2(cof) M + log_2(cof) S
    # 2^(e-f)-exponentiation consists of (e-f) S
    cost_count["M"] += 5 * (1 + ceil(log(cofactor, 2)))
    cost_count["S"] += 5 * (e - f + ceil(log(cofactor, 2)))

    # We don't include the cost of the final discrete logs (irrelevant to the performance of pairings)
    # for optimized discrete logs in this context, see eprint 2023/753
    D = ZZ(1 << f)
    r1 = tRQ.log(tPQ, order=D)
    r2 = -tRP.log(tPQ, order=D)
    s1 = tSQ.log(tPQ, order=D)
    s2 = -tSP.log(tPQ, order=D)
    return r1, r2, s1, s2

def dlog_2pt_tate_biext_power3(
    PQ_basis, RS_basis, diff_points, A24, e, f, cofactor, cost_count
):
    """
    Given bases <P, Q> of E[3^e] and <R, S> of E[3^f] with e >= f, compute
    scalars r1, r2, s1, s2 such that:
    - R = 3^(e-f) * (r1*P + r2*Q),
    - S = 3^(e-f) * (s1*P + s2*Q)

    Assumes the curve and all the points are defined over F_{p^2}.
    """

    # A24 = (A+2)/4, where A is the montgomery coefficient of the curve we're working on
    # PQ_basis, RS_basis and diff_points look like this:
    Px, Qx, PQx, ixP, ixQ, ixPQ, e = PQ_basis
    Rx, Sx, RSx, ixR, ixS, ixRS, f = RS_basis
    RPx, RQx, SPx, SQx, ixRP, ixRQ, ixSP, ixSQ = diff_points

    # for any A,B:
    # AB is the point A+B in x-only coordinates
    # ixAB = 1/x(A+B)

    # compute ladder for P, Q
    P3e, P3eQ = ladder_offset(Px, Qx, PQx, ixP, ixQ, ixPQ, A24, 3**e, cost_count, extdegs=[2, 2, 2, 2])

    ### we could instead compute [cof]Q and [cof]Q+P, but this still requires a ladder,
    ### and then we'd need to normalize the result again... don't think we'd spare much.

    # we're computing an even pairing:
    # compute multiladder [3**(f-1)]R, [3**(f-1)]R + P, [3**(f-1)]R + Q
    # and the same for S
    R3f, R3fP, R3fQ = multiladder_offset(
        Rx, [Px, Qx], [RPx, RQx], ixR, [ixP, ixQ], [ixRP, ixRQ], A24, 3**f, cost_count, extdegs=[2, 2, 2, 2]
    )
    S3f, S3fP, S3fQ = multiladder_offset(
        Sx, [Px, Qx], [SPx, SQx], ixS, [ixP, ixQ], [ixSP, ixSQ], A24, 3**f, cost_count, extdegs=[2, 2, 2, 2]
    )

    # the nonreduced pairings are:
    # tPQ = z(P3eQ) / x(P3e)
    # tRP = z(R3fP) / x(R3f)
    # tRQ = z(R3fQ) / x(R3f)
    # tSP = z(S3fP) / x(S3f)
    # tSQ = z(S3fQ) / x(S3f)

    # they should be raised to the power cof * 3**(e3-f) * (p-1)
    # in order to do the p-1 exponentiation, we compute (a/b)**(p-1) = a**(p-1) / b**(p-1) = Frob(a)b / Frob(b)a.
    # note Frob(x) = x^p = x.conjugate() in Fp2
    As = [P3eQ[1], R3fP[1], R3fQ[1], S3fP[1], S3fQ[1]]
    Bs = [P3e[0], R3f[0], R3f[0], S3f[0], S3f[0]]
    nums = [A.conjugate() * B for A, B in zip(As, Bs)]
    dens = [B.conjugate() * A for A, B in zip(As, Bs)]
    cost_count["M"] += 10
    idens = batched_inversion(dens, cost_count, extdeg=2)

    tPQ, tRP, tRQ, tSP, tSQ = [
        (num * iden) ** (cofactor * 3**(e - f)) for num, iden in zip(nums, idens)
    ]

    # cof-exponentiation costs (worst-case) log_2(cof) M + log_2(cof) S
    # 3^(e-f)-exponentiation consists of (e-f) S
    cost_count["M"] += 5 * (1 + ceil(log(cofactor * 3**(e - f), 2)))
    cost_count["S"] += 5 * (e - f + ceil(log(cofactor * 3**(e - f), 2)))

    # We don't include the cost of the final discrete logs (irrelevant to the performance of pairings)
    # for optimized discrete logs in this context, see eprint 2023/753
    D = 3**f
    r1 = discrete_log(tRQ, tPQ, D)
    r2 = -discrete_log(tRP, tPQ, D)
    s1 = discrete_log(tSQ, tPQ, D)
    s2 = -discrete_log(tSP, tPQ, D)
    return r1, r2, s1, s2


def tate_pairing_fp_odd(P, Q, PQ, ixP, ixQ, ixPQ, A24, n, cost_count):
    # return the square of the non-reduced Tate pairing e(P,Q) of order n
    # using the biextension pairing algorithm from our paper
    ##
    # Assumptions:
    # P, Q are in Fp
    # P is an n-torsion point
    # PQ = P-Q might be in Fp2
    # P, Q, PQ are of the form (x, 1)
    # ixP, ixQ, ixPQ are the inverses of x(P), x(Q), x(P-Q)

    nP, nPQ = ladder_offset(
        P, Q, PQ, ixP, ixQ, ixPQ, A24, n, cost_count, extdegs=[1, 1, 2, 1]
    )
    if nP[1] != 0:
        raise ValueError("P is not an n-torsion point")
    return nPQ[1]
    # we'd need to divide by nP.x(), but it belongs to Fp -> killed in the final exponentiation


def tate_pairing_fp_even(P, Q, PQ, ixP, ixQ, ixPQ, A24, N, cost_count):
    # return the non-reduced Tate pairing e(P,Q) of even order N
    # using the biextension even pairing algorithm from our paper
    ##
    # Assumptions:
    # N even
    # P, Q are in Fp
    # P is an n-torsion point
    # PQ = P-Q might be in Fp2
    # P, Q, PQ are of the form (x, 1)
    # ixP, ixQ, ixPQ are the inverses of x(P), x(Q), x(P-Q)

    n = N // 2
    nP, nPQ = ladder_offset(
        P, Q, PQ, ixP, ixQ, ixPQ, A24, n, cost_count, extdegs=[1, 1, 2, 1]
    )
    NPQ = translate(nPQ, nP, cost_count, extdeg=[2, 1])
    if NPQ[0] != Q[0] * NPQ[1]:
        raise ValueError("P is not an n-torsion point")
    cost_count["Mm"] += 1
    return NPQ[0]
    # we'd need to divide by NP.x(), but it belongs to Fp -> killed in the final exponentiation


def tate_pairing_fp(P, Q, PQ, ixP, ixQ, ixPQ, A24, N, cost_count):
    # wrapper for the two cases
    if N % 2 == 0:
        return tate_pairing_fp_even(P, Q, PQ, ixP, ixQ, ixPQ, A24, N, cost_count)
    else:
        return tate_pairing_fp_odd(P, Q, PQ, ixP, ixQ, ixPQ, A24, N, cost_count)
