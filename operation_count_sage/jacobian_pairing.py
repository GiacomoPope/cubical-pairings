from utilities import batched_inversion, torsion_basis, int_to_naf, lucas_sequences
from sage.all import ZZ, ceil, log, GF, discrete_log

##################################
# Compare against state-of-the-art:
# i.e., to our knowledge, eprint 2024/575
##################################

# Modified Jacobian coordinates:
# a point on a Weierstrass curve EE = y^2 = x^3 + a*x + b is represented via
# (X,Y,Z,T) = (xZ^2, yZ^3, Z, aZ^4)
# where x,y are the usual Weierstrass coordinates


# Miller line evaluation for a doubling step
def line_fun_double(xx, yy, zz, X, Y, lam, f, cost_count, extdeg=2):
    # xx,yy,zz,_ = [2m]P in modified Jacobian coordinates
    # X,Y = Q in Weierstrass coordinates
    # lam = slope of the line tangent to E at [m]P
    # f accumulates the line evaluations: return f_{2m,P}(Q)
    # P lies in Fp, Q = (x,iy) with x,y in Fp.
    # f lies in F_{p^2}.
    t1 = zz**2  # 1S
    t2 = t1 * zz  # 1M
    t3 = t1 * X  # 1M
    t3 = t3 - xx
    t4 = t3 * lam  # 1M     # t4 = lam * (X*zz^2 - xx)
    t2 = t2 * Y  # 1M
    t4 = t4 - t2
    t4 = t4 - yy  # t4 = lam * (X*zz^2 - xx) - (Y*zz^3 + yy)
    f = f**2 * t4  # 1M+1S

    if extdeg == 2:
        den = t3 * zz  # 1M             # den = (X*zz^2-xx)*zz
        f = f * den.conjugate()  # 1M   # should be f/den
    # we can replace inversions by *(.).conjugate() because of the final exponentiation

    # Update the costs
    if extdeg == 1:
        cost_count["M"] += 1
        cost_count["S"] += 1
        cost_count["m"] += 4
        cost_count["s"] += 1
        cost_count["a"] += 2
    elif extdeg == 2:
        cost_count["M"] += 7
        cost_count["S"] += 2
        cost_count["A"] += 3

    return f


# Miller line evaluation for a quadrupling step
def line_fun_qua(xx2, yy2, zz2, X, Y, lam1, lam2, f, cost_count, extdeg=2):
    # xx2,yy2,zz2,_ = [2m]P in modified Jacobian coordinates
    # X,Y = Q in Weierstrass coordinates
    # lam1 = slope of the line tangent to E at [m]P
    # lam2 = slope of the line tangent to E at [2m]P
    # f accumulates the line evaluations: return f_{4m,P}(Q)
    # if extdeg == 1, then P is an Fp-point, while Q = (x, iy) with x,y in Fp.
    # f lies in F_{p^2}.

    t1 = zz2**2  # 1S
    t2 = t1 * zz2  # 1M
    t5 = t1 * X  # 1M
    t5 = t5 - xx2
    t6 = t5 * lam1  # 1M    # t6 = lam1 * (X*zz^2-xx)
    t7 = Y * t2  # 1M
    t7 = t7 + yy2  # t7 = Y*zz^3+yy
    t6 = t6 - t7
    f = f**2  # 1S
    f = f * t6  # 1M
    f = f**2  # 1S
    t5 = t5 * lam2  # 1M    # t5 = lam2 * (X*zz^2-xx)
    yy2 = yy2 + yy2
    t6 = yy2 * t7  # 1M     # t6 = 2yy * (Y*zz^3+yy)
    t5 = t5 + t6
    if extdeg == 2:
        f = f * yy2  # 1M
        t5 = t5 * t2  # 1M
    f = -f * t5.conjugate()  # 1M
    ## we can replace inverse by conjugate thanks to the final exponentiation

    # Update the costs
    if extdeg == 1:
        cost_count["M"] += 2
        cost_count["Mm"] += 1  # t6 = yy2 * t7
        cost_count["S"] += 2
        cost_count["m"] += 6
        cost_count["s"] += 1
        cost_count["a"] += 4  # t7 = t7 + yy2 is for free
    elif extdeg == 2:
        cost_count["M"] += 10
        cost_count["S"] += 3
        cost_count["A"] += 5
    return f


# Miller line evaluation for a double-and-addition or double-and-subtraction step in modified Jacobian coordinates
# Cost: 2M+1S+5m
def line_fun_dbladd(
    xx, yy, zz, x1, y1, tmp, tmp1, tmp2, tmp6, tmp0, X, Y, f, cost_count, extdeg=2
):
    """
    For a complete description, see eprint 2024/575

    xx,yy,zz,_ = P in modified Jacobian coordinates
    X,Y = Q in Weierstrass coordinates
    if extdeg == 1, then P is an Fp-point, while Q = (x, iy) with x,y in Fp. Variables tmp, tmpn lie in Fp
    f lies in F_{p^2}.
    """

    # t1 = tmp*(X*zz^2-xx)-(Y*zz^3+yy)*tmp2
    # t12 = tmp6*(X*zz^2-xx)+(Y*zz^3+yy)*tmp2
    t1 = tmp1
    t3 = t1 * X - xx  # 1m 1a
    t4 = tmp * t3  # 1m
    t5 = Y * tmp0  # 1m
    t6 = t5 + yy  # free for extdeg=1
    t7 = t4 - t6  # 1a
    t8 = t7
    t10 = tmp6 * t3  # 1m
    t11 = t6 * tmp2  # 1Mm
    t12 = t10 + t11  # 1a
    if t12 != 0:
        t1 = t8 * t12.conjugate()  # 1M
        f = f**2 * t1  # 1S+1M
    else:
        f = f.parent()(1)

    # Update costs
    if extdeg == 1:
        cost_count["M"] += 2
        cost_count["S"] += 1
        cost_count["Mm"] += 1
        cost_count["m"] += 4
        cost_count["a"] += 3
    elif extdeg == 2:
        cost_count["M"] += 7
        cost_count["S"] += 1
        cost_count["A"] += 4
    return f


# point doubling in modified Jacobian coordinates
# Cost:3m+5s
def DBL(xx, yy, zz, tt, cost_count, extdeg=2):
    t1 = xx**2  # 1S
    t2 = 2 * yy**2  # 1S
    t3 = t2**2  # 1S
    t4 = 2 * t3
    t5 = (xx + t2) ** 2 - t1 - t3  # 1S
    lam = 3 * t1 + tt
    x3 = lam**2 - 2 * t5  # 1S
    y3 = lam * (t5 - x3) - t4  # 1M
    z3 = 2 * yy * zz  # 1M
    t3 = 2 * t4 * tt  # 1M

    # Update the costs
    if extdeg == 1:
        cost_count["m"] += 3
        cost_count["s"] += 5
        cost_count["a"] += 14
    elif extdeg == 2:
        cost_count["M"] += 3
        cost_count["S"] += 5
        cost_count["A"] += 14

    return x3, y3, z3, t3, lam


# point addition in modified Jacobian coordinates
# we assume Z2==1, i.e. the second point is in Weierstrass coordinates
# Cost: 8m+6s
def ADD(xx, yy, zz, x1, y1, AA, cost_count, extdeg=2):
    t1 = zz**2  # 1S
    t2 = x1 * t1 - xx  # 1M
    t3 = t2**2  # 1S
    t4 = 4 * t3
    t5 = t2 * t4  # 1M
    t0 = zz * t1  # 1M
    t6 = y1 * t0 - yy  # 1M
    t7 = 2 * t6
    t8 = t4 * xx  # 1M
    x3 = t7**2 - t5 - 2 * t8  # 1S
    y3 = t7 * (t8 - x3) - 2 * yy * t5  # 2M
    z3 = (zz + t2) ** 2 - t1 - t3  # 1S
    t3 = AA * z3**4  # 2S+1M

    if extdeg == 1:
        cost_count["m"] += 8
        cost_count["s"] += 6
        cost_count["a"] += 11
    elif extdeg == 2:
        cost_count["M"] += 8
        cost_count["S"] += 6
        cost_count["A"] += 11
    return x3, y3, z3, t3, t1, t2, t6, t0

def TPL(
    xx, yy, zz, AA,
    cost_count, extdeg=2
):

    XX = xx**2  # 1s
    YY = yy**2  # 1s
    ZZ = zz**2  # 1s
    Y4 = YY**2  # 1s
    T0 = ZZ**2  # 1s
    T1 = AA * T0    # 1m
    # T1 = tt = AA * zz^4
    M = T1 + 3 * XX # 1a 1c
    MM = M**2 # 1s
    T2 = 6 * ((xx + YY)**2 - XX - Y4) # 1s 1c 3a
    E = T2 - MM # 1a
    EE = E**2 # 1s
    T = 16 * Y4 # 1c
    T3 = (M + E)**2 - EE - MM # 1s 3a
    U = T3 - T # 1a
    T4 = 4 * YY * U # 1m 1c
    X3 = 4 * (xx * EE - T4) # 1m 1a 1c
    T5 = T - U # 1
    T6 = E * EE # 1m
    T7 = U * T5 # 1m
    Y3 = 8 * yy * (T7 - T6) # 1m 1a 1c
    Z3 = (zz + E)**2 - ZZ - EE # 1s 3a

    aux = [E, M, YY, ZZ, U, T3]
    # Update costs
    if extdeg == 1:
        cost_count["m"] += 6
        cost_count["s"] += 10
        cost_count["a"] += 31
    elif extdeg == 2:
        cost_count["M"] += 6
        cost_count["S"] += 10
        cost_count["A"] += 31
    
    return X3, Y3, Z3, aux

def line_fun_tpl(
    xx, yy, zz, X, Y, aux, f, cost_count, extdeg=2
):
    """
    For a complete description, see eprint 2024/575

    xx,yy,zz,_ = P in modified Jacobian coordinates
    X, Y = Q in Weierstrass coordinates

    if extdeg == 1, then P is an Fp-point, while Q = (x, iy) with x,y in Fp. aux=[E,M,YY] lies in Fp
    if extdeg == 2, P and Q are Fp2 points.
    f lies in F_{p^2}.
    """
    E, M, YY, ZZ, U, T3 = aux
    U = U - T3/2 # 2a # division by 2 = conditional_add p, shift right
    T1 = ZZ * X - xx # 1m 1a
    T2 = E * T1 # 1m
    T3 = M * T1 # 1m
    T4 = yy * ZZ # 1m
    T5 = T4 * zz # 1m
    T6 = Y * T5 # 1m
    T7 = T2 * (T3 - 2 * T6 + 2 * YY) # 1m 4a 
    # extdeg=1 -> T3, YY in Fp. T6 in iFp. cost 3a. (2T6 + the rest) is for free
    T8 = E * (T6 + YY) # 1m 1a
    # extdeg=1 -> (T6 + YY) for free
    T1 = U * T1 - 2 * T8 # 1m 2a
    T2 = ZZ * T1 # 1m
    f = f**3 * T7 * T2.conjugate() # 1s 3m

    # Update costs
    if extdeg == 1:
        cost_count["m"] += 6
        cost_count["Mm"] += 4 # T6, T7, T8, T2=ZZ*T1
        cost_count["S"] += 1
        cost_count["M"] += 3
        cost_count["a"] += 7
        cost_count["A"] += 1
    elif extdeg == 2:
        cost_count["M"] += 13
        cost_count["S"] += 1
        cost_count["A"] += 10

    return f

def miller_modified_power2(X, Y, xP, yP, AA, exp, cost_count, extdeg=2):
    # (xP, yP) = P in Weierstrass coordinates
    f = 1
    xx, yy, zz, tt = xP, yP, 1, AA
    for i in range(0, (exp - 1) // 2):
        xx, yy, zz, tt, lam = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        xx2, yy2, zz2, tt, lam2 = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        f = line_fun_qua(xx, yy, zz, X, Y, lam, lam2, f, cost_count, extdeg=extdeg)
        xx, yy, zz = xx2, yy2, zz2
    if ((exp - 1) % 2) == 1:
        xx, yy, zz, tt, lam = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        f = line_fun_double(xx, yy, zz, X, Y, lam, f, cost_count, extdeg=extdeg)

    # now we have just one final step to make (the final tangent line is vertical)
    f = f**2
    if extdeg == 2:
        zzsq = zz**2
        f = f * (X * zzsq - xx) * zzsq.conjugate()

    # Cost of the final step
    if extdeg == 1:
        cost_count["S"] += 1
    elif extdeg == 2:
        cost_count["M"] += 3
        cost_count["S"] += 2
        cost_count["A"] += 1

    return f


def miller_modified_power2_dbl_only(X, Y, xP, yP, AA, exp, cost_count, extdeg=2):
    f = 1
    xx, yy, zz, tt = xP, yP, 1, AA
    for i in range(0, exp - 1):
        xx, yy, zz, tt, lam = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        f = line_fun_double(xx, yy, zz, X, Y, lam, f, cost_count, extdeg=extdeg)

    # now we have just one final step to make (the final tangent line is vertical)
    f = f**2
    if extdeg == 2:
        zzsq = zz**2
        f = f * (X * zzsq - xx) * zzsq.conjugate()

    # Cost of the final step
    if extdeg == 1:
        cost_count["S"] += 1
    elif extdeg == 2:
        cost_count["M"] += 3
        cost_count["S"] += 2
        cost_count["A"] += 1

    return f


def miller_modified_power2_multi(Xs, Ys, xP, yP, AA, exp, cost_count, extdeg=2):
    fs = [1] * len(Xs)
    assert len(Xs) == len(Ys)
    xx, yy, zz, tt = xP, yP, 1, AA
    for i in range(0, (exp - 1) // 2):
        xx, yy, zz, tt, lam = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        xx2, yy2, zz2, tt, lam2 = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        fs = [
            line_fun_qua(xx, yy, zz, X, Y, lam, lam2, f, cost_count, extdeg=extdeg)
            for X, Y, f in zip(Xs, Ys, fs)
        ]
        xx, yy, zz = xx2, yy2, zz2
    if ((exp - 1) % 2) == 1:
        xx, yy, zz, tt, lam = DBL(xx, yy, zz, tt, cost_count, extdeg=extdeg)
        fs = [
            line_fun_double(xx, yy, zz, X, Y, lam, f, cost_count, extdeg=extdeg)
            for X, Y, f in zip(Xs, Ys, fs)
        ]

    # now we have just one final step to make (the final tangent line is vertical)
    fs = [f**2 for f in fs]
    if extdeg == 2:
        zzsq = zz**2
        fs = [f * (X * zzsq - xx) * zzsq.conjugate() for X, f in zip(Xs, fs)]

    # Cost of the final step
    if extdeg == 1:
        cost_count["S"] += len(Xs)
    elif extdeg == 2:
        cost_count["M"] += 3 * len(Xs)
        cost_count["S"] += len(Xs)
        cost_count["A"] += len(Xs)

    return fs

def miller_modified_power3_multi(Xs, Ys, xP, yP, AA, exp, cost_count, extdeg=2):
    fs = [1] * len(Xs)
    assert len(Xs) == len(Ys)
    xx, yy, zz = xP, yP, 1
    for i in range(0, exp - 1):
        xx3, yy3, zz3, aux = TPL(xx, yy, zz, AA, cost_count, extdeg=extdeg)
        fs = [
            line_fun_tpl(xx, yy, zz, X, Y, aux, f, cost_count, extdeg=extdeg)
            for X, Y, f in zip(Xs, Ys, fs)
        ]
        xx, yy, zz = xx3, yy3, zz3

    # now we have just one final step to make (the final tangent line is the line passing through (x0, y0) = (xx/zz^2,yy/zz^3) with multiplicity 3
    # final step: f <- f^3 * line(X,Y)
    # the line function is
    # line(X,Y) = lambda * (X - x0) - (Y - y0)
    # = [(3*xx^2 + a*zz^4)(X * zz^2 - xx) - (Y * 2 * yy * zz^3 - 2 * yy^2)] / [2 * yy * zz^3]
    fs = [f**3 for f in fs]
    YY = yy**2
    ZZ = zz**2
    zz4 = ZZ**2
    M = AA * zz4 + 3 * xx**2
    
    zz3 = ZZ * zz
    den = 2 * yy * zz3
    t1s = [X * ZZ - xx for X in Xs]
    t2s = [Y * den for Y in Ys]
    nums = [M * t1 + 2 * YY - t2 for t1, t2 in zip(t1s, t2s)]
    fs = [f * num for f, num in zip(fs, nums)] 
    if extdeg == 2:
        
        fs = [f * den.conjugate() for f in fs]

    # Cost of the final step
    if extdeg == 1:
        cost_count["m"] += 3 * len(Xs) + 3
        cost_count["M"] += 2 * len(Xs)
        cost_count["S"] += len(Xs)
        cost_count["s"] += 4
        cost_count["a"] += 2 * len(Xs) + 4
    elif extdeg == 2:
        cost_count["M"] += 6 * len(Xs) + 3
        cost_count["S"] += len(Xs) + 4
        cost_count["A"] += 4 * len(Xs) + 3

    return fs


def miller_modified_odd_naf_trace_fp(X, Y, xx, yy, zz, tt, AA, exp, fin, p, cost_count):
    """
    xx,yy,zz,tt = P in modified Jacobian coordinates
    P is an exp-torsion point
    X,Y = Q in Weierstrass coordinates
    P is an Fp-point, while Q = (x, iy) with x,y in Fp.
    f lies in F_{p^2}.

    Compute f_{m,P}(Q) for m = exp
    """
    Fp2 = GF(p**2, "i", modulus=[1, 0, 1])
    f = Fp2(1)
    naf_string = int_to_naf(exp)
    l = len(naf_string)
    flag = 0

    # Note that zz == 1
    x1 = xx
    y1 = yy
    i = 1
    while i != l:
        xx, yy, zz, tt, tmp = DBL(xx, yy, zz, tt, cost_count, extdeg=1)
        if naf_string[i] == 1 and i != l - 1:
            if xx != x1 * zz**2:
                xx2, yy2, zz2, tt, tmp1, tmp2, tmp6, tmp0 = ADD(
                    xx, yy, zz, x1, y1, AA, cost_count, extdeg=1
                )
                f = line_fun_dbladd(
                    xx,
                    yy,
                    zz,
                    x1,
                    y1,
                    tmp,
                    tmp1,
                    tmp2,
                    tmp6,
                    tmp0,
                    X,
                    Y,
                    f,
                    cost_count,
                    extdeg=1,
                )
                xx = xx2
                yy = yy2
                zz = zz2
                i = i + 1
            else:
                flag = 1
        else:
            if naf_string[i] == -1 and i != l - 1:
                if xx != x1 * zz**2:
                    xx2, yy2, zz2, tt, tmp1, tmp2, tmp6, tmp0 = ADD(
                        xx, yy, zz, x1, -y1, AA, cost_count, extdeg=1
                    )
                    f = line_fun_dbladd(
                        xx,
                        yy,
                        zz,
                        x1,
                        y1,
                        tmp,
                        tmp1,
                        tmp2,
                        tmp6,
                        tmp0,
                        X,
                        Y,
                        f,
                        cost_count,
                        extdeg=1,
                    )
                    xx = xx2
                    yy = yy2
                    zz = zz2
                    i = i + 1
            else:
                if naf_string[i] == 0 and naf_string[i + 1] == 0 and i <= l - 2:
                    xx2, yy2, zz2, tt, tmp2 = DBL(xx, yy, zz, tt, cost_count, extdeg=1)
                    f = line_fun_qua(
                        xx, yy, zz, X, Y, tmp, tmp2, f, cost_count, extdeg=1
                    )
                    flag = 0
                    i = i + 2
                    xx = xx2
                    yy = yy2
                    zz = zz2
                else:
                    f = line_fun_double(xx, yy, zz, X, Y, tmp, f, cost_count, extdeg=1)
                    zzs = (
                        zz**2
                    )  # not counted in the cost: can be returned by line_fun_double for free
                    i = i + 1
    label = 0
    if naf_string[l - 1] == 1:
        if yy == -zzs * y1 * zz:
            label = 1
    else:
        if yy == zzs * y1 * zz:
            label = 1
    cost_count["m"] += 2

    return f, label

def dlog_2pt_tate_clz(PQ_basis, RS_basis, AA, e, f, cofactor, cost_count):
    P, Q, e = PQ_basis
    R, S, f = RS_basis
    EE = P.curve()

    # Ensure basis have expected order
    P.set_order(multiple=(1 << e))
    Q.set_order(multiple=(1 << e))
    R.set_order(multiple=(1 << e))
    S.set_order(multiple=(1 << e))
    assert P.order() == (1 << e)
    assert Q.order() == (1 << e)
    assert R.order() == (1 << f)
    assert S.order() == (1 << f)

    # embedding degree 1:
    # the Tate pairing must be computed as e_{2^e}(P,Q) = f_{2^e,P}(Q+T)/f_{2^e,P}(T)
    # for some random point T

    # 1. compute random point T and sum it to evaluation points P,Q
    T, _ = torsion_basis(EE, cofactor, 2**e)
    PT = P + T
    QT = Q + T
    # cost per addition in Jacobian coordinates: 8M + 2S + 13A, counting Z=1
    cost_count["M"] += 16
    cost_count["S"] += 4
    cost_count["A"] += 26
    # compute pairings

    tPQnum, tPQden = miller_modified_power2_multi(
        [QT.x(), T.x()], [QT.y(), T.y()], P.x(), P.y(), AA, e, cost_count, extdeg=2
    )
    tRPnum, tRQnum, tRden = miller_modified_power2_multi(
        [PT.x(), QT.x(), T.x()],
        [PT.y(), QT.y(), T.y()],
        R.x(),
        R.y(),
        AA,
        f,
        cost_count,
        extdeg=2
    )
    tSPnum, tSQnum, tSden = miller_modified_power2_multi(
        [PT.x(), QT.x(), T.x()],
        [PT.y(), QT.y(), T.y()],
        S.x(),
        S.y(),
        AA,
        f,
        cost_count,
        extdeg=2
    )

    # the nonreduced pairings are:
    # tAB = tABnum / tABden

    # they should be raised to the power cof * 2**(e-f) * (p-1)
    # in order to do the p-1 exponentiation, we compute (a/b)**(p-1) = a**(p-1) / b**(p-1) = Frob(a)b / Frob(b)a.
    # note Frob(x) = x^p = x.conjugate() in Fp2
    As = [tPQnum, tRPnum, tRQnum, tSPnum, tSQnum]
    Bs = [tPQden, tRden, tRden, tSden, tSden]
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

def dlog_2pt_tate_clz_power3(PQ_basis, RS_basis, AA, e, f, cofactor, cost_count):
    P, Q, e = PQ_basis
    R, S, f = RS_basis
    EE = P.curve()
    Dfull = 3**e
    D = 3**f

    # Ensure basis have expected order
    P.set_order(multiple=Dfull)
    Q.set_order(multiple=Dfull)
    R.set_order(multiple=Dfull)
    S.set_order(multiple=Dfull)
    assert P.order() == Dfull
    assert Q.order() == Dfull
    assert R.order() == D
    assert S.order() == D

    # embedding degree 1:
    # the Tate pairing must be computed as e_{2^e}(P,Q) = f_{2^e,P}(Q+T)/f_{2^e,P}(T)
    # for some random point T

    # 1. compute random point T and sum it to evaluation points P,Q
    T, _ = torsion_basis(EE, cofactor, Dfull)
    PT = P + T
    QT = Q + T
    # cost per addition in Jacobian coordinates: 8M + 2S + 13A, counting Z=1
    cost_count["M"] += 16
    cost_count["S"] += 4
    cost_count["A"] += 26
    # compute pairings

    tPQnum, tPQden = miller_modified_power3_multi(
        [QT.x(), T.x()], [QT.y(), T.y()], P.x(), P.y(), AA, e, cost_count, extdeg=2
    )
    tRPnum, tRQnum, tRden = miller_modified_power3_multi(
        [PT.x(), QT.x(), T.x()],
        [PT.y(), QT.y(), T.y()],
        R.x(),
        R.y(),
        AA,
        f,
        cost_count,
        extdeg=2
    )
    tSPnum, tSQnum, tSden = miller_modified_power3_multi(
        [PT.x(), QT.x(), T.x()],
        [PT.y(), QT.y(), T.y()],
        S.x(),
        S.y(),
        AA,
        f,
        cost_count,
        extdeg=2
    )

    # the nonreduced pairings are:
    # tAB = tABnum / tABden

    # they should be raised to the power cof * 3**(e-f) * (p-1)
    # in order to do the p-1 exponentiation, we compute (a/b)**(p-1) = a**(p-1) / b**(p-1) = Frob(a)b / Frob(b)a.
    # note Frob(x) = x^p = x.conjugate() in Fp2
    As = [tPQnum, tRPnum, tRQnum, tSPnum, tSQnum]
    Bs = [tPQden, tRden, tRden, tSden, tSden]
    nums = [A.conjugate() * B for A, B in zip(As, Bs)]
    dens = [B.conjugate() * A for A, B in zip(As, Bs)]
    cost_count["M"] += 10
    idens = batched_inversion(dens, cost_count, extdeg=2)

    tPQ, tRP, tRQ, tSP, tSQ = [
        (num * iden) ** (cofactor * 3**(e - f)) for num, iden in zip(nums, idens)
    ]
    # cof-exponentiation costs (worst-case) log_2(cof) M + log_2(cof) S
    # 3^(e-f)-exponentiation consists of (e-f) S
    cost_count["M"] += 5 * (1 + ceil(log(cofactor, 2)))
    cost_count["S"] += 5 * (e - f + ceil(log(cofactor, 2)))

    # sanity checks
    assert tPQ == P.tate_pairing(Q, 3**e, 2)**(3**(e - f))
    assert tRQ == R.tate_pairing(Q, 3**f, 2)
    assert tRP == R.tate_pairing(P, 3**f, 2)
    assert tSP == S.tate_pairing(P, 3**f, 2)
    assert tSQ == S.tate_pairing(Q, 3**f, 2)

    # We don't include the cost of the final discrete logs (irrelevant to the performance of pairings)
    # for optimized discrete logs in this context, see eprint 2023/753
    # r1 = tRQ.log(tPQ, order=D)
    # r2 = -tRP.log(tPQ, order=D)
    # s1 = tSQ.log(tPQ, order=D)
    # s2 = -tSP.log(tPQ, order=D)
    r1 = discrete_log(tRQ, tPQ, 3**f)
    r2 = -discrete_log(tRP, tPQ, 3**f)
    s1 = discrete_log(tSQ, tPQ, 3**f)
    s2 = -discrete_log(tSP, tPQ, 3**f)

    return r1, r2, s1, s2