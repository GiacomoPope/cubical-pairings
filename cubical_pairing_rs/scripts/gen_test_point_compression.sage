proof.all(False)


def field_to_hex(x):
    p = x.base_ring().characteristic()
    n = (p.nbits() + 7) // 8

    re = int(x[0]).to_bytes(n, byteorder="little").hex()
    im = int(x[1]).to_bytes(n, byteorder="little").hex()

    return re + im


def ele_for_rust(s, name):
    r1 = f'pub static {name}_str: &str = "{s}";'
    r2 = f"pub static ({name}, _) = Fq::decode(&hex::decode({name}_str).unwrap());"
    return r1, r2


def field_to_rust(x, name):
    s = field_to_hex(x)
    r1, _ = ele_for_rust(s, name)
    print(r1)


def point_to_rust(P, name):
    PX, PY = P[0], P[1]

    PX_str = field_to_hex(PX)
    PY_str = field_to_hex(PY)

    s1, _ = ele_for_rust(PX_str, f"{name}X")
    r1, _ = ele_for_rust(PY_str, f"{name}Y")

    print(s1)
    print(r1)


def int_to_little(d, name):
    d_list = list(int(d).to_bytes((int(d).bit_length() + 7) // 8, byteorder="little"))
    print(f"pub static {name}: [u8; {len(d_list)}] = {d_list};")
    print(f"pub static {name}_BITLEN: usize = {ZZ(d).nbits()};")
    return None

def random_odd(D):
    a = randint(0, D)
    if a % 2:
        return a
    return a - 1

# Compute indices of powers of g to save for DLP with order n.
def simu_dlp_n_inner(dd, base, lg):
    dd.append(base)
    if lg == 1:
        return
    lg0 = lg >> 1
    lg1 = lg - lg0
    simu_dlp_n_inner(dd, base + lg1, lg0)
    simu_dlp_n_inner(dd, base + lg0, lg1)

def simu_dlp_n(lg):
    dd = []
    simu_dlp_n_inner(dd, 0, lg)
    dd.sort()
    dd2 = []
    for i in range(0, len(dd)):
        v = dd[i]
        j = len(dd2)
        if j == 0 or v != dd2[j - 1]:
            dd2.append(v)
    return dd2

def torsion_basis(E, n, cofactor=1):
    while True:
        P = E.random_point()
        P = P * cofactor
        n_1 = n // 2
        if n_1 * P == E(0):
            continue
        n_2 = n // 3
        if n_2 * P == E(0):
            continue
        break

    assert P.order() == n

    while True:
        Q = E.random_point()
        Q = Q * cofactor
        n_1 = n // 2
        if n_1 * Q == E(0):
            continue
        n_2 = n // 3
        if n_2 * Q == E(0):
            continue

        assert Q.order() == n

        ePQ = P.weil_pairing(Q, n)
        if ePQ**n_1 == 1:
            continue
        if ePQ**n_2 == 1:
            continue

        break

    return P, Q


def gen_data(params):
    ea = params["ea"]
    f = params["f"]

    A = ZZ(2**ea)
    p = f * A - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])

    A0_re = params["A0_re"]
    A0_im = params["A0_im"]

    A0 = F([A0_re, A0_im])
    C0 = 1

    E0 = EllipticCurve(F, [0, A0 / C0, 0, 1, 0])
    E0.set_order((f * A) ** 2, num_checks=0)

    set_random_seed()
    P, Q = torsion_basis(E0, A, cofactor=f)
    # Clear the cofactor
    P, Q = f * P, f * Q
    assert P.order() == 2**ea
    assert Q.order() == 2**ea

    # Set some points of smaller order
    e_sml = ea // 2
    a = random_odd(2**e_sml)
    b = random_odd(2**e_sml)
    c = random_odd(2**e_sml)
    d = random_odd(2**e_sml)

    R = 2**(ea // 2) * (a * P + b * Q)
    S = 2**(ea // 2) * (c * P + d * Q)
    assert R.order() == 2**e_sml
    assert S.order() == 2**e_sml

    # Point data
    if f == 1:
        print(f"// Characteristic: p = 2^{ea} - 1")
    else:
        print(f"// Characteristic: p = {f}*2^{ea} - 1")

    print(f"pub static f: usize = {ea};")
    print(f"pub static e: usize = {e_sml};")
    print()

    print("// Supersingular elliptic curve with Montgomery coeff A0")
    field_to_rust(A0, "A0")
    print()

    print("// <P, Q> are the generators of the curve of order N = 2^f")
    point_to_rust(P, "P")
    point_to_rust(Q, "Q")
    print("// P - Q")
    point_to_rust(P - Q, "PmQ")

    print("// <R, S> are the generators of the curve of order N = 2^e")
    point_to_rust(R, "R")
    point_to_rust(S, "S")
    print("// R - S")
    point_to_rust(R - S, "RmS")
    print()

    # dlog data
    dlog_table = simu_dlp_n(e_sml)
    print("// dlog_table")
    print(f"pub static dlog_table: [usize; {len(dlog_table)}] = {dlog_table};")
    print()

    # Tate data
    dA = (p + 1) // A
    print("// dA = (p + 1) // 2^f")
    int_to_little(dA, "dA")
    print()

    # Compression coefficients we solve for
    int_to_little((a % 2**e_sml), "a")
    int_to_little((-b % 2**e_sml), "b")
    int_to_little((c % 2**e_sml), "c")
    int_to_little((-d % 2**e_sml), "d")



# Primes of the form f * 2^ea - 1
params_1 = {
    "ea": 248,
    "A0_re": 892069292378472517671184093833649435641713298664819259761443361395066011226,
    "A0_im": 2219778098328978246023655261682712071088970775865360962858747592170518415252,
    "f": 5,
}
params_3 = {
    "ea": 376,
    "A0_re": 8087337920290289723003894903718453498620293880852672619403229753674558041137231639497198572980590428497598949898513,
    "A0_im": 4414845428708609956976540215075779401216981780829760904179078553405713283117350468838413838469351402879603005537808,
    "f": 65,
}
params_5 = {
    "ea": 500,
    "A0_re": 76724271569422367173588225530238450404110494734092521370894306217215785109574008324947678991767602396343828379924103029757881684088981207158322047685976,
    "A0_im": 57546393961204530533516766509314029906250597156156989975299314705864190228394792916820208414583735350872378285215412226170867616047967419494900935379496,
    "f": 27,
}
# gen_data(params_1)
gen_data(params_3)
# gen_data(params_5)
