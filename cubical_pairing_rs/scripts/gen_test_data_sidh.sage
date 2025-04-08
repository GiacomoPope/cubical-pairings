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
    eb = params["eb"]

    A = ZZ(2**ea)
    B = ZZ(3**eb)
    p = A * B - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])

    A0_re = params["A0_re"]
    A0_im = params["A0_im"]

    A0 = F([A0_re, A0_im])
    C0 = 1

    E0 = EllipticCurve(F, [0, A0 / C0, 0, 1, 0])
    E0.set_order((A * B) ** 2, num_checks=0)

    set_random_seed(0)
    P, Q = torsion_basis(E0, A * B)

    P_even = B * P
    Q_even = B * Q

    P_odd = A * P
    Q_odd = A * Q

    # Point data
    print(f"// Characteristic: p = 2^{ea} * 3^{eb} - 1")

    print(f"pub static ea: usize = {ea};")
    print(f"pub static eb: usize = {eb};")
    print()

    print("// Supersingular elliptic curve with Montgomery coeff A0")
    field_to_rust(A0, "A0")
    print()

    print(f"// <P, Q> are the generators of the curve of order N = 2^{ea} * 3^{eb}")
    point_to_rust(P, "P")
    point_to_rust(Q, "Q")
    print("// P - Q")
    point_to_rust(P - Q, "PmQ")
    print()

    print("// <P2, Q2> are the generators of E0[A]")
    point_to_rust(P_even, "P2")
    point_to_rust(Q_even, "Q2")
    print("// P2 - Q2")
    point_to_rust(P_even - Q_even, "P2mQ2")
    print()

    print("// <P3, Q3> are the generators of E0[B]")
    point_to_rust(P_odd, "P3")
    point_to_rust(Q_odd, "Q3")
    print("// P3 - Q3")
    point_to_rust(P_odd - Q_odd, "P3mQ3")
    print()

    print("// [A]P + Q")
    point_to_rust(A * P + Q, "AP_Q")
    print()

    print("// [B]P + Q")
    point_to_rust(B * P + Q, "BP_Q")
    print()

    # Order data
    print(f"// B = 3^{eb}")
    int_to_little(B, "B")
    print(f"// N = 2^{ea} * 3^{eb}")
    int_to_little(A * B, "N")
    print()

    assert A * P + Q == (2**ea) * P + Q == P_odd + Q

    # Tate data
    dA = (p + 1) // A
    dB = (p + 1) // B
    dN = (p + 1) // (A * B)
    print("// dX = (p + 1) // X")
    int_to_little(dA, "dA")
    int_to_little(dB, "dB")
    int_to_little(dN, "dN")
    print()

    tate_even = P.tate_pairing(Q, A * B, 2)
    tate_odd = P_odd.tate_pairing(Q_odd, B, 2)**2
    tate_exp2 = P_even.tate_pairing(Q_even, A, 2)

    print("// Tate pairings for order N, A, B: e_X(Q, P)")
    field_to_rust(tate_even, "tate_pairing")
    field_to_rust(tate_odd, "tate_pairing_B")
    field_to_rust(tate_exp2, "tate_pairing_A")
    print()

    weil_even = P.weil_pairing(Q, A * B)
    weil_odd = P_odd.weil_pairing(Q_odd, B)**2
    weil_exp2 = P_even.weil_pairing(Q_even, A)

    print("// Weil pairings for order N, A, B: e_X(Q, P)")
    field_to_rust(weil_even, "weil_pairing")
    field_to_rust(weil_odd, "weil_pairing_B")
    field_to_rust(weil_exp2, "weil_pairing_A")
    print()


params_one = {
    "ea": 216,
    "eb": 137,
    "A0_re": 1489012386689001439522300804648263213802201861511044849633226372238834666277740589462210996176870689046702546639505171240702941853,
    "A0_im": 17663889362144477580916585759944911743158483023725007140786215218276755897899045981657178300860367326649717635169412074973557252477,
}
params_three = {
    "ea": 305,
    "eb": 192,
    "A0_re": 554917942424378647726098346198556632856145542832893246192491776847667700420026451774953542222717515334389069653754312443832045513483840253674583072478945068996120806710460057849880061,
    "A0_im": 235732700173042386372290354964724057533819329507651563907467905662720643660628280001865665559548402704441406672937648641621346278277867417177498644612381803097315826114604799848367203,
}
params_five = {
    "ea": 372,
    "eb": 239,
    "A0_re": 7340221077745800903372940483520025152604782541972823594355133892260005495472877675032833078688559820718084701079247134017423577349510237026133415607364152476607718590403950030369126815596149616631501160852844457680880097381214,
    "A0_im": 3246082702913096380124014081123689021953141035417712366780826029403784999421024442220533602325865806287332165940482725118509557492218860762666633419613871768614097854344767137095633588313931660735073952748080310523301225786799,
}

gen_data(params_one)
# gen_data(params_three)
# gen_data(params_five)
