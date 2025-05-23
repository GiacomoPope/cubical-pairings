"""
This file, given a prime p = 3 mod 4 computes all the following constants needed
for the macros fp_gen.rs and fp2_gen.rs

//      BITLEN (usize)            modulus length in bits
//      MODULUS ([u64; N])        modulus p (little-endian)
//      HALF_MODULUS ([u64; N])   (p + 1)/2 (little-endian)
//      R_VAL ([u64; N])          2^(64*N) mod p
//      DR_VAL ([u64; N])         2*R_VAL mod p = 2^(64*N+1) mod p
//      TR_VAL ([u64; N])         3*R_VAL mod p = 3*2^(64*N) mod p
//      QR_VAL ([u64; N])         4*R_VAL mod p = 2^(64*N+2) mod p
//      R2_VAL ([u64; N])         R_VAL^2 mod p = 2^(2*64*N) mod p
//      P0I (u64)                 -1/p mod 2^64
//      TFIXDIV_VAL               corrective factor for division
//      TDEC_VAL                  2^(64*(2*N-1)) mod p
//      SQRT_EH, SQRT_EL          encoding of (p + 1)/4
//      FOURTH_ROOT_EH, FOURTH_ROOT_EL          encoding of (p + 1)/8
//      P1 (u64)                  floor(p / 2^(BITLEN - 32))
//      P1DIV_M (u64)             1 + floor((2^32 - P1)*2^64 / P1)
//      NQR_RE_VAL                NQR_RE + i is a non-square in GF(p^2)
"""

proof.all(False)

# SIDH primes
# ea, eb = 216, 137
# ea, eb = 305, 192
ea, eb = 372, 239

A = 2**ea
B = 3**eb
p = A * B - 1

# ells = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
# 	103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
# 	223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
# 	349, 353, 359, 367, 373, 587]
# p = 4  * prod(ells) - 1 # CSIDH-512 prime

# SQIsign primes
# p = 5*2**248 - 1
# p = 65*2**376 - 1
# p = 27*2**500 - 1

assert p.is_prime()

def to_little_u64(n):
    res = []
    while n:
        tmp = n % 2**64
        res.append(hex(tmp))
        n >>= 64
    l = fmt_little_u64(res)
    return fmt_list(l)

def fmt_little_u64(res):
    out = []
    for r in res:
        num = r[2:]
        new_num = num.zfill(16)
        new_num = new_num.upper()
        new_num = "0x" + new_num
        out.append(new_num)
    return out

def fmt_list(l):
    l = str(l)
    l = l.replace("[", "")
    l = l.replace("]", "")
    l = l.replace("'", "")
    return l


BITLEN = p.nbits()
N = ceil(BITLEN / 64)
MODULUS = to_little_u64(p)
P_PLUS_ONE = to_little_u64(p+1)
HALF_MODULUS = to_little_u64((p + 1) // 2)
R_VAL = to_little_u64(2**(64*N) % p)
MINUS_R_VAL = to_little_u64(-2**(64*N) % p)
DR_VAL = to_little_u64(2**(64*N+1) % p)
TR_VAL = to_little_u64(3*2**(64*N) % p)
QR_VAL = to_little_u64(2**(64*N + 2) % p)
R2_VAL = to_little_u64(2**(2*64*N) % p)
P0I = inverse_mod(-p, 2**64)

"""
let n1 = floor((2*BITLEN - 34) / 31)
let n2 = 2*BITLEN - 31*n1 - 2
TFIXDIV_VAL = 2^(33*n1 + 64 - n2 + 2*64*N) mod p
"""
n1 = floor((2*BITLEN - 34) // 31)
n2 = 2*BITLEN - 31*n1 - 2
TFIXDIV_VAL = 2**(33*n1 + 64 - n2 + 2*64*N) % p
TFIXDIV_VAL = to_little_u64(TFIXDIV_VAL)

TDEC_VAL = to_little_u64(2**(64*(2*N-1)) % p)

e = (p + 1)//4
SQRT_EL = BITLEN
while True:
    mod = 2**(5*SQRT_EL)
    if e % mod == 0:
        break
    SQRT_EL -= 1

SQRT_EH = e // 2^(5*SQRT_EL)
SQRT_EH = SQRT_EH.digits(2**5)

e = (p + 1)//8
FOURTH_ROOT_EL = BITLEN
while True:
    mod = 2**(5*FOURTH_ROOT_EL)
    if e % mod == 0:
        break
    FOURTH_ROOT_EL -= 1

FOURTH_ROOT_EH = e // 2^(5*FOURTH_ROOT_EL)
FOURTH_ROOT_EH = FOURTH_ROOT_EH.digits(2**5)

P1 = floor(p // 2**(BITLEN - 32))
P1DIV_M = 1 + ((2**32 - P1)*2**64 // P1)


F.<i> = GF(p**2, modulus=[1,0,1])
while True:
    NQR_RE_VAL = randint(0, p)
    if not F([NQR_RE_VAL, 1]).is_square():
        break

assert not F([NQR_RE_VAL, 1]).is_square()
NQR_RE_VAL = to_little_u64(2**(64*N) * NQR_RE_VAL % p)



str = f"""
pub mod Fp{BITLEN} {{
    const N: usize = {N};
    const BITLEN: usize = {BITLEN};
    const MODULUS: [u64; N] = [
        {MODULUS}
    ];
    const HALF_MODULUS: [u64; N] = [
        {HALF_MODULUS}
    ];
    const R_VAL: [u64; N] = [
        {R_VAL}
    ];
    const MINUS_R_VAL: [u64; N] = [
        {MINUS_R_VAL}
    ];
    const DR_VAL: [u64; N] = [
        {DR_VAL}
    ];
    const TR_VAL: [u64; N] = [
        {TR_VAL}
    ];
    const QR_VAL: [u64; N] = [
        {QR_VAL}
    ];
    const R2_VAL: [u64; N] = [
        {R2_VAL}
    ];
    const P0I: u64 = {P0I};
    const TFIXDIV_VAL: [u64; N] = [
        {TFIXDIV_VAL}
    ];
    const TDEC_VAL: [u64; N] = [
        {TDEC_VAL}
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; {len(SQRT_EH)}] = [
        {str(SQRT_EH).replace('[', '').replace(']', '')}
    ];
    const SQRT_EL: usize = {SQRT_EL};
    const P1: u64 = {P1};
    const P1DIV_M: u64 = {P1DIV_M};

    crate::finitefield::fp_gen::define_fp_core!{{}}

    #[cfg(test)]
    mod tests {{
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }}
}}

pub mod Fp{BITLEN}Ext {{
    use super::Fp{BITLEN}::Fp;
    const NQR_RE: Fp = Fp::new([
        {NQR_RE_VAL}
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {{}}
    #[cfg(test)]
    mod tests {{
        crate::finitefield::fp2_gen::define_fp2_tests! {{}}
    }}
}}
"""

print(str)
