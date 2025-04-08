// ========================================================
// Definitions of base fields GF(p) = Z / pZ for the level
// 1, 3 and 5 parameters from the SQIsign NIST submission.
// Constants defined are for the macro and generated from
// gen_fp.sage
// ========================================================

// The three fields we need are
// level one: p = 5*2^248 - 1
// level two: p = 65*3^376 - 1
// level three: p = 27 * 2^500 - 1

// p = 5*2^248 - 1
pub mod Fp248 {
    const N: usize = 4;
    const BITLEN: usize = 251;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x04FFFFFFFFFFFFFF,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0280000000000000,
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000033,
        0x0000000000000000,
        0x0000000000000000,
        0x0100000000000000,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFCC,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x03FFFFFFFFFFFFFF,
    ];
    const DR_VAL: [u64; N] = [
        0x0000000000000066,
        0x0000000000000000,
        0x0000000000000000,
        0x0200000000000000,
    ];
    const TR_VAL: [u64; N] = [
        0x0000000000000099,
        0x0000000000000000,
        0x0000000000000000,
        0x0300000000000000,
    ];
    const QR_VAL: [u64; N] = [
        0x00000000000000CC,
        0x0000000000000000,
        0x0000000000000000,
        0x0400000000000000,
    ];
    const R2_VAL: [u64; N] = [
        0x3333333333333D70,
        0x3333333333333333,
        0x3333333333333333,
        0x0333333333333333,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x49BA5E3BCD35A858,
        0xF7CED916872B020C,
        0x72B020C49BA5E353,
        0x025E353F7CED9168,
    ];
    const TDEC_VAL: [u64; N] = [
        0x3333333333333333,
        0x3333333333333333,
        0x3333333333333333,
        0x0100000000000033,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 1] = [10];
    const SQRT_EL: usize = 49;
    const P1: u64 = 2684354559;
    const P1DIV_M: u64 = 11068046455220847252;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }
}

pub mod Fp383 {
    const N: usize = 6;
    const BITLEN: usize = 383;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x40FFFFFFFFFFFFFF,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x2080000000000000,
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000003,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x3D00000000000000,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFFC,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x03FFFFFFFFFFFFFF,
    ];
    const DR_VAL: [u64; N] = [
        0x0000000000000007,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x3900000000000000,
    ];
    const TR_VAL: [u64; N] = [
        0x000000000000000B,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x3500000000000000,
    ];
    const QR_VAL: [u64; N] = [
        0x000000000000000F,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x3100000000000000,
    ];
    const R2_VAL: [u64; N] = [
        0x3F03F03F03F03F13,
        0x03F03F03F03F03F0,
        0xF03F03F03F03F03F,
        0x3F03F03F03F03F03,
        0x03F03F03F03F03F0,
        0x1D3F03F03F03F03F,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0xD2220C7F2ED57A49,
        0x135335144B641933,
        0x0C221D3394164B45,
        0x46498DFB36D4EF28,
        0x995A997AD0D06037,
        0x106EE23EE8BBD4A1,
    ];
    const TDEC_VAL: [u64; N] = [
        0x03F03F03F03F03F0,
        0xF03F03F03F03F03F,
        0x3F03F03F03F03F03,
        0x03F03F03F03F03F0,
        0xF03F03F03F03F03F,
        0x1000000000000003,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 3] = [16, 0, 1];
    const SQRT_EL: usize = 74;
    const P1: u64 = 2181038079;
    const P1DIV_M: u64 = 17879151965019966409;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }
}

pub mod Fp505 {
    const N: usize = 8;
    const BITLEN: usize = 505;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x01AFFFFFFFFFFFFF,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x00D8000000000000,
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000097,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0130000000000000,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFF68,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x007FFFFFFFFFFFFF,
    ];
    const DR_VAL: [u64; N] = [
        0x000000000000012F,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x00B0000000000000,
    ];
    const TR_VAL: [u64; N] = [
        0x00000000000001C7,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0030000000000000,
    ];
    const QR_VAL: [u64; N] = [
        0x000000000000025E,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0160000000000000,
    ];
    const R2_VAL: [u64; N] = [
        0xED097B425ED0F19A,
        0x097B425ED097B425,
        0x7B425ED097B425ED,
        0x425ED097B425ED09,
        0x5ED097B425ED097B,
        0xD097B425ED097B42,
        0x97B425ED097B425E,
        0x0045ED097B425ED0,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x7593784F6B9C46CD,
        0x0F4EE8861E907F90,
        0xA3E3F57062D8C2AD,
        0xAFFC765B1EFE5FCD,
        0x9B307952438DCE57,
        0x9EBD729237116A4C,
        0x277FEA5B96DBB9D1,
        0x005BFB5EABF79FC1,
    ];
    const TDEC_VAL: [u64; N] = [
        0x097B425ED097B425,
        0x7B425ED097B425ED,
        0x425ED097B425ED09,
        0x5ED097B425ED097B,
        0xD097B425ED097B42,
        0x97B425ED097B425E,
        0xB425ED097B425ED0,
        0x0190000000000097,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 2] = [24, 6];
    const SQRT_EL: usize = 99;
    const P1: u64 = 3623878655;
    const P1DIV_M: u64 = 3416063723386606284;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }
}

// ========================================================
// Definitions of extension fields above the base fields
// GF(p^2) with modulus x^2 + 1 = 0 (using p = 3 mod 4)
// ========================================================

// Level one finite field GF(p^2)
pub mod Fp248Ext {
    use super::Fp248::Fp;
    const NQR_RE: Fp = Fp::new([
        0xAF32723ECE082DD0,
        0xACEA4D9296C5A390,
        0x59B3715FE4B31D02,
        0x03A726359E3001D0,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp383Ext {
    use super::Fp383::Fp;
    const NQR_RE: Fp = Fp::new([
        0xCE4E6F30066908CF,
        0x1570E4E1D9BF237A,
        0x692405E29806D8A7,
        0x743FB8658A28A44E,
        0xAF00638959A256D7,
        0x0FB71531825DBA2D,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
pub mod Fp505Ext {
    use super::Fp505::Fp;
    const NQR_RE: Fp = Fp::new([
        0x3C476591EDCD38EC,
        0x5F55831FD7D64F0C,
        0xBF926122A7726ED3,
        0x8C049BE702B220E5,
        0xA643D3441EB2B4F0,
        0xB3346CFC6D1880C3,
        0xE052C27120DE73C0,
        0x016EC47392F3B353,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
