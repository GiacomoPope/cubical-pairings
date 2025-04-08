pub mod Fp434 {
    const N: usize = 7;
    const BITLEN: usize = 434;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFDC1767AE2FFFFFF,
        0x7BC65C783158AEA3,
        0x6CFC5FD681C52056,
        0x0002341F27177344,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0xFEE0BB3D71800000,
        0x3DE32E3C18AC5751,
        0x367E2FEB40E2902B,
        0x00011A0F938BB9A2,
    ];
    const R_VAL: [u64; N] = [
        0x000000000000742C,
        0x0000000000000000,
        0x0000000000000000,
        0xB90FF404FC000000,
        0xD801A4FB559FACD4,
        0xE93254545F77410C,
        0x0000ECEEA7BD2EDA,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFF8BD3,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x44B18275E6FFFFFF,
        0xA3C4B77CDBB901CF,
        0x83CA0B82224DDF49,
        0x000147307F5A4469,
    ];
    const DR_VAL: [u64; N] = [
        0x000000000000E858,
        0x0000000000000000,
        0x0000000000000000,
        0x721FE809F8000000,
        0xB00349F6AB3F59A9,
        0xD264A8A8BEEE8219,
        0x0001D9DD4F7A5DB5,
    ];
    const TR_VAL: [u64; N] = [
        0x0000000000015C85,
        0x0000000000000000,
        0x0000000000000000,
        0x2D6E659411000000,
        0x0C3E9279CF8657DA,
        0x4E9A9D269CA0A2D0,
        0x000092ACD020194C,
    ];
    const QR_VAL: [u64; N] = [
        0x000000000001D0B1,
        0x0000000000000000,
        0x0000000000000000,
        0xE67E59990D000000,
        0xE4403775252604AE,
        0x37CCF17AFC17E3DC,
        0x00017F9B77DD4827,
    ];
    const R2_VAL: [u64; N] = [
        0x28E55B65DCD69B30,
        0xACEC7367768798C2,
        0xAB27973F8311688D,
        0x175CC6AF8D6C7C0B,
        0xABCD92BF2DDE347E,
        0x69E16A61C7686D9A,
        0x000025A89BCDD12A,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x6EC1344E70FF5EE9,
        0xB1FB01EE3A3DA2AD,
        0x8318449945BCD2B3,
        0x4056D65B06D4D6BE,
        0x344EE62DA62887A9,
        0x2F088C9F003B84A8,
        0x0000C9758E47EDE0,
    ];
    const TDEC_VAL: [u64; N] = [
        0xACEC7367768798C2,
        0xAB27973F8311688D,
        0x032D272B1D6C7C0B,
        0x68468FD91D3B1B99,
        0xD2085CC1A8F7FE9F,
        0x3D12564C40A5CB75,
        0x00005A1E47612BD0,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 45] = [
        16, 17, 11, 15, 22, 11, 16, 27, 31, 17, 26, 21, 24, 10, 12, 16, 7, 14, 25, 24, 27, 19, 21,
        0, 18, 2, 7, 16, 22, 30, 23, 24, 15, 22, 17, 8, 19, 27, 5, 14, 18, 15, 16, 6, 2,
    ];
    const SQRT_EL: usize = 42;
    const P1: u64 = 2366097861;
    const P1DIV_M: u64 = 15037992048227349790;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }
}

pub mod Fp610 {
    const N: usize = 10;
    const BITLEN: usize = 610;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x6E01FFFFFFFFFFFF,
        0xB1784DE8AA5AB02E,
        0x9AE7BF45048FF9AB,
        0xB255B2FA10C4252A,
        0x819010C251E7D88C,
        0x000000027BF6A768,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x3701000000000000,
        0xD8BC26F4552D5817,
        0x4D73DFA28247FCD5,
        0x592AD97D08621295,
        0x40C8086128F3EC46,
        0x000000013DFB53B4,
    ];
    const R_VAL: [u64; N] = [
        0x00000000670CC8E6,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x9A34000000000000,
        0x4D99C2BD28717A3F,
        0x0A4A1839A323D41C,
        0xD2B62215D06AD1E2,
        0x1369026E862CAF3D,
        0x000000010894E964,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFF98F33719,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xD3CDFFFFFFFFFFFF,
        0x63DE8B2B81E935EE,
        0x909DA70B616C258F,
        0xDF9F90E440595348,
        0x6E270E53CBBB294E,
        0x000000017361BE04,
    ];
    const DR_VAL: [u64; N] = [
        0x00000000CE1991CC,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x3468000000000000,
        0x9B33857A50E2F47F,
        0x149430734647A838,
        0xA56C442BA0D5A3C4,
        0x26D204DD0C595E7B,
        0x000000021129D2C8,
    ];
    const TR_VAL: [u64; N] = [
        0x0000000135265AB3,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x609A000000000000,
        0x3754FA4ECEF9BE90,
        0x83F68967E4DB82A9,
        0xC5CCB347607C507B,
        0xB8AAF689409E352C,
        0x000000009DC814C3,
    ];
    const QR_VAL: [u64; N] = [
        0x000000019C332399,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0xFACE000000000000,
        0x84EEBD0BF76B38CF,
        0x8E40A1A187FF56C5,
        0x9882D55D30E7225D,
        0xCC13F8F7C6CAE46A,
        0x00000001A65CFE27,
    ];
    const R2_VAL: [u64; N] = [
        0xE75F5D201A197727,
        0xE0B85963B627392E,
        0x6BC1707818DE493D,
        0xDC7F419940D1A0C5,
        0x7358030979EDE54A,
        0x84F4BEBDEED75A5C,
        0x7ECCA66E13427B47,
        0xC5BB4E65280080B3,
        0x7019950F516DA19A,
        0x000000008E290FF3,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x54AB783CEEFDBBFC,
        0xBF3FF1BFC1DC7E4B,
        0x65AA12D3E2EC2CC0,
        0xC9E3641B740CA026,
        0x7EBEB686F85AA8F6,
        0xD5B0AF649E60A215,
        0xA97EA657ED228146,
        0x189B08459E89E4DF,
        0xF022C2062407CFB6,
        0x000000006CA6AE86,
    ];
    const TDEC_VAL: [u64; N] = [
        0xE0B85963B627392E,
        0x6BC1707818DE493D,
        0xDC7F419940D1A0C5,
        0x23A6030979EDE54A,
        0x3B5B21A27082B8C4,
        0xF5C4EF5C42713031,
        0x3DB71867A73C3941,
        0x4BE9B30B7DEB9228,
        0xC72D47CF2708B75A,
        0x000000023EC878EF,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 62] = [
        8, 0, 14, 19, 11, 0, 11, 13, 9, 21, 8, 15, 19, 16, 23, 24, 14, 21, 25, 31, 3, 9, 16, 2, 29,
        23, 7, 23, 6, 21, 18, 18, 16, 24, 16, 16, 30, 5, 27, 10, 9, 22, 12, 4, 22, 15, 30, 8, 9,
        24, 16, 0, 4, 3, 8, 20, 29, 20, 22, 31, 30, 4,
    ];
    const SQRT_EL: usize = 60;
    const P1: u64 = 2667424218;
    const P1DIV_M: u64 = 11255379038027277338;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }
}

pub mod Fp751 {
    const N: usize = 12;
    const BITLEN: usize = 751;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xEEAFFFFFFFFFFFFF,
        0xE3EC968549F878A8,
        0xDA959B1A13F7CC76,
        0x084E9867D6EBE876,
        0x8562B5045CB25748,
        0x0E12909F97BADC66,
        0x00006FE5D541F71C,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x7758000000000000,
        0x71F64B42A4FC3C54,
        0x6D4ACD8D09FBE63B,
        0x04274C33EB75F43B,
        0x42B15A822E592BA4,
        0x0709484FCBDD6E33,
        0x000037F2EAA0FB8E,
    ];
    const R_VAL: [u64; N] = [
        0x00000000000249AD,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x8310000000000000,
        0x5527B1E4375C6C66,
        0x697797BF3F4F24D0,
        0xC89DB7B2AC5C4E2E,
        0x4CA4B439D2076956,
        0x10F7926C7512C7E9,
        0x00002D5B24BCE5E2,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFDB652,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x6B9FFFFFFFFFFFFF,
        0x8EC4E4A1129C0C42,
        0x711E035AD4A8A7A6,
        0x3FB0E0B52A8F9A48,
        0x38BE00CA8AAAEDF1,
        0xFD1AFE3322A8147D,
        0x0000428AB0851139,
    ];
    const DR_VAL: [u64; N] = [
        0x000000000004935A,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0620000000000000,
        0xAA4F63C86EB8D8CD,
        0xD2EF2F7E7E9E49A0,
        0x913B6F6558B89C5C,
        0x99496873A40ED2AD,
        0x21EF24D8EA258FD2,
        0x00005AB64979CBC4,
    ];
    const TR_VAL: [u64; N] = [
        0x000000000006DD08,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x9A80000000000000,
        0x1B8A7F275C1CCC8A,
        0x61D12C23A9F5A1FA,
        0x518A8EB02E290214,
        0x608B67A91963E4BC,
        0x24D426A5C77D7B55,
        0x0000182B98F4BA8A,
    ];
    const QR_VAL: [u64; N] = [
        0x00000000000926B5,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x1D90000000000000,
        0x70B2310B937938F1,
        0xCB48C3E2E944C6CA,
        0x1A284662DA855042,
        0xAD301BE2EB6B4E13,
        0x35CBB9123C90433E,
        0x00004586BDB1A06C,
    ];
    const R2_VAL: [u64; N] = [
        0x233046449DAD4058,
        0xDB010161A696452A,
        0x5E36941472E3FD8E,
        0xF40BFE2082A2E706,
        0x4932CCA8904F8751,
        0x1F735F1F1EE7FC81,
        0xA24F4D80C1048E18,
        0xB56C383CCDB607C5,
        0x441DD47B735F9C90,
        0x5673ED2C6A6AC82A,
        0x06C905261132294B,
        0x000041AD830F1F35,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0xD1EABC64CC4EFA21,
        0x873B72C73A9B339D,
        0xAD1F614319E93501,
        0xBDB7608700747D23,
        0x325D78AA668F2BF8,
        0xD212376F2C70A082,
        0xAF7D8F367315A5D6,
        0xCC3F9E167761967A,
        0xF00F9956B17A232A,
        0x25D594B94082E085,
        0x268832D4EEF03F84,
        0x000068E3F72EEF89,
    ];
    const TDEC_VAL: [u64; N] = [
        0xDB010161A696452A,
        0x5E36941472E3FD8E,
        0xF40BFE2082A2E706,
        0x4932CCA8904F8751,
        0x2BF35F1F1EE7FC81,
        0xE701CBDCF7E380C6,
        0xB5AD6BE3623D6D87,
        0xE09CC87DE8B9A717,
        0x785070DF7EB808F2,
        0xDA95D958B6D2447A,
        0x88813CAD423F06A8,
        0x00000F6185F6D76F,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 76] = [
        11, 23, 3, 21, 24, 3, 30, 19, 20, 2, 26, 18, 12, 31, 24, 13, 7, 6, 31, 30, 19, 16, 6, 22,
        25, 10, 10, 27, 22, 3, 26, 23, 14, 11, 31, 12, 24, 20, 19, 16, 0, 4, 29, 10, 18, 5, 23, 8,
        16, 26, 10, 12, 5, 20, 25, 24, 13, 29, 30, 18, 31, 4, 4, 5, 1, 7, 16, 3, 23, 15, 16, 10,
        29, 18, 31, 13,
    ];
    const SQRT_EL: usize = 74;
    const P1: u64 = 3754666627;
    const P1DIV_M: u64 = 2654506818854534766;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests!(true);
    }
}

pub mod Fp434Ext {
    use super::Fp434::Fp;
    const NQR_RE: Fp = Fp::new([
        0x084FFA2A435F0EF9,
        0xDFC67AB5FEE99F99,
        0x9128233843682A94,
        0x33D19614142B2766,
        0xE782A1CA42C4E740,
        0xA84983D1E40248E1,
        0x0000BF9243337C09,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp610Ext {
    use super::Fp610::Fp;
    const NQR_RE: Fp = Fp::new([
        0x3766EB743BCFE364,
        0xBFE71485E58953C6,
        0xF848117C05A975ED,
        0xEA4B1EFADC56A518,
        0xC0C06297112A2E23,
        0x27F0443E4A7418D1,
        0x9021BED2483B9118,
        0x02BD492CBF18D3BC,
        0xA2EC1D766B473968,
        0x000000005DC6E474,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp751Ext {
    use super::Fp751::Fp;
    const NQR_RE: Fp = Fp::new([
        0xA38600AAD6A0F10A,
        0x8D018C7D5D062B70,
        0x26E2DB70F4A088FD,
        0x1F49D5B0EB270D8B,
        0x83837C1C593A065F,
        0x0DE4E7F46AE68EEC,
        0x0C75D344F2E5B297,
        0xCBB77DAC1BBE16D2,
        0x6E25DB27BB360AF0,
        0x38E6A237E0A9A7C6,
        0x44863537701E71B7,
        0x000019F570B0592E,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
