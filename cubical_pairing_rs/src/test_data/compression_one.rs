// Characteristic: p = 5*2^248 - 1
pub static f: usize = 248;
pub static e: usize = 124;

// Supersingular elliptic curve with Montgomery coeff A0
pub static A0_str: &str = "5afe2b2360b0f02d30eb4c7e9180f5fe1b9049a0a56c82bf3819b10fb1e4f801940705ebcd0d4bc103dc2a456cd7b12ccff444e30cafff6c54963cc19159e804";

// <P, Q> are the generators of the curve of order N = 2^f
pub static PX_str: &str = "54e7eead54a6966230bd4f2cbb1f21e98bfc0c92e29873d7692a32a3082003030f7a39fda4a2c5fdcfd603811e849899ac2a8fa1546ec59e5dcab3ce81efbd02";
pub static PY_str: &str = "2347788a7bf788f352af5f62ba816d9da1d1b4ed81ec370ad88eb9a86d26aa038b0967c455c5de871472c557c4a4d18c2631f47eb6d1d0ab6c0f1ce301907f00";
pub static QX_str: &str = "454eead047b869fe64009c88b4b07d41ac8e36669eb561336258f8f95925df039622947853adb064403b42f487de1b9f8284f144b135cd2849d4c7a5e73a9b01";
pub static QY_str: &str = "a9e4a2d7e6f8a27bb151dfd9c66857a06161f15d4213bb9251dc0492181fae0320a1d2f264ab452b9b2a7d9556b82c759359579343a22c6a1acb5e21844ea702";
// P - Q
pub static PmQX_str: &str = "9e869f279df7d63583bc37328cf538e1b8818f34b9693715cf88620b8ebe8102a32d186ba5227c8276da60c2f952ced4623f6526b299396eab0226b5645e0201";
pub static PmQY_str: &str = "8395af7ff3acd188fe92cff3084d6e64290f7a30d1d23af639b1cc95abde5601b8b89714a5f149030f8b619335ad3d69d0b6a2ddd87cc5ab618514f88a169802";
// <R, S> are the generators of the curve of order N = 2^e
pub static RX_str: &str = "c361fc21e6a92d7f7d10d9f373a5ca7efb158881408645cd7ca605c0f1c17404ff4a7be793aecea93cb7b0b76df3c95dce8e9af09bf0ab29b3e38d2b1d378402";
pub static RY_str: &str = "4afc0d16fd610e6abf56e2689bf99802b12c8db5181e92130812c070821c4f04ff72fc0d6e562f635c36ed70e19a2aa7877bfe57241a1a418aa67a068e8afb02";
pub static SX_str: &str = "3a4d0e7834f5ebbaceb49cd83a0747b226eea16252a0066c4db9383eb11b76022a6b25cb60d7507334e51a305b7bac12b2feec3dd65231905a08138d6f73fc02";
pub static SY_str: &str = "5edf0a7232d16c7fec2ed906aaa2b49667152193eff280989a9a26b725f1280347bf7a55cd915fea5a0c2e4cf2043d393adb73fa3895b710d83630d1c72ff102";
// R - S
pub static RmSX_str: &str = "881733702997e9667a5263611949fd88779d917b09a5bc20756f7d3931b30d0368488845c39c79bac41c2fca46e92cc39b7ea9a6a8a3f9aebca406032c55eb00";
pub static RmSY_str: &str = "68879c33de07f3c4447158a1716c03c9775397446705a3a1657d808cc9c31c02e47c3409bad03a981ac602274cccaf6493c0b21a7c3bb540c66ff1822aceb603";

// dlog_table
pub static dlog_table: [usize; 11] = [0, 62, 93, 108, 109, 116, 117, 120, 121, 122, 123];

// dA = (p + 1) // 2^f
pub static dA: [u8; 1] = [5];
pub static dA_BITLEN: usize = 3;

// Randomized variables to find with dlogs
pub static a: [u8; 16] = [
    231, 236, 12, 226, 96, 63, 8, 148, 235, 149, 78, 41, 234, 78, 223, 10,
];
pub static b: [u8; 16] = [
    29, 36, 137, 134, 164, 32, 173, 196, 36, 87, 201, 117, 70, 176, 216, 10,
];
pub static c: [u8; 16] = [
    13, 1, 11, 179, 59, 21, 131, 227, 185, 72, 84, 22, 159, 90, 10, 14,
];
pub static d: [u8; 16] = [
    237, 22, 11, 123, 195, 142, 41, 178, 254, 241, 222, 40, 213, 125, 223, 5,
];
