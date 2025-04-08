// Characteristic: p = 5*2^248 - 1
pub static ea: usize = 248;

// Supersingular elliptic curve with Montgomery coeff A0
pub static A0_str: &str = "5afe2b2360b0f02d30eb4c7e9180f5fe1b9049a0a56c82bf3819b10fb1e4f801940705ebcd0d4bc103dc2a456cd7b12ccff444e30cafff6c54963cc19159e804";

// <P, Q> are the generators of the curve of order N = 2^ea
pub static PX_str: &str = "618d5b89152bdbf31078931105c5c6b3c449519ad498f273b058c54d611f2904df88997980684f4a4f8db581c6cce1aca5e17a2d1e114e8346e0a5881b90bf01";
pub static PY_str: &str = "fb4498e35cd0fd13e003358c034b3c3080491d247ac88ab3d9f0deaa4f7d0a02225936be0069440af3767be81fa02a9c0fe6b979f777b3cddb1595c9b8828203";
pub static QX_str: &str = "094276e052db1ce67a71192db1fb0e1c4b41d9bb224fc4e77e7c18e04a584b001f7b31453b2ca17ad36b309b91b7a7c380b442a7e28448446b22ec3513889703";
pub static QY_str: &str = "69760c89b852b6df423e4e34438d7c9c9a4527d6fced725443f3c50ba931f202c93fd34ea045a4e12259554748ded777020d3b736cc2b23e0921b36da925cf01";
// P - Q
pub static PmQX_str: &str = "7b621910b3f2610e4de0b4fd07c9c444f417ea5b9b7ea91f93b9faa6bc8d3503bfdbc3174c32855427e7c675e16ef89b9ba80b911d582823ee83359becdd3800";
pub static PmQY_str: &str = "2d31321e87745c1caa2c62d3fcab0e9f22ee054913c388e206f124f3fd0e2e00b6c1b1bf1413018c99a607b2f606f6265f0d05e0d46a903bdb99e6db19cbed00";
// dA = (p + 1) // A
pub static dA: [u8; 1] = [5];
pub static dA_BITLEN: usize = 3;

// Tate pairings for order: e_2^ea(P, Q)
pub static tate_pairing_str: &str = "7558fbdd7877502f500af6cb1a0aaf240545ca8ee0514fb0111a6b097d2fa60131292a781c0418806da296196dfa1e294bbc5a82f0a7a4efe26e20837fe81f03";

// Weil pairings for order: e_2^ea(P, Q)
pub static weil_pairing_str: &str = "83b9cd40d9f7cf693b171ff34ce87cbcbb32fb194615790b2a35c252979f170470e6468c7a1f06963a33f70375fb3a4900c9df8ba90abda41322d82baadf2003";
