// Characteristic: p = 2^216 * 3^137 - 1
pub static ea: usize = 216;
pub static eb: usize = 137;

// Supersingular elliptic curve with Montgomery coeff A0
pub static A0_str: &str = "9d9ab33208c854c7c5d70b7429e8ed3c8101cf22464be069336d00de3adfd028826eaa63ef47728a9272d2fe76215171f0fa84bb5e22007d6199a0c10dfae6185579efe044ce22560c3bcedbf594753a52d7568d0857ea804273d176e6e60f1d5a14749e03051105ddd8ccb99701";

// <P, Q> are the generators of the curve of order N = 2^216 * 3^137
pub static PX_str: &str = "6a372af6de8ee92fe4f1d52990ec0c84648e05254cf87a961fd1d390f4638d48660c2e4ca17b8d5a93fcb59a9efc0028b39db5e9ad2601c2ce85c2eff902a44e2dd4c35a3f8d004967ff20e4c2f976e91373cc648a27fd6bd6ea4ae235b93af5d3968e125948255415d7362e3b01";
pub static PY_str: &str = "72c6cde6735de52e3900e3e7e00263df51e44c04ba8621b9663e55588b88c38cad71bc082d1c36132e611c0ea1da5689d92dd820f3da0052d85443e810dfe61b24e175ef3574ca6d3c152029344011f6dc1ca8f1911e11e6e2f3e155450369c06919a2cc5fd640e89e13fcbb4101";
pub static QX_str: &str = "3a0faf61f347560686560e692187bfcaff6e7889e442b3ded053c9349fc2e4af4d31cfa7252e023985af73d1064fb80006b3087e490a027feba57e3dc27884b92b0ed7cf9919ab1e02a506d6983edf06433959eb56af2752393aab619f9db4d18ff4f89546befa619d8a297a1300";
pub static QY_str: &str = "adb67c68c3eb82843b7e6dd0768e5f6926565b7d1abdebee615a09247810bf256d5d8670bb70c803b18e676588e16b809aaa02bceb20008b8b95e3d9249cc9504e0b29abe2946a0a44981df837e2adc65ea2ee97a4a2d1e523427ae4c9a51a9d98912dd2785ed14b974fa0795600";
// P - Q
pub static PmQX_str: &str = "f2a676ca9e354357c6647e4928ccbf53c9189361bf4c77daea2d9bb42ccaeb475ad4d7f9542d54de90d3df099d59d26d05d471ca754b001ca3948b8f1b226e9a09e9eca4181d079ecf003488105d29db60d509567c40de455fea0297f70c6b5d3064e8158c1dcee308797c2ae301";
pub static PmQY_str: &str = "db1fccd5f547cff0cd19dc22e0d37783080e5e51615d1cbe28c4af84e7c2333bf1557af6b39a6512ee6fedb4b5c173a4132c4df416190057b858a04a291ec77bd61e3b03fd3ab51a54d2be17fe8ef31b6923949543cc4acaffa2a900847021a8aeae6836068425d28abd6828c300";

// <P2, Q2> are the generators of E0[A]
pub static P2X_str: &str = "3412971c46b2e48d9fc16de5a7111ca678228033058378a68525cacca39f9041e8a9dba3b7457ce0072e4ad30761bc9baf2e6b474b6600f30e7495bbe4ae5c336d22615250311380208dee97b12f22904ded0c322bab96d65b1543713ba1a1b7cd78f0d44293c98b12ad5b4c4b00";
pub static P2Y_str: &str = "bc43139f42d6ec58f97315fbeac5c2536d608d1656d6a1467040b1457c286a63c89a68048c2dc161013c3387098671981b6e72d5442d011deaefb5a9afe6fb28a757cd4c8d47eded03a6ca2f84e0d8a3cfb193fc4b90054e56240e3c82e6ba2916c1a6e0b10e0d6dc5be37e2b801";
pub static Q2X_str: &str = "5a79e05fd69c957a6bacc4459574de4400e7c366644a2c11b4dbda6b479728aff2280061864a3beed2856a260115ed8f08e99f55fe5b0039e8522f8b892484a9bb86e063c6316261e962be80dcf13988a749f007f399a84c8a98031f758da061f764c14a38ac99e98638815f0200";
pub static Q2Y_str: &str = "028219fbb767d8eeff6d237d14ee7ee1e34c2baabdbbb3970d80ef681c0e6bc9b1d1835a59b38c84a55c501cd3f96ee2c39479be229200a6f5d8e10b68097cd5f7ec48be2d74fd5d9240fa94be65518546b3161f6d153ec8eb0cf8a1d9943f79ea36b4aac360e2727e644e02ad01";
// P2 - Q2
pub static P2mQ2X_str: &str = "640b80005988fd480611b498ddee1e34b5d8b9b27ec6c04ddab0cd3c8c183f596f2ef98bd6dcba97e7dff02dc9c4d63816fd05a74820019e9274d996c135e8bb12f46177fc48fd0246dee712e696c0cf52aacb7f6dcd04b8ca3c5080a02e39c8813a050956d4c659c4fb476cb001";
pub static P2mQ2Y_str: &str = "efb40f74bf78c97748352c83d309807210e1ffc6459b8252ab5cf535246c7c2edd8a0e62aca89e79047eb6c0a8664fb81c0d85f11887008599531bfbc1ee4a11ad1b6c1a5dbbf907548a3ecae83c6842f3397512c0944f3fa4cf730b946c94b5a3a51e4ea76367f40d0347c0c100";

// <P3, Q3> are the generators of E0[B]
pub static P3X_str: &str = "9811b1f089f1c6a24361fbc7e0ae603581bc635a97ca89ec314df1318beaea2ea9d651d79f7bcc6136b829ba32bd7f32b2c3443a5412000c0427996562ed0b1792d366605e5e8ff49d06587131d0455a5e37e6dc26c46146e80388ced69fd69a258eed48e90df9a438efec9a0502";
pub static P3Y_str: &str = "0f048bb216001acdf022e22c3b692b4c649d35dbfd570190f58567317ac8f1bf8e2978380f4a45062944632425c1161a8be74326a44c0125e234ca416499e95630a00837d574dffd5270c4c903c6ab972b9c8de7d4552c51dbefef3d815f78436b63158c7957e56552170d0da801";
pub static Q3X_str: &str = "ff868d58f50206db4405ddec9a254c7ab11d8fec33342a42c7ac792a283aff045313379361d407e69fb668382604123b863e1edbecce011001f855b9f581eab459ad4426588dfb4952a2925d68483d26f173b925037a9b200dd2a1acfe2a47e6ed7a9a25b9bbce78e9f429057300";
pub static Q3Y_str: &str = "f416949a43e26c2f6776087d1f105c8d85c2e44bb14f8a1944ecc3e3fc89f0628a8dd517bd92fff61738af7879416cd0c0d56611bb9e01bb9845143b12ab4081f309ba2f5238c0e3d473db87b933b4406d87077fc1b3f79990a3b057474e8be55c1b8054be72a2b93488d590d900";
// P3 - Q3
pub static P3mQ3X_str: &str = "044a12aa7153ed1af74775442eb8a54ae7405f9f465d810494895b421ff69df622ce8d134a65441d6acf818693367857358f5b88dfa201e972a0c223b307552ee19574d5e35041bba32a1631a8a10a9bd19bbd8625030b40eeedc7a7691fcfc4398065a5f53eef947d3b6a58c901";
pub static P3mQ3Y_str: &str = "a18e70a646119cf09f5ac0b20c9ec824cf1bcb9157fcfcf8dccf97f924940002098754396f89266d14b3d5ccbc9ab2c123b537afcf4d0103cb5f729b6badf7e735cf9d135ac5d7317f0347910657ac2c477ed72b5c0ee7d7b23d39f26582c17a5b3d7190631b3d6b1272c4322902";

// [A]P + Q
pub static AP_QX_str: &str = "139165d5f73c15240bc3b8a70d25b1de1337af6ac0e4b9cfafbcb67cf2df6c08ca2641ebd3a251674fffa0ded8c931d882cd0c18313a0161dd34892ca9ab14cf34ba504ee7fae92dcb68477482ef7c271921176e517c3f7cc853f043ca1e1b32fecc6bcfcd127b9cc0c000dbde00";
pub static AP_QY_str: &str = "09d46f1a3c5ea0d6c15dc72a5d44956ebf925085b2f7b37f98ef683250cc5a09305906c7f9a2eaf25aca2b272764f298622454faa5d90022cb5952b14e1b472f172874db3f3235f2a31b606a07f04541de2f4d69715afe4dac71ad68ab7ba75067444361995aee829d9a7bae1101";

// [B]P + Q
pub static BP_QX_str: &str = "ef83e5079133d3c3fc30b6a0bbbb29cdf4dda1d23a61c564ce5c45c71a369a6ab947e15ecf419f22c415729c7073464b77bf7585b95100fc60a327bdd4a7a821f9c6a8fcf7a2ed34a8bbf385d916854779130c34a798ed3e075ecf21db57f411d11f6b929581206f3150f3ab4000";
pub static BP_QY_str: &str = "81b493efdb74d72b28554a01b4d1c64462474e83b76be1d8d7a6afeefd971fef158e22781d3c43e52111c81e8d7fc350c14f3cf0a4d201a4f5d4156be24e31ca70d0c14f1a04b22b1564d08d214c43bccac5bfe14c37ad835a3cf644f119c672719aed206b61eda91119b1180102";

// B = 3^137
pub static B: [u8; 28] = [
    227, 122, 118, 193, 253, 163, 174, 88, 49, 120, 92, 198, 123, 86, 32, 197, 129, 214, 95, 252,
    108, 68, 115, 23, 39, 31, 52, 2,
];
pub static B_BITLEN: usize = 218;
// N = 2^216 * 3^137
pub static N: [u8; 55] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 227, 122, 118,
    193, 253, 163, 174, 88, 49, 120, 92, 198, 123, 86, 32, 197, 129, 214, 95, 252, 108, 68, 115,
    23, 39, 31, 52, 2,
];
pub static N_BITLEN: usize = 434;

// dX = (p + 1) // X
pub static dA: [u8; 28] = [
    227, 122, 118, 193, 253, 163, 174, 88, 49, 120, 92, 198, 123, 86, 32, 197, 129, 214, 95, 252,
    108, 68, 115, 23, 39, 31, 52, 2,
];
pub static dA_BITLEN: usize = 218;
pub static dB: [u8; 28] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
];
pub static dB_BITLEN: usize = 217;
pub static dN: [u8; 1] = [1];
pub static dN_BITLEN: usize = 1;

// Tate pairings for order N, A, B: e_X(Q, P)
pub static tate_pairing_str: &str = "03a5fc5baf13feee0158f691bbdf53a904f23ae2097d06f376f8ff89dbec74beb144c9884b6283dcdd0f37d0b5f07b1a8bfa1089b3b60197282228e36be9b1f8b94e38f11ce7a2c49d1e15776bc1d81eed893486a9fb540fd1e06bbd59720856d1ee535090ec25ebfca7d4a36f01";
pub static tate_pairing_B_str: &str = "665ed36f2b5f5c83ac7b9b0bc3d33d82812a6da55486efca7f76b74b6ddc4718e8d82cf11ad6fa863f12727f8879736b52cc4ba6f932015925f6e6260b3789c6f5146228156cc4ee96b8d88251cb135bbee94f79465932b0692604ccf8690edafb1c65285b9a2585e03f6454d300";
pub static tate_pairing_A_str: &str = "337355edb6fdcc2fcde4bb95cdc383b9cacf126b40887372ed5e84f3303e91b610946413d14c42007824583ca4282c62678010fa971e00db7220bdf3fec4dd8493d5bc6d9121ae92d139a44b20147c2c9a75017609c0165064e9494b5ec36c96c1e8207f6c50ad583b376bb2e200";

// Weil pairings for order N, A, B: e_X(Q, P)
pub static weil_pairing_str: &str = "03a5fc5baf13feee0158f691bbdf53a904f23ae2097d06f376f8ff89dbec74beb144c9884b6283dcdd0f37d0b5f07b1a8bfa1089b3b60168d7ddd71c94164e0746b1c70ee3185d3b62e1ea88943e27e11276aef4ccc5a894dd77c5ba025473004fd62d86cf0f4759766f527bc400";
pub static weil_pairing_B_str: &str = "13125102956fd55f99301ee8743c8139a6d361d893236464bc56f8292aa2653bbeb01fac2daaaa6a7344f5836c3d873bd76796d7a3750038230b4047ead67767ae158126ace065c2d3239c5eb3b3454ca62b318d55b8f975faa6ffa9426187d5efed64529b95ad0abdce17e3c800";
pub static weil_pairing_A_str: &str = "f2f160e0590892f8e64107b1ca763b0e995a778762683ade3ff618f38161e76ff85dfcf1ca0271239ed8baec56c49bcdbb85bdd84bc301484c3355325156eb01e282fa47f24e0aedc0ea6c9a87a761499112c37cfeb3a2a9e96550760f1f039918d3e0974cdbe8e5433361fc3100";
