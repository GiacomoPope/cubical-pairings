#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::too_many_arguments)]
#[allow(unused_macros)]
macro_rules! static_assert {
    ($condition:expr) => {
        let _ = &[()][1 - ($condition) as usize];
    };
}

pub mod csidh_fields;
pub mod cubical_pairings;
mod finitefield;
pub mod sidh_fields;
pub mod sqisign_fields;

pub mod kummer434 {
    pub type Fq = crate::sidh_fields::Fp434Ext::Fp2;
    crate::cubical_pairings::define_cubical_core!(Fq);
}

pub mod kummer610 {
    pub type Fq = crate::sidh_fields::Fp610Ext::Fp2;
    crate::cubical_pairings::define_cubical_core!(Fq);
}

pub mod kummer751 {
    pub type Fq = crate::sidh_fields::Fp751Ext::Fp2;
    crate::cubical_pairings::define_cubical_core!(Fq);
}

// Uses macro arithmetic
pub mod kummer248 {
    pub type Fq = crate::sqisign_fields::Fp248Ext::Fp2;
    crate::cubical_pairings::define_cubical_core!(Fq);
}

pub mod kummer383 {
    pub type Fq = crate::sqisign_fields::Fp383Ext::Fp2;
    crate::cubical_pairings::define_cubical_core!(Fq);
}

pub mod kummer505 {
    pub type Fq = crate::sqisign_fields::Fp505Ext::Fp2;
    crate::cubical_pairings::define_cubical_core!(Fq);
}

pub mod test_data;

#[cfg(test)]
mod test_cubical_pairings_sidh;
#[cfg(test)]
mod test_cubical_pairings_sqisign;

#[cfg(test)]
mod test_point_compression;
