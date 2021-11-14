/// Default is no debug printing (0-4).
pub const DEBUG_LEVEL: i32 = 0;

pub struct Control {
    /// "dense" if degree > dense * sqrt(n)
    dense: f64,
    /// Do aggressive absorption.
    aggressive: bool,
}

pub fn default_control_settings() -> Control {
    return Control {
        dense: 10.0,
        aggressive: true,
    };
}

#[derive(Debug)]
pub struct Info {
    /// Return value of order and l_order.
    pub status: Status,
    /// A is n-by-n.
    pub n: i32,
    /// Number of nonzeros in A.
    pub nz: i32,
    /// Symmetry of pattern (true is sym., false is unsym.)
    pub symmetry: bool,
    /// Number of entries on diagonal.
    pub nz_diag: i32,
    /// nz in A+A'.
    pub nz_a_plus_at: i32,
    /// Number of "dense" rows/columns in A.
    pub n_dense: i32,
}

#[derive(Debug)]
pub enum Status {
    OK,
    // OutOfMemory,
    /// Input arguments are not valid.
    Invalid,
    /// Input matrix is OK for order, but columns were not sorted, and/or
    /// duplicate entries were present. AMD had to do extra work before
    /// ordering the matrix. This is a warning, not an error.
    OkButJumbled,
}
