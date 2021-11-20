/// Copyright (c) 1996-2015 Timothy A. Davis, Patrick R. Amestoy and Iain S. Duff.
/// Copyright (c) 2011-2021 Richard Lincoln.
/// All Rights Reserved.
mod aat;
pub mod amd;
mod amd_1;
mod amd_2;
pub mod control;
mod dump;
pub mod info;
mod internal;
mod post_tree;
mod postorder;
mod preprocess;
mod valid;

use aat::aat;
use amd::*;
use amd_1::amd_1;
use internal::*;
use preprocess::preprocess;
use std::cmp::max;
use valid::valid;

pub fn order(
    n: i32,
    a_p: &[i32],
    a_i: &[i32],
    control: &Control,
) -> Result<(Vec<i32>, Info), Status> {
    let mut info = Info {
        status: Status::OK,
        n,
        nz: 0,
        symmetry: false,
        nz_diag: 0,
        nz_a_plus_at: 0,
        n_dense: 0,
        n_cmp_a: 0,
        lnz: 0,
        n_div: 0,
        n_mult_subs_ldl: 0,
        n_mult_subs_lu: 0,
        d_max: 0,
    };

    if n < 0 {
        return Err(Status::Invalid);
    }
    if n == 0 {
        let p: Vec<i32> = Vec::new();
        return Ok((p, info));
    }
    let nz: i32 = a_p[n as usize];
    info.nz = nz;
    if nz < 0 {
        return Err(Status::Invalid);
    }

    // Check the input matrix.
    let status = valid(n, n, a_p, a_i);
    if status == Status::Invalid {
        return Err(Status::Invalid);
    }

    let (c_p, c_i) = if status == Status::OkButJumbled {
        // Sort the input matrix and remove duplicate entries.
        debug1_println!("Matrix is jumbled");
        preprocess(n, a_p, a_i) // R = A'.
    } else {
        // Order the input matrix as-is. No need to compute R = A' first.
        (a_p.to_vec(), a_i.to_vec())
    };

    // Determine the symmetry and count off-diagonal nonzeros in A+A'.
    let (nzaat, mut len) = aat(n, &c_p, &c_i, &mut info);
    debug1_print!("nzaat: {}\n", nzaat);
    debug_assert!((max(nz - n, 0) <= nzaat) && (nzaat <= 2 * nz));

    let iwlen = nzaat + (nzaat / 5) + n; // Space for matrix + elbow room.
    debug1_print!("iwlen {}\n", iwlen);

    // Order the matrix.
    let (p, _p_inv) = amd_1(n, &c_p, &c_i, &mut len, iwlen, &control, &mut info);

    info.status = status;

    return Ok((p, info));
}
