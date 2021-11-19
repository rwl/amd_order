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

pub fn order(n: i32, a_p: &[i32], a_i: &[i32], control: Control) -> (Vec<i32>, Info) {
    let mut p: Vec<i32> = vec![0; n as usize];
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
        info.status = Status::Invalid;
        return (p, info);
    }
    if n == 0 {
        return (p, info);
    }
    let nz: i32 = a_p[n as usize];
    info.nz = nz;
    if nz < 0 {
        info.status = Status::Invalid;
        return (p, info);
    }
    // TODO: check if n or nz will cause size_t overflow.

    // Check the input matrix: OK, INVALID, or OK_BUT_JUMBLED.
    let status = valid(n, n, a_p, a_i);

    if status == Status::Invalid {
        info.status = Status::Invalid;
        return (p, info); // Matrix is invalid.
    }

    let mut len: Vec<i32> = vec![0; n as usize];
    let mut p_inv: Vec<i32> = vec![0; n as usize];

    // Order the input matrix as-is. No need to compute R = A' first.
    // let c_p: &[i32] = a_p;
    // let c_i: &[i32] = a_i;

    let (c_p, c_i) = if status == Status::OkButJumbled {
        // Sort the input matrix and remove duplicate entries.
        debug1_println!("Matrix is jumbled");

        // let r_p: Vec<i32> = vec![0; n as usize + 1];
        // let r_i: Vec<i32> = vec![0; max(nz as usize, 1)];

        // Use Len and Pinv as workspace to create R = A'.
        // preprocess(n, a_p, a_i, &r_p, &r_i, &len, &p_inv);
        // let (r_p, r_i) = preprocess(n, a_p, a_i);
        // let (c_p, c_i) = preprocess(n, a_p, a_i);
        // c_p = &r_p;
        // c_i = &r_i;
        preprocess(n, a_p, a_i)
    } else {
        // Order the input matrix as-is. No need to compute R = A' first.
        // c_p = a_p;
        // c_i = a_i;
        (a_p.to_vec(), a_i.to_vec())
    };

    // Determine the symmetry and count off-diagonal nonzeros in A+A'.

    let nzaat = aat(n, &c_p, &c_i, &mut len, &mut p, &mut info);
    debug1_print!("nzaat: {}\n", nzaat);
    debug_assert!((max(nz - n, 0) <= nzaat) && (nzaat <= 2 * nz));

    // Allocate workspace for matrix, elbow room, and 6 size-n vectors.

    // Space for matrix + elbow room.
    let mut slen = match nzaat.checked_add(nzaat / 5) {
        Some(v) => v,
        None => {
            info.status = Status::OutOfMemory;
            return (p, info);
        }
    };
    // size-n elbow room, 6 size-n work
    for _ in 0..7 {
        slen = match slen.checked_add(n) {
            Some(v) => v,
            None => {
                info.status = Status::OutOfMemory;
                return (p, info);
            }
        };
    }
    // check for overflow
    if !((slen as usize) < usize::MAX / std::mem::size_of::<i32>()) {
        info.status = Status::OutOfMemory;
        return (p, info);
    }
    // S[i] for Int i must be OK
    if !(slen < i32::MAX) {
        info.status = Status::OutOfMemory;
        return (p, info);
    }
    debug1_print!("slen {}\n", slen);

    let mut s: Vec<i32> = vec![0; slen as usize];

    // Order the matrix.

    amd_1(
        n, &c_p, &c_i, &mut p, &mut p_inv, &mut len, slen, &mut s, control, &mut info,
    );

    info.status = status;

    return (p, info);
}
