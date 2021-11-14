/// Copyright (c) 1996-2015 Timothy A. Davis, Patrick R. Amestoy and Iain S. Duff.
/// Copyright (c) 2011-2021 Richard Lincoln.
/// All Rights Reserved.
mod amd;
mod preprocess;
mod valid;

use amd::*;
use preprocess::*;
use std::cmp::max;
use valid::*;

pub fn order(n: i32, a_p: Vec<i32>, a_i: Vec<i32>, control: Control) -> (Vec<i32>, Info) {
    let p: Vec<i32> = vec![0; n as usize];
    let mut info = Info {
        status: Status::OK,
        n,
        nz: 0,
        symmetry: false,
        nz_diag: 0,
        nz_a_plus_at: 0,
        n_dense: 0,
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

    if matches!(status, Status::Invalid) {
        info.status = Status::Invalid;
        return (p, info); // Matrix is invalid.
    }

    let len: Vec<i32> = vec![0; n as usize];
    let p_inv: Vec<i32> = vec![0; n as usize];

    let c_p: Vec<i32>;
    let c_i: Vec<i32>;

    if matches!(status, Status::OkButJumbled) {
        // Sort the input matrix and remove duplicate entries.
        if DEBUG_LEVEL >= 1 {
            println!("Matrix is jumbled")
        }
        let r_p: Vec<i32> = vec![0; n as usize + 1];
        let r_i: Vec<i32> = vec![0; max(nz as usize, 1)];

        // Use Len and Pinv as workspace to create R = A'.
        preprocess(n, a_p, a_i, r_p, r_i, len, p_inv);
        c_p = r_p;
        c_i = r_i;
    } else {
        // Order the input matrix as-is. No need to compute R = A' first.
        c_p = a_p;
        c_i = a_i;
    }

    return (p, info);
}
