use crate::amd::*;
use crate::internal::*;
use crate::valid::*;
use std::cmp::max;

pub fn preprocess(
    n: i32,
    // nz: i32,
    a_p: &[i32],
    a_i: &[i32],
    // r_p: &[i32],
    // r_i: &[i32],
    // w: &[i32],
    // flag: &[i32],
) -> (Vec<i32>, Vec<i32>) {
    debug_assert!(valid(n, n, a_p, a_i) != Status::Invalid);

    let mut w: Vec<i32> = vec![0; n as usize];
    let mut flag: Vec<i32> = vec![0; n as usize];

    // Count the entries in each row of A (excluding duplicates).

    for i in 0..n {
        w[i as usize] = 0; // # of nonzeros in row i (excl duplicates)
        flag[i as usize] = EMPTY; // flag[i] = j if i appears in column j.
    }
    for j in 0..n {
        let p2 = a_p[j as usize + 1];
        for p in a_p[j as usize]..p2 {
            let i = a_i[p as usize];
            if flag[i as usize] != j {
                // Row index i has not yet appeared in column j.
                w[i as usize] += 1; // One more entry in row i.
                flag[i as usize] = j; // Flag row index i as appearing in col j.
            }
        }
    }

    // Compute the row pointers for R.

    let nz: i32 = a_p[n as usize];
    let mut r_p: Vec<i32> = vec![0; n as usize + 1];
    let mut r_i: Vec<i32> = vec![0; max(nz as usize, 1)];

    r_p[0] = 0;
    for i in 0..n {
        r_p[i as usize + 1] = r_p[i as usize] + w[i as usize];
    }
    for i in 0..n {
        w[i as usize] = r_p[i as usize];
        flag[i as usize] = EMPTY
    }

    // Construct the row form matrix R.

    // R = row form of pattern of A.
    for j in 0..n {
        let p2 = a_p[j as usize + 1];
        for p in a_p[j as usize]..p2 {
            let i = a_i[p as usize];
            if flag[i as usize] != j {
                // Row index i has not yet appeared in column j.
                r_i[w[i as usize] as usize] = j; // Put col j in row i.
                w[i as usize] += 1;
                flag[i as usize] = j; // Flag row index i as appearing in col j.
            }
        }
    }

    if DEBUG_LEVEL != 0 {
        debug_assert!(valid(n, n, &r_p, &r_i) == Status::OK);
        for j in 0..n {
            debug_assert!(w[j as usize] == r_p[j as usize + 1])
        }
    }

    return (r_p, r_i);
}
