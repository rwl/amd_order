use crate::amd::*;
use crate::internal::*;
use crate::valid::*;
use std::cmp::max;

pub fn preprocess(n: usize, a_p: &[usize], a_i: &[usize]) -> (Vec<usize>, Vec<usize>) {
    debug_assert!(valid(n, n, a_p, a_i) != Status::Invalid);

    let mut w: Vec<usize> = vec![0; n];
    let mut flag: Vec<isize> = vec![0; n];

    // Count the entries in each row of A (excluding duplicates).

    for i in 0..n {
        w[i] = 0; // # of nonzeros in row i (excl duplicates)
        flag[i] = EMPTY; // flag[i] = j if i appears in column j.
    }
    for j in 0..n {
        let p2 = a_p[j + 1];
        for p in a_p[j]..p2 {
            let i = a_i[p];
            if flag[i] != j as isize {
                // Row index i has not yet appeared in column j.
                w[i] += 1; // One more entry in row i.
                flag[i] = j as isize; // Flag row index i as appearing in col j.
            }
        }
    }

    // Compute the row pointers for R.

    let nz: usize = a_p[n];
    let mut r_p: Vec<usize> = vec![0; n + 1];
    let mut r_i: Vec<usize> = vec![0; max(nz, 1)];

    r_p[0] = 0;
    for i in 0..n {
        r_p[i + 1] = r_p[i] + w[i];
    }
    for i in 0..n {
        w[i] = r_p[i];
        flag[i] = EMPTY
    }

    // Construct the row form matrix R.

    // R = row form of pattern of A.
    for j in 0..n {
        let p2 = a_p[j + 1];
        for p in a_p[j]..p2 {
            let i = a_i[p];
            if flag[i] != j as isize {
                // Row index i has not yet appeared in column j.
                r_i[w[i]] = j; // Put col j in row i.
                w[i] += 1;
                flag[i] = j as isize; // Flag row index i as appearing in col j.
            }
        }
    }

    debug_assert!(valid(n, n, &r_p, &r_i) == Status::OK);
    for j in 0..n {
        debug_assert!(w[j] == r_p[j + 1])
    }

    (r_p, r_i)
}
