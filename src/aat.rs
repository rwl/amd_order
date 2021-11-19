use crate::amd::*;
use crate::internal::EMPTY;
use crate::internal::*;
use crate::valid::valid;

pub fn aat(
    n: i32,
    a_p: &[i32],
    a_i: &[i32],
    len: &mut [i32],
    t_p: &mut [i32],
    info: &mut Info,
) -> i32 {
    #[cfg(feature = "debug1")]
    for k in 0..n {
        t_p[k as usize] = EMPTY
    }
    debug_assert!(valid(n, n, a_p, a_i) == Status::OK);

    // Clear the info array, if it exists.
    info.n = 0;
    info.nz = 0;
    info.symmetry = false;
    info.nz_diag = 0;
    info.nz_a_plus_at = 0;
    info.n_dense = 0;
    // for i := 0; i < INFO; i++ {
    //     info[i] = empty
    // }
    info.status = Status::OK;

    for k in 0..n {
        len[k as usize] = 0
    }

    let mut nzdiag: i32 = 0;
    let mut nzboth: i32 = 0;
    let nz: i32 = a_p[n as usize];

    for k in 0..n {
        // let p: i32;
        // let mut pj: i32;
        let p1 = a_p[k as usize];
        let p2 = a_p[k as usize + 1];
        debug2_print!("\nAAT Column: {} p1: {} p2: {}\n", k, p1, p2);

        // Construct A+A'.
        let mut p = p1;
        while p < p2 {
            // Scan the upper triangular part of A.
            let j = a_i[p as usize];
            if j < k {
                // Entry A(j,k) is in the strictly upper triangular part,
                // add both A(j,k) and A(k,j) to the matrix A+A'.
                len[j as usize] += 1;
                len[k as usize] += 1;
                debug3_print!("    upper ({},{}) ({},{})\n", j, k, k, j);
                p += 1;
            } else if j == k {
                // Skip the diagonal.
                p += 1;
                nzdiag += 1;
                break;
            } else {
                // j > k
                // First entry below the diagonal.
                break;
            }

            // Scan lower triangular part of A, in column j until reaching
            // row k. Start where last scan left off.
            debug_assert!(t_p[j as usize] != EMPTY);
            debug_assert!(
                a_p[j as usize] <= t_p[j as usize] && t_p[j as usize] <= a_p[j as usize + 1]
            );

            let pj2 = a_p[j as usize + 1];
            let mut pj = t_p[j as usize];
            while pj < pj2 {
                let i = a_i[pj as usize];
                if i < k {
                    // A(i,j) is only in the lower part, not in upper.
                    // add both A(i,j) and A(j,i) to the matrix A+A'.
                    len[i as usize] += 1;
                    len[j as usize] += 1;
                    debug3_print!("    lower ({},{}) ({},{})\n", i, j, j, i);
                    pj += 1;
                } else if i == k {
                    // Entry A(k,j) in lower part and A(j,k) in upper.
                    pj += 1;
                    nzboth += 1;
                    break;
                } else {
                    // i > k
                    // Consider this entry later, when k advances to i.
                    break;
                }
            }
            t_p[j as usize] = pj;
        }
        // Tp[k] points to the entry just below the diagonal in column k.
        t_p[k as usize] = p;
    }

    // Clean up, for remaining mismatched entries.
    for j in 0..n {
        for pj in t_p[j as usize]..a_p[j as usize + 1] {
            let i = a_i[pj as usize];
            // A(i,j) is only in the lower part, not in upper.
            // add both A(i,j) and A(j,i) to the matrix A+A'.
            len[i as usize] += 1;
            len[j as usize] += 1;
            debug3_print!("    lower cleanup ({},{}) ({},{})\n", i, j, j, i);
        }
    }

    // Compute the symmetry of the nonzero pattern of A.

    // Given a matrix A, the symmetry of A is:
    //  B = tril (spones (A), -1) + triu (spones (A), 1) ;
    //  sym = nnz (B & B') / nnz (B) ;
    //  or 1 if nnz (B) is zero.

    let sym: f64 = if nz == nzdiag {
        1.0
    } else {
        (2.0 * nzboth as f64) / (nz - nzdiag) as f64
    };

    let mut nzaat: i32 = 0;
    for k in 0..n {
        nzaat += len[k as usize];
    }

    debug1_print!("AMD nz in A+A', excluding diagonal (nzaat) = {}\n", nzaat);
    debug1_print!(
        "   nzboth: {} nz: {} nzdiag: {} symmetry: {}\n",
        nzboth,
        nz,
        nzdiag,
        sym
    );

    info.status = Status::OK;
    info.n = n;
    info.nz = nz;
    info.symmetry = sym != 0.0; // Symmetry of pattern of A.
    info.nz_diag = nzdiag; // Nonzeros on diagonal of A.
    info.nz_a_plus_at = nzaat; // Nonzeros in A+A'.

    nzaat
}
