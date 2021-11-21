use crate::amd::*;
use crate::amd_2::amd_2;
use crate::internal::*;
use crate::valid::valid;

pub fn amd_1(
    n: usize,
    a_p: &[usize],
    a_i: &[usize],
    mut len: &mut [usize],
    iwlen: usize,
    control: &Control,
    mut info: &mut Info,
) -> (Vec<usize>, Vec<usize>) {
    // Construct the matrix for amd_2.

    debug_assert!(n > 0);

    let mut p_e: Vec<isize> = vec![0; n];
    let mut s_p: Vec<usize> = vec![0; n];
    let mut t_p: Vec<usize> = vec![0; n];
    let mut i_w: Vec<isize> = vec![0; iwlen];

    debug_assert!(valid(n, n, a_p, a_i) == Status::OK);

    // Construct the pointers for A+A'.

    let mut pfree: usize = 0;
    for j in 0..n {
        p_e[j] = pfree as isize;
        s_p[j] = pfree;
        pfree += len[j];
    }

    // Note that this restriction on iwlen is slightly more restrictive than
    // what is strictly required in amd_2. amd_2 can operate with no elbow
    // room at all, but it will be very slow. For better performance, at
    // least size-n elbow room is enforced.
    debug_assert!(iwlen >= pfree + n);

    #[cfg(feature = "debug1")]
    for p in 0..iwlen {
        i_w[p] = EMPTY;
    }

    for k in 0..n {
        debug1_print!("Construct row/column k = {} of A+A'\n", k);
        let p1 = a_p[k];
        let p2 = a_p[k + 1];

        // Construct A+A'.
        let mut p = p1;
        while p < p2 {
            // Scan the upper triangular part of A.
            let j = a_i[p];
            debug_assert!(/*j >= 0 &&*/ j < n);
            if j < k {
                // Entry A(j,k) in the strictly upper triangular part.
                #[cfg(feature = "debug1")]
                if j == (n - 1) {
                    debug_assert!(s_p[j] < pfree);
                } else {
                    debug_assert!((s_p[j] as isize) < p_e[j + 1]);
                }
                #[cfg(feature = "debug1")]
                if k == (n - 1) {
                    debug_assert!(s_p[k] < pfree);
                } else {
                    debug_assert!((s_p[k] as isize) < p_e[k + 1]);
                }

                i_w[s_p[j]] = k as isize;
                s_p[j] += 1;

                i_w[s_p[k]] = j as isize;
                s_p[k] += 1;

                p += 1;
            } else if j == k {
                // Skip the diagonal.
                p += 1;
                break;
            } else {
                // j > k
                // First entry below the diagonal.
                break;
            }

            // Scan lower triangular part of A, in column j until reaching
            // row k. Start where last scan left off.
            debug_assert!(a_p[j] <= t_p[j] && t_p[j] <= a_p[j + 1]);
            let mut pj = t_p[j];
            let pj2 = a_p[j + 1];
            while pj < pj2 {
                let i = a_i[pj];
                debug_assert!(/*i >= 0 &&*/ i < n);
                if i < k {
                    // A (i,j) is only in the lower part, not in upper.
                    #[cfg(feature = "debug1")]
                    if i == n - 1 {
                        debug_assert!(s_p[i] < pfree);
                    } else {
                        debug_assert!((s_p[i] as isize) < p_e[i + 1]);
                    }
                    #[cfg(feature = "debug1")]
                    if j == n - 1 {
                        debug_assert!(s_p[j] < pfree);
                    } else {
                        debug_assert!((s_p[j] as isize) < p_e[j + 1]);
                    }

                    i_w[s_p[i]] = j as isize;
                    s_p[i] += 1;

                    i_w[s_p[j]] = i as isize;
                    s_p[j] += 1;

                    pj += 1;
                } else if i == k {
                    // Entry A(k,j) in lower part and A(j,k) in upper.
                    pj += 1;
                    break;
                } else {
                    // i > k
                    // Consider this entry later, when k advances to i.
                    break;
                }
            }
            t_p[j] = pj;
        }
        t_p[k] = p;
    }

    // Clean up, for remaining mismatched entries.
    for j in 0..n {
        for pj in t_p[j]..a_p[j + 1] {
            let i = a_i[pj];

            debug_assert!(/*i >= 0 &&*/ i < n);
            // A(i,j) is only in the lower part, not in upper.
            #[cfg(feature = "debug1")]
            if i == n - 1 {
                debug_assert!(s_p[i] < pfree);
            } else {
                debug_assert!((s_p[i] as isize) < p_e[i + 1]);
            }
            #[cfg(feature = "debug1")]
            if j == n - 1 {
                debug_assert!(s_p[j] < pfree);
            } else {
                debug_assert!((s_p[j] as isize) < p_e[j + 1]);
            }

            i_w[s_p[i]] = j as isize;
            s_p[i] += 1;

            i_w[s_p[j]] = i as isize;
            s_p[j] += 1;
        }
    }

    #[cfg(feature = "debug1")]
    for j in 0..n - 1 {
        debug_assert!((s_p[j] as isize) == p_e[j + 1])
    }
    debug_assert!(s_p[n - 1] == pfree);

    // Tp and Sp no longer needed.

    // Order the matrix.
    let (_nv, p_inv, p, _e_len) = amd_2(
        n, &mut p_e, &mut i_w, &mut len, iwlen, pfree, control, &mut info,
    );

    (p, p_inv)
}
