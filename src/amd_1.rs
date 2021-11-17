use crate::amd::*;
use crate::amd_2::amd_2;
use crate::internal::EMPTY;
use crate::valid::valid;

pub fn amd_1(
    n: i32,
    a_p: &[i32],
    a_i: &[i32],
    p: &mut [i32],
    p_inv: &mut [i32],
    mut len: &mut [i32],
    slen: i32,
    s: &mut [i32],
    // mut s: &Vec<i32>,
    control: Control,
    mut info: &mut Info,
) {
    // Construct the matrix for amd_2.

    debug_assert!(n > 0);

    let iwlen = slen - 6 * n;
    let (p_e, s) = s.split_at_mut(n as usize);
    let (n_v, s) = s.split_at_mut(n as usize);
    let (head, s) = s.split_at_mut(n as usize);
    let (e_len, s) = s.split_at_mut(n as usize);
    let (degree, s) = s.split_at_mut(n as usize);
    let (w, s) = s.split_at_mut(n as usize);
    let (i_w, _s) = s.split_at_mut(iwlen as usize);
    // let mut p_e: Vec<i32> = vec![0; n];
    // let mut n_v: Vec<i32> = vec![0; n];
    // let mut head: Vec<i32> = vec![0; n];
    // let mut e_len: Vec<i32> = vec![0; n];
    // let mut degree: Vec<i32> = vec![0; n];
    // let mut w: Vec<i32> = vec![0; n];
    // let mut i_w: Vec<i32> = vec![0; iwlen];

    debug_assert!(valid(n, n, a_p, a_i) == Status::OK);

    // Construct the pointers for A+A'.

    let s_p = n_v; // Use Nv and W as workspace for Sp and Tp.
    let t_p = w;
    let mut pfree: i32 = 0;
    for j in 0..n {
        p_e[j as usize] = pfree;
        s_p[j as usize] = pfree;
        pfree += len[j as usize];
    }

    // Note that this restriction on iwlen is slightly more restrictive than
    // what is strictly required in amd_2. amd_2 can operate with no elbow
    // room at all, but it will be very slow. For better performance, at
    // least size-n elbow room is enforced.
    debug_assert!(iwlen >= pfree + n);

    if DEBUG_LEVEL != 0 {
        for p in 0..iwlen {
            i_w[p as usize] = EMPTY;
        }
    }

    for k in 0..n {
        if DEBUG_LEVEL >= 1 {
            print!("Construct row/column k = {} of A+A'\n", k);
        }
        // var p, pj int
        let p1 = a_p[k as usize];
        let p2 = a_p[k as usize + 1];

        // Construct A+A'.
        let mut p = p1;
        while p < p2 {
            // Scan the upper triangular part of A.
            let j = a_i[p as usize];
            if DEBUG_LEVEL != 0 {
                debug_assert!(j >= 0 && j < n);
            }
            if j < k {
                if DEBUG_LEVEL != 0 {
                    // Entry A(j,k) in the strictly upper triangular part.
                    if j == (n - 1) {
                        debug_assert!(s_p[j as usize] < pfree);
                    } else {
                        debug_assert!(s_p[j as usize] < p_e[j as usize + 1]);
                    }
                    if k == (n - 1) {
                        debug_assert!(s_p[k as usize] < pfree);
                    } else {
                        debug_assert!(s_p[k as usize] < p_e[k as usize + 1]);
                    }
                }
                i_w[s_p[j as usize] as usize] = k;
                s_p[j as usize] += 1;

                i_w[s_p[k as usize] as usize] = j;
                s_p[k as usize] += 1;

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
            if DEBUG_LEVEL != 0 {
                debug_assert!(
                    a_p[j as usize] <= t_p[j as usize] && t_p[j as usize] <= a_p[j as usize + 1]
                );
            }
            let mut pj = t_p[j as usize];
            let pj2 = a_p[j as usize + 1];
            while pj < pj2 {
                let i = a_i[pj as usize];
                if DEBUG_LEVEL != 0 {
                    debug_assert!(i >= 0 && i < n);
                }
                if i < k {
                    if DEBUG_LEVEL != 0 {
                        // A (i,j) is only in the lower part, not in upper.
                        if i == n - 1 {
                            debug_assert!(s_p[i as usize] < pfree);
                        } else {
                            debug_assert!(s_p[i as usize] < p_e[i as usize + 1]);
                        }
                        if j == n - 1 {
                            debug_assert!(s_p[j as usize] < pfree);
                        } else {
                            debug_assert!(s_p[j as usize] < p_e[j as usize + 1]);
                        }
                    }
                    i_w[s_p[i as usize] as usize] = j;
                    s_p[i as usize] += 1;

                    i_w[s_p[j as usize] as usize] = i;
                    s_p[j as usize] += 1;

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
            t_p[j as usize] = pj;
        }
        t_p[k as usize] = p;
    }

    // Clean up, for remaining mismatched entries.
    for j in 0..n {
        for pj in t_p[j as usize]..a_p[j as usize + 1] {
            let i = a_i[pj as usize];

            if DEBUG_LEVEL != 0 {
                debug_assert!(i >= 0 && i < n);
                // A(i,j) is only in the lower part, not in upper.
                if i == n - 1 {
                    debug_assert!(s_p[i as usize] < pfree);
                } else {
                    debug_assert!(s_p[i as usize] < p_e[i as usize + 1]);
                }
                if j == n - 1 {
                    debug_assert!(s_p[j as usize] < pfree);
                } else {
                    debug_assert!(s_p[j as usize] < p_e[j as usize + 1]);
                }
            }

            i_w[s_p[i as usize] as usize] = j;
            s_p[i as usize] += 1;

            i_w[s_p[j as usize] as usize] = i;
            s_p[j as usize] += 1;
        }
    }

    if DEBUG_LEVEL != 0 {
        for j in 0..(n - 1) as usize {
            debug_assert!(s_p[j] == p_e[j + 1])
        }
        debug_assert!(s_p[n as usize - 1] == pfree);
    }

    // Tp and Sp no longer needed.

    let mut n_v = s_p;
    let w = t_p;

    // Order the matrix.

    amd_2(
        n, p_e, i_w, &mut len, iwlen, pfree, &mut n_v, p_inv, p, head, e_len, degree, w, control,
        &mut info,
    )
}
