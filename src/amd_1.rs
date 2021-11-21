use crate::amd::*;
use crate::amd_2::amd_2;
use crate::internal::*;
use crate::valid::valid;

pub fn amd_1(
    n: i32,
    a_p: &[i32],
    a_i: &[i32],
    mut len: &mut [i32],
    iwlen: i32,
    control: &Control,
    mut info: &mut Info,
) -> (Vec<i32>, Vec<i32>) {
    // Construct the matrix for amd_2.

    debug_assert!(n > 0);

    let mut p_e: Vec<isize> = vec![0; n as usize];
    let mut s_p: Vec<i32> = vec![0; n as usize];
    let mut t_p: Vec<i32> = vec![0; n as usize];
    let mut i_w: Vec<isize> = vec![0; iwlen as usize];

    debug_assert!(valid(n, n, a_p, a_i) == Status::OK);

    // Construct the pointers for A+A'.

    let mut pfree: i32 = 0;
    for j in 0..n {
        p_e[j as usize] = pfree as isize;
        s_p[j as usize] = pfree;
        pfree += len[j as usize];
    }

    // Note that this restriction on iwlen is slightly more restrictive than
    // what is strictly required in amd_2. amd_2 can operate with no elbow
    // room at all, but it will be very slow. For better performance, at
    // least size-n elbow room is enforced.
    debug_assert!(iwlen >= pfree + n);

    #[cfg(feature = "debug1")]
    for p in 0..iwlen {
        i_w[p as usize] = EMPTY;
    }

    for k in 0..n {
        debug1_print!("Construct row/column k = {} of A+A'\n", k);
        let p1 = a_p[k as usize];
        let p2 = a_p[k as usize + 1];

        // Construct A+A'.
        let mut p = p1;
        while p < p2 {
            // Scan the upper triangular part of A.
            let j = a_i[p as usize];
            debug_assert!(j >= 0 && j < n);
            if j < k {
                // Entry A(j,k) in the strictly upper triangular part.
                #[cfg(feature = "debug1")]
                if j == (n - 1) {
                    debug_assert!(s_p[j as usize] < pfree);
                } else {
                    debug_assert!(s_p[j as usize] < p_e[j as usize + 1]);
                }
                #[cfg(feature = "debug1")]
                if k == (n - 1) {
                    debug_assert!(s_p[k as usize] < pfree);
                } else {
                    debug_assert!(s_p[k as usize] < p_e[k as usize + 1]);
                }

                i_w[s_p[j as usize] as usize] = k as isize;
                s_p[j as usize] += 1;

                i_w[s_p[k as usize] as usize] = j as isize;
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
            debug_assert!(
                a_p[j as usize] <= t_p[j as usize] && t_p[j as usize] <= a_p[j as usize + 1]
            );
            let mut pj = t_p[j as usize];
            let pj2 = a_p[j as usize + 1];
            while pj < pj2 {
                let i = a_i[pj as usize];
                debug_assert!(i >= 0 && i < n);
                if i < k {
                    // A (i,j) is only in the lower part, not in upper.
                    #[cfg(feature = "debug1")]
                    if i == n - 1 {
                        debug_assert!(s_p[i as usize] < pfree);
                    } else {
                        debug_assert!(s_p[i as usize] < p_e[i as usize + 1]);
                    }
                    #[cfg(feature = "debug1")]
                    if j == n - 1 {
                        debug_assert!(s_p[j as usize] < pfree);
                    } else {
                        debug_assert!(s_p[j as usize] < p_e[j as usize + 1]);
                    }

                    i_w[s_p[i as usize] as usize] = j as isize;
                    s_p[i as usize] += 1;

                    i_w[s_p[j as usize] as usize] = i as isize;
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

            debug_assert!(i >= 0 && i < n);
            // A(i,j) is only in the lower part, not in upper.
            #[cfg(feature = "debug1")]
            if i == n - 1 {
                debug_assert!(s_p[i as usize] < pfree);
            } else {
                debug_assert!(s_p[i as usize] < p_e[i as usize + 1]);
            }
            #[cfg(feature = "debug1")]
            if j == n - 1 {
                debug_assert!(s_p[j as usize] < pfree);
            } else {
                debug_assert!(s_p[j as usize] < p_e[j as usize + 1]);
            }

            i_w[s_p[i as usize] as usize] = j as isize;
            s_p[i as usize] += 1;

            i_w[s_p[j as usize] as usize] = i as isize;
            s_p[j as usize] += 1;
        }
    }

    #[cfg(feature = "debug1")]
    for j in 0..(n - 1) as usize {
        debug_assert!(s_p[j] == p_e[j + 1])
    }
    debug_assert!(s_p[n as usize - 1] == pfree);

    // Tp and Sp no longer needed.

    // Order the matrix.
    let (_nv, p_inv, p, _e_len) = amd_2(
        n, &mut p_e, &mut i_w, &mut len, iwlen, pfree, &control, &mut info,
    );

    (p, p_inv)
}
