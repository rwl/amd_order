use crate::amd::*;
use crate::internal::*;
use crate::valid::valid;
use std::cmp::{max, min};

fn clear_flag(wflg: i32, wbig: i32, w: &mut [i32], n: i32) -> i32 {
    if wflg < 2 || wflg >= wbig {
        for x in 0..n as usize {
            if w[x] != 0 {
                w[x] = 1;
            }
        }
        return 2;
    }
    // At this point, W[0..n-1] < wflg holds.
    return wflg;
}

pub fn amd_2(
    n: i32,
    Pe: &mut [i32],
    Iw: &[i32],
    Len: &[i32],
    iwlen: i32,
    pfree: i32,
    Nv: &mut [i32],
    Next: &mut [i32],
    Last: &mut [i32],
    Head: &mut [i32],
    Elen: &mut [i32],
    degree: &mut [i32],
    W: &mut [i32],
    control: Control,
    info: &mut Info,
) {
    let hash: u32; // unsigned, so that hash % n is well defined.

    // Any parameter (Pe[...] or pfree) or local variable starting with "p" (for
    // Pointer) is an index into Iw, and all indices into Iw use variables starting
    // with "p". The only exception to this rule is the iwlen input argument.

    // Initializations

    // Note that this restriction on iwlen is slightly more restrictive than
    // what is actually required in amd_2. amd_2 can operate with no elbow
    // room at all, but it will be slow. For better performance, at least
    // size-n elbow room is enforced.
    if DEBUG_LEVEL != 0 {
        debug_assert!(iwlen >= pfree + n);
        debug_assert!(n > 0);
    }

    // Initialize output statistics.
    let mut lnz = 0.0; // The number of nonzeros in L (excluding the diagonal).
    let mut ndiv = 0.0; // Number of divisions for LU or LDL' factorizations.
    let mut nms_lu = 0.0; // Number of multiply-subtract pairs for LU factorization.
    let mut nms_ldl = 0.0; // Number of multiply-subtract pairs for LDL' factorization.
    let mut dmax = 1.0; // The largest number of entries in any column of L, including the diagonal.
                        // Current supervariable being eliminated, and the current
                        // element created by eliminating that supervariable.
    let mut me = EMPTY;

    let mut mindeg = 0; // Current minimum degree.
    let mut ncmpa = 0; // Number of garbage collections.
    let mut nel = 0; // Number of pivots selected so far.
    let mut lemax = 0; // Largest |Le| seen so far (called dmax in Fortran version).

    // Get control parameters.
    let aggressive = if control.aggressive { 1 } else { 0 };
    // "dense" degree ratio.
    let alpha = control.dense;
    // Note: if alpha is NaN, this is undefined:
    let mut dense = if alpha < 0.0 {
        // Only remove completely dense rows/columns.
        (n - 2)
    } else {
        (alpha * (n as f64).sqrt()) as i32
    };
    dense = max(16, dense);
    let dense = min(n, dense);
    if DEBUG_LEVEL >= 1 {
        print!("\n\nAMD (debug), alpha {}, aggr. {}\n", alpha, aggressive);
    }

    for i in 0..n as usize {
        Last[i] = EMPTY;
        Head[i] = EMPTY;
        Next[i] = EMPTY;
        // if separate Hhead array is used for hash buckets:
        //   Hhead[i] = EMPTY
        Nv[i] = 1;
        W[i] = 1;
        Elen[i] = 0;
        degree[i] = Len[i];
    }

    if DEBUG_LEVEL >= 1 {
        print!("\n======Nel {} initial\n", nel)
        // TODO: dump(n, Pe, Iw, Len, iwlen, pfree, Nv, Next, Last, Head, Elen, degree, W, -1);
    }

    // INT_MAX - n for the int version, UF_long_max - n for the
    // int64 version. wflg is not allowed to be >= wbig.
    let wbig = std::i32::MAX - n;
    // Used for flagging the W array. See description of Iw.
    let wflg = clear_flag(0, wbig, W, n);

    // Initialize degree lists and eliminate dense and empty rows.

    let mut ndense = 0; // Number of "dense" rows/columns.

    for i in 0..n as usize {
        let deg = degree[i]; // The degree of a variable or element.
        if DEBUG_LEVEL != 0 {
            debug_assert!(deg >= 0 && deg < n);
        }
        if deg == 0 {
            // We have a variable that can be eliminated at once because
            // there is no off-diagonal non-zero in its row. Note that
            // Nv [i] = 1 for an empty variable i. It is treated just
            // the same as an eliminated element i.

            Elen[i] = flip(1);
            nel += 1;
            Pe[i] = EMPTY;
            W[i] = 0;
        } else if deg > dense {
            // Dense variables are not treated as elements, but as unordered,
            // non-principal variables that have no parent. They do not take
            // part in the postorder, since Nv [i] = 0. Note that the Fortran
            // version does not have this option.

            if DEBUG_LEVEL >= 1 {
                print!("Dense node {} degree {}\n", i, deg);
            }
            ndense += 1;
            Nv[i] = 0;
            // do not postorder this node
            Elen[i] = EMPTY;
            nel += 1;
            Pe[i] = EMPTY
        } else {
            // Place i in the degree list corresponding to its degree.

            let inext = Head[deg as usize]; // The entry in a link list following i.
            if DEBUG_LEVEL != 0 {
                debug_assert!(inext >= EMPTY && inext < n);
            }
            if inext != EMPTY {
                Last[inext as usize] = i as i32;
            }
            Next[i] = inext;
            Head[deg as usize] = i as i32;
        }
    }
}
