#[cfg(feature = "debug1")]
use crate::internal::*;

#[cfg(feature = "debug1")]
pub fn dump(
    n: i32,
    pe: &[i32],
    iw: &[i32],
    len: &[i32],
    iwlen: i32,
    pfree: i32,
    nv: &[i32],
    next: &[i32],
    last: &[i32],
    head: &[i32],
    e_len: &[i32],
    degree: &[i32],
    w: &[i32],
    nel: i32,
) {
    debug_assert!(pfree <= iwlen);
    debug3_print!("\nAMD dump, pfree: {}\n", pfree);
    for i in 0..n {
        let pei = pe[i as usize];
        let elen = e_len[i as usize];
        let nvi = nv[i as usize];
        let leni = len[i as usize];
        let wi = w[i as usize];

        if elen >= EMPTY {
            if nvi == 0 {
                debug3_print!("\nI {}: nonprincipal:    \n", i);
                debug_assert!(elen == EMPTY);
                if pei == EMPTY {
                    debug3_println!(" dense node");
                    debug_assert!(wi == 1);
                } else {
                    debug_assert!(pei < EMPTY);
                    debug3_print!(" i {} -> parent {}\n", i, flip(pe[i as usize]));
                }
            } else {
                debug3_print!("\nI {}: active principal supervariable:\n", i);
                debug3_print!("   nv(i): {}  Flag: {}\n", nvi, nvi < 0);
                debug_assert!(elen >= 0);
                debug_assert!(nvi > 0 && pei >= 0);
                let mut p = pei;
                debug3_print!("   e/s: ");
                if elen == 0 {
                    debug3_print!(" : ");
                }
                debug_assert!(pei + leni <= pfree);
                for k in 0..leni {
                    let j = iw[p as usize];
                    debug3_print!("  {}", j);
                    debug_assert!(j >= 0 && j < n);
                    if k == elen - 1 {
                        debug3_print!(" : ");
                    }
                    p += 1;
                }
                debug3_println!();
            }
        } else {
            #[cfg(feature = "debug3")]
            let e = i;
            if wi == 0 {
                debug3_print!("\nE {}: absorbed element: w {}\n", e, wi);
                debug_assert!(nvi > 0 && pei < 0);
                debug3_print!(" e {} -> parent {}\n", e, flip(pe[e as usize]));
            } else {
                debug3_print!("\nE {}: unabsorbed element: w {}\n", e, wi);
                debug_assert!(nvi > 0 && pei >= 0);
                let mut p = pei;
                debug3_print!(" : ");
                debug_assert!(pei + leni <= pfree);
                for _k in 0..leni {
                    let j = iw[p as usize];
                    debug3_print!("  {}", j);
                    debug_assert!(j >= 0 && j < n);
                    p += 1;
                }
                debug3_println!();
            }
        }
    }

    // This routine cannot be called when the hash buckets are non-empty.
    debug3_println!("\nDegree lists:");
    if nel >= 0 {
        let mut cnt = 0;
        for deg in 0..n {
            if head[deg as usize] == EMPTY {
                continue;
            }
            let mut ilast = EMPTY;
            debug3_print!("{}: \n", deg);
            let mut i = head[deg as usize];
            while i != EMPTY {
                debug3_print!(
                    "   {} : next {} last {} deg {}\n",
                    i,
                    next[i as usize],
                    last[i as usize],
                    degree[i as usize]
                );
                debug_assert!(
                    i >= 0 && i < n && ilast == last[i as usize] && deg == degree[i as usize]
                );
                cnt += nv[i as usize];
                ilast = i;

                i = next[i as usize];
            }
            debug3_print!("\n");
        }
        debug_assert!(cnt == n - nel);
    }
}

// #[cfg(not(feature = "debug1"))]
// pub fn dump(
//     _n: i32,
//     _pe: &[i32],
//     _iw: &[i32],
//     _len: &[i32],
//     _iwlen: i32,
//     _pfree: i32,
//     _nv: &[i32],
//     _next: &[i32],
//     _last: &[i32],
//     _head: &[i32],
//     _e_len: &[i32],
//     _degree: &[i32],
//     _w: &[i32],
//     _nel: i32,
// ) {
// }
