use crate::amd::*;
use crate::internal::*;

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
    if DEBUG_LEVEL < 3 {
        return;
    }

    assert!(pfree <= iwlen);
    print!("\nAMD dump, pfree: {}\n", pfree);
    for i in 0..n {
        let pei = pe[i as usize];
        let elen = e_len[i as usize];
        let nvi = nv[i as usize];
        let leni = len[i as usize];
        let wi = w[i as usize];

        if elen >= EMPTY {
            if nvi == 0 {
                print!("\nI {}: nonprincipal:    \n", i);
                assert!(elen == EMPTY);
                if pei == EMPTY {
                    println!(" dense node");
                    assert!(wi == 1);
                } else {
                    assert!(pei < EMPTY);
                    print!(" i {} -> parent {}\n", i, flip(pe[i as usize]));
                }
            } else {
                print!("\nI {}: active principal supervariable:\n", i);
                print!("   nv(i): {}  Flag: {}\n", nvi, nvi < 0);
                assert!(elen >= 0);
                assert!(nvi > 0 && pei >= 0);
                let mut p = pei;
                print!("   e/s: ");
                if elen == 0 {
                    print!(" : ");
                }
                assert!(pei + leni <= pfree);
                for k in 0..leni {
                    let j = iw[p as usize];
                    print!("  {}", j);
                    assert!(j >= 0 && j < n);
                    if k == elen - 1 {
                        print!(" : ");
                    }
                    p += 1;
                }
                println!();
            }
        } else {
            let e = i;
            if wi == 0 {
                print!("\nE {}: absorbed element: w {}\n", e, wi);
                assert!(nvi > 0 && pei < 0);
                print!(" e {} -> parent {}\n", e, flip(pe[e as usize]));
            } else {
                print!("\nE {}: unabsorbed element: w {}\n", e, wi);
                assert!(nvi > 0 && pei >= 0);
                let mut p = pei;
                print!(" : ");
                assert!(pei + leni <= pfree);
                for _k in 0..leni {
                    let j = iw[p as usize];
                    print!("  {}", j);
                    assert!(j >= 0 && j < n);
                    p += 1;
                }
                println!();
            }
        }
    }

    // This routine cannot be called when the hash buckets are non-empty.
    println!("\nDegree lists:");
    if nel >= 0 {
        let mut cnt = 0;
        for deg in 0..n {
            if head[deg as usize] == EMPTY {
                continue;
            }
            let mut ilast = EMPTY;
            print!("{}: \n", deg);
            let mut i = head[deg as usize];
            while i != EMPTY {
                print!(
                    "   {} : next {} last {} deg {}\n",
                    i, next[i as usize], last[i as usize], degree[i as usize]
                );
                assert!(i >= 0 && i < n && ilast == last[i as usize] && deg == degree[i as usize]);
                cnt += nv[i as usize];
                ilast = i;

                i = next[i as usize];
            }
            //debug3("\n")
        }
        assert!(cnt == n - nel);
    }
}
