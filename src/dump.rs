use crate::amd::*;
use crate::internal::*;

pub fn dump(
    n: i32,
    Pe: &[i32],
    Iw: &[i32],
    Len: &[i32],
    iwlen: i32,
    pfree: i32,
    Nv: &[i32],
    Next: &[i32],
    Last: &[i32],
    Head: &[i32],
    Elen: &[i32],
    degree: &[i32],
    W: &[i32],
    nel: i32,
) {
    if DEBUG_LEVEL < 3 {
        return;
    }

    assert!(pfree <= iwlen);
    print!("\nAMD dump, pfree: {}\n", pfree);
    for i in 0..n {
        let pe = Pe[i];
        let elen = Elen[i];
        let nv = Nv[i];
        let len = Len[i];
        let w = W[i];

        if elen >= EMPTY {
            if nv == 0 {
                print!("\nI {}: nonprincipal:    \n", i);
                assert!(elen == EMPTY);
                if pe == EMPTY {
                    println!(" dense node");
                    assert!(w == 1);
                } else {
                    assert!(pe < EMPTY);
                    print!(" i {} -> parent {}\n", i, flip(Pe[i]));
                }
            } else {
                print!("\nI {}: active principal supervariable:\n", i);
                print!("   nv(i): {}  Flag: {}\n", nv, nv < 0);
                assert!(elen >= 0);
                assert!(nv > 0 && pe >= 0);
                let mut p = pe;
                print!("   e/s: ");
                if elen == 0 {
                    print!(" : ");
                }
                assert!(pe + len <= pfree);
                for k in 0..len {
                    let j = Iw[p];
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
            if w == 0 {
                print!("\nE {}: absorbed element: w {}\n", e, w);
                assert!(nv > 0 && pe < 0);
                print!(" e {} -> parent {}\n", e, flip(Pe[e]));
            } else {
                print!("\nE {}: unabsorbed element: w {}\n", e, w);
                assert!(nv > 0 && pe >= 0);
                let mut p = pe;
                print!(" : ");
                assert!(pe + len <= pfree);
                for k in 0..len {
                    let j = Iw[p];
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
            if Head[deg] == EMPTY {
                continue;
            }
            let mut ilast = EMPTY;
            print!("{}: \n", deg);
            let mut i = Head[deg];
            while i != EMPTY {
                print!(
                    "   {} : next {} last {} deg {}\n",
                    i, Next[i], Last[i], degree[i]
                );
                assert!(i >= 0 && i < n && ilast == Last[i] && deg == degree[i]);
                cnt += Nv[i];
                ilast = i;

                i = Next[i];
            }
            //debug3("\n")
        }
        assert!(cnt == n - nel);
    }
}
