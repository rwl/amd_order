use crate::amd::*;
use crate::dump::dump;
use crate::internal::*;
use crate::postorder::postorder;
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
    mut pfree: i32,
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
    // let mut pfree = pfree_;
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
        print!("\n======Nel {} initial\n", nel);
        dump(
            n, Pe, Iw, Len, iwlen, pfree, Nv, Next, Last, Head, Elen, degree, W, -1,
        );
    }

    // INT_MAX - n for the int version, UF_long_max - n for the
    // int64 version. wflg is not allowed to be >= wbig.
    let wbig = std::i32::MAX - n;
    // Used for flagging the W array. See description of Iw.
    let mut wflg = clear_flag(0, wbig, W, n);

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

    // While (selecting pivots) do.
    while nel < n {
        if DEBUG_LEVEL >= 1 {
            print!("\n======Nel {}\n", nel);
            if DEBUG_LEVEL >= 2 {
                dump(
                    n, Pe, Iw, Len, iwlen, pfree, Nv, Next, Last, Head, Elen, degree, W, nel,
                );
            }
        }

        // Get pivot of minimum degree.

        // Find next supervariable for elimination.
        if DEBUG_LEVEL != 0 {
            debug_assert!(mindeg >= 0 && mindeg < n);
        }
        let mut deg = mindeg; // The degree of a variable or element.
        while deg < n {
            me = Head[deg as usize];
            if me != EMPTY {
                break;
            }
            deg += 1;
        }
        mindeg = deg;
        if DEBUG_LEVEL != 0 {
            debug_assert!(me >= 0 && me < n);
            if DEBUG_LEVEL >= 1 {
                print!("=================me: {}\n", me);
            }
        }

        // Remove chosen variable from link list.
        let inext = Next[me]; // The entry in a link list following i.
        if DEBUG_LEVEL != 0 {
            debug_assert!(inext >= EMPTY && inext < n);
        }
        if inext != EMPTY {
            Last[inext] = EMPTY;
        }
        Head[deg] = inext;

        // me represents the elimination of pivots nel to nel+Nv[me]-1.
        // place me itself as the first in this set.
        let elenme = Elen[me]; // The length, Elen [me], of element list of pivotal variable.
        let nvpiv = Nv[me]; // Number of pivots in current element.
        if DEBUG_LEVEL != 0 {
            debug_assert!(nvpiv > 0);
        }
        nel += nvpiv;

        // Construct new element.

        // At this point, me is the pivotal supervariable. It will be
        // converted into the current element. Scan list of the pivotal
        // supervariable, me, setting tree pointers and constructing new list
        // of supervariables for the new element, me. p is a pointer to the
        // current position in the old list.

        // Flag the variable "me" as being in Lme by negating Nv[me].
        Nv[me] = -nvpiv;
        let mut degme = 0; // Size, |Lme|, of the current element, me (= degree[me]).
        if DEBUG_LEVEL != 0 {
            debug_assert!(Pe[me] >= 0 && Pe[me] < iwlen);
        }

        let mut pme1: i32; // The current element, me, is stored in Iw[pme1...pme2].
        let mut pme2: i32; // The end of the current element.
        if elenme == 0 {
            // Construct the new element in place.
            pme1 = Pe[me];
            pme2 = pme1 - 1;

            for p in pme1..=pme1 + Len[me] - 1 {
                let i = Iw[p];
                if DEBUG_LEVEL != 0 {
                    debug_assert!(i >= 0 && i < n && Nv[i] >= 0);
                }
                let nvi = Nv[i]; // The number of variables in a supervariable i (= Nv[i])
                if nvi > 0 {
                    // i is a principal variable not yet placed in Lme.
                    // store i in new list

                    // Flag i as being in Lme by negating Nv[i].
                    degme += nvi;
                    Nv[i] = -nvi;
                    pme2 += 1;
                    Iw[pme2] = i;

                    // Remove variable i from degree list.
                    let ilast = Last[i]; // The entry in a link list preceding i.
                    inext = Next[i];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(ilast >= EMPTY && ilast < n);
                        debug_assert!(inext >= EMPTY && inext < n);
                    }
                    if inext != EMPTY {
                        Last[inext] = ilast;
                    }
                    if ilast != EMPTY {
                        Next[ilast] = inext;
                    } else {
                        // i is at the head of the degree list.
                        if DEBUG_LEVEL != 0 {
                            debug_assert!(degree[i] >= 0 && degree[i] < n);
                        }
                        Head[degree[i]] = inext;
                    }
                }
            }
        } else {
            // Construct the new element in empty space, Iw[pfree ...]
            let p = Pe[me];
            pme1 = pfree;
            // Number of variables in variable list of pivotal variable.
            let slenme = Len[me] - elenme;

            for knt1 in 1..=elenme + 1 {
                let e: i32;
                let pj: i32;
                let ln: i32;
                if knt1 > elenme {
                    // Search the supervariables in me.
                    e = me;
                    pj = p;
                    ln = slenme;
                    if DEBUG_LEVEL >= 2 {
                        print!("Search sv: {} {} {}\n", me, pj, ln);
                    }
                } else {
                    // Search the elements in me.
                    e = Iw[p];
                    p += 1;
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(e >= 0 && e < n);
                    }
                    pj = Pe[e];
                    ln = Len[e];
                    if DEBUG_LEVEL != 0 {
                        if DEBUG_LEVEL >= 2 {
                            print!("Search element e {} in me {}\n", e, me);
                        }
                        debug_assert!(Elen[e] < EMPTY && W[e] > 0 && pj >= 0);
                    }
                }
                if DEBUG_LEVEL != 0 {
                    debug_assert!(ln >= 0 && (ln == 0 || (pj >= 0 && pj < iwlen)));
                }

                // search for different supervariables and add them to the
                // new list, compressing when necessary. this loop is
                // executed once for each element in the list and once for
                // all the supervariables in the list.

                for knt2 in 1..=ln {
                    let i = Iw[pj];
                    pj += 1;
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(i >= 0 && i < n && (i == me || Elen[i] >= EMPTY));
                    }
                    // The number of variables in a supervariable i (= Nv[i]).
                    let nvi = Nv[i];
                    if DEBUG_LEVEL >= 2 {
                        print!(": {} {} {} {}\n", i, Elen[i], Nv[i], wflg);
                    }

                    if nvi > 0 {
                        // Compress Iw, if necessary.
                        if pfree >= iwlen {
                            if DEBUG_LEVEL >= 1 {
                                println!("GARBAGE COLLECTION");
                            }

                            // Prepare for compressing Iw by adjusting pointers
                            // and lengths so that the lists being searched in
                            // the inner and outer loops contain only the
                            // remaining entries.

                            Pe[me] = p;
                            Len[me] -= knt1;
                            // Check if nothing left of supervariable me.
                            if Len[me] == 0 {
                                Pe[me] = EMPTY;
                            }
                            Pe[e] = pj;
                            Len[e] = ln - knt2;
                            // Nothing left of element e.
                            if Len[e] == 0 {
                                Pe[e] = EMPTY;
                            }

                            ncmpa += 1; // One more garbage collection.

                            // Store first entry of each object in Pe
                            // flip the first entry in each object
                            for j in 0..n {
                                let pn = Pe[j];
                                if pn >= 0 {
                                    if DEBUG_LEVEL != 0 {
                                        debug_assert!(pn >= 0 && pn < iwlen);
                                    }
                                    Pe[j] = Iw[pn];
                                    Iw[pn] = flip(j);
                                }
                            }

                            // psrc/pdst point to source/destination
                            let mut psrc = 0;
                            let mut pdst = 0;
                            let pend = pme1 - 1;

                            while psrc <= pend {
                                // Search for next flip'd entry.
                                let j = flip(Iw[psrc]);
                                psrc += 1;
                                if j >= 0 {
                                    if DEBUG_LEVEL >= 2 {
                                        print!("Got object j: {}\n", j);
                                    }
                                    Iw[pdst] = Pe[j];
                                    Pe[j] = pdst;
                                    pdst += 1;
                                    let lenj = Len[j];
                                    // Copy from source to destination.
                                    for knt3 in 0..=lenj - 2 {
                                        Iw[pdst] = Iw[psrc];
                                        pdst += 1;
                                        psrc += 1;
                                    }
                                }
                            }

                            // Move the new partially-constructed element.
                            let p1 = pdst;
                            psrc = pme1;
                            while psrc <= pfree - 1 {
                                Iw[pdst] = Iw[psrc];
                                pdst += 1;
                                psrc += 1;
                            }
                            pme1 = p1;
                            pfree = pdst;
                            pj = Pe[e];
                            p = Pe[me];
                        }

                        // i is a principal variable not yet placed in Lme
                        // store i in new list.

                        // Flag i as being in Lme by negating Nv[i].
                        degme += nvi;
                        Nv[i] = -nvi;
                        Iw[pfree] = i;
                        pfree += 1;
                        if DEBUG_LEVEL >= 2 {
                            print!("     s: {}     nv {}\n", i, Nv[i]);
                        }

                        // Remove variable i from degree link list.

                        let ilast = Last[i]; // The entry in a link list preceding i.
                        inext = Next[i];
                        if DEBUG_LEVEL != 0 {
                            debug_assert!(ilast >= EMPTY && ilast < n);
                            debug_assert!(inext >= EMPTY && inext < n);
                        }
                        if inext != EMPTY {
                            Last[inext] = ilast;
                        }
                        if ilast != EMPTY {
                            Next[ilast] = inext;
                        } else {
                            // i is at the head of the degree list.
                            if DEBUG_LEVEL != 0 {
                                debug_assert!(degree[i] >= 0 && degree[i] < n);
                            }
                            Head[degree[i]] = inext;
                        }
                    }
                }

                if e != me {
                    // Set tree pointer and flag to indicate element e is
                    // absorbed into new element me (the parent of e is me).
                    if DEBUG_LEVEL >= 1 {
                        print!(" Element {} => {}\n", e, me);
                    }
                    Pe[e] = flip(me);
                    W[e] = 0;
                }
            }

            pme2 = pfree - 1;
        }

        // me has now been converted into an element in Iw[pme1..pme2]

        // degme holds the external degree of new element.
        degree[me] = degme;
        Pe[me] = pme1;
        Len[me] = pme2 - pme1 + 1;
        if DEBUG_LEVEL != 0 {
            debug_assert!(Pe[me] >= 0 && Pe[me] < iwlen);
        }

        Elen[me] = flip(nvpiv + degme);
        // flip(Elen(me)) is now the degree of pivot (including diagonal part).

        if DEBUG_LEVEL >= 2 {
            print!("New element structure: length={}\n", pme2 - pme1 + 1);
            if DEBUG_LEVEL >= 3 {
                for pme in pme1..=pme2 {
                    print!(" {}", Iw[pme]);
                }
                println!();
            }
        }

        // Make sure that wflg is not too large.

        // With the current value of wflg, wflg+n must not cause integer overflow.

        wflg = clear_flag(wflg, wbig, W, n);

        // compute(W [e] - wflg) = |Le\Lme| for all elements.

        // Scan 1:  compute the external degrees of previous elements with
        // respect to the current element. That is:
        //       (W [e] - wflg) = |Le \ Lme|
        // for each element e that appears in any supervariable in Lme. The
        // notation Le refers to the pattern (list of supervariables) of a
        // previous element e, where e is not yet absorbed, stored in
        // Iw [Pe [e] + 1 ... Pe [e] + Len [e]]. The notation Lme
        // refers to the pattern of the current element (stored in
        // Iw [pme1..pme2]).  If aggressive absorption is enabled, and
        // (W [e] - wflg) becomes zero, then the element e will be absorbed
        // in Scan 2.

        if DEBUG_LEVEL >= 2 {
            print!("me: ");
        }
        for pme in pme1..=pme2 {
            let i = Iw[pme];
            if DEBUG_LEVEL != 0 {
                debug_assert!(i >= 0 && i < n);
            }
            let eln = Elen[i]; // The length, Elen[...], of an element list.
            if DEBUG_LEVEL >= 3 {
                print!("{} Elen {}: \n", i, eln);
            }
            if eln > 0 {
                // Note that Nv[i] has been negated to denote i in Lme:
                let nvi = -Nv[i];
                if DEBUG_LEVEL != 0 {
                    debug_assert!(nvi > 0 && Pe[i] >= 0 && Pe[i] < iwlen);
                }
                let wnvi = wflg - nvi;
                for p in Pe[i]..=Pe[i] + eln - 1 {
                    let e = Iw[p];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(e >= 0 && e < n);
                    }
                    let we = W[e];
                    if DEBUG_LEVEL >= 4 {
                        print!("    e {} we {} ", e, we);
                    }
                    if we >= wflg {
                        // Unabsorbed element e has been seen in this loop.
                        if DEBUG_LEVEL >= 4 {
                            print!("    unabsorbed, first time seen");
                        }
                        we -= nvi;
                    } else if we != 0 {
                        // e is an unabsorbed element.
                        // This is the first we have seen e in all of Scan 1.
                        if DEBUG_LEVEL >= 4 {
                            print!("    unabsorbed");
                        }
                        we = degree[e] + wnvi;
                    }
                    if DEBUG_LEVEL >= 4 {
                        println!();
                    }
                    W[e] = we;
                }
            }
        }
        if DEBUG_LEVEL >= 3 {
            println!();
        }

        // Degree update and element absorption.

        // Scan 2:  for each i in Lme, sum up the degree of Lme (which is
        // degme), plus the sum of the external degrees of each Le for the
        // elements e appearing within i, plus the supervariables in i.
        // Place i in hash list.

        for pme in pme1..=pme2 {
            let i = Iw[pme];
            if DEBUG_LEVEL != 0 {
                debug_assert!(i >= 0 && i < n && Nv[i] < 0 && Elen[i] >= 0);
                if DEBUG_LEVEL >= 2 {
                    print!("Updating: i {} {} {}\n", i, Elen[i], Len[i]);
                }
            }
            let p1 = Pe[i];
            let p2 = p1 + Elen[i] - 1;
            let pn = p1;
            hash = 0;
            deg = 0;
            if DEBUG_LEVEL != 0 {
                debug_assert!(p1 >= 0 && p1 < iwlen && p2 >= -1 && p2 < iwlen);
            }

            // scan the element list associated with supervariable i .

            // UMFPACK/MA38-style approximate degree:
            if aggressive != 0 {
                for p in p1..=p2 {
                    let e = Iw[p];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(e >= 0 && e < n);
                    }
                    let we = W[e];
                    if we != 0 {
                        // e is an unabsorbed element.
                        let dext = we - wflg; // External degree, |Le \ Lme|, of some element e.
                        if dext > 0 {
                            deg += dext;
                            Iw[pn] = e;
                            pn += 1;
                            hash += e as u32;
                            if DEBUG_LEVEL >= 4 {
                                print!(" e: {} hash = {}\n", e, hash);
                            }
                        } else {
                            // External degree of e is zero, absorb e into me.
                            if DEBUG_LEVEL != 0 {
                                print!(" Element {} => {} (aggressive)\n", e, me);
                                debug_assert!(dext == 0);
                            }
                            Pe[e] = flip(me);
                            W[e] = 0;
                        }
                    }
                }
            } else {
                for p in p1..=p2 {
                    let e = Iw[p];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(e >= 0 && e < n);
                    }
                    let we = W[e];
                    if we != 0 {
                        // e is an unabsorbed element.
                        let dext = we - wflg;
                        if DEBUG_LEVEL != 0 {
                            debug_assert!(dext >= 0);
                        }
                        deg += dext;
                        Iw[pn] = e;
                        pn += 1;
                        hash += e as u32;
                        if DEBUG_LEVEL >= 4 {
                            print!(" e: {} hash = {}\n", e, hash);
                        }
                    }
                }
            }

            // Count the number of elements in i (including me):
            Elen[i] = pn - p1 + 1;

            // Scan the supervariables in the list associated with i.

            // The bulk of the AMD run time is typically spent in this loop,
            // particularly if the matrix has many dense rows that are not
            // removed prior to ordering.
            let p3 = pn;
            let p4 = p1 + Len[i];
            for p in p2 + 1..p4 {
                let j = Iw[p];
                if DEBUG_LEVEL != 0 {
                    debug_assert!(j >= 0 && j < n);
                }
                let nvj = Nv[j];
                if nvj > 0 {
                    // j is unabsorbed, and not in Lme.
                    // Add to degree and add to new list.
                    deg += nvj;
                    Iw[pn] = j;
                    pn += 1;
                    hash += j as u32;
                    if DEBUG_LEVEL >= 4 {
                        print!("  s: {} hash {} Nv[j]= {}\n", j, hash, nvj);
                    }
                }
            }

            // Update the degree and check for mass elimination.

            if DEBUG_LEVEL != 0 {
                // With aggressive absorption, deg==0 is identical to the
                // Elen [i] == 1 && p3 == pn test, below.
                debug_assert!(implies(
                    aggressive != 0,
                    (deg == 0) == (Elen[i] == 1 && p3 == pn)
                ));
            }

            if Elen[i] == 1 && p3 == pn {
                // Mass elimination

                // There is nothing left of this node except for an edge to
                // the current pivot element. Elen [i] is 1, and there are
                // no variables adjacent to node i. Absorb i into the
                // current pivot element, me. Note that if there are two or
                // more mass eliminations, fillin due to mass elimination is
                // possible within the nvpiv-by-nvpiv pivot block. It is this
                // step that causes AMD's analysis to be an upper bound.
                //
                // The reason is that the selected pivot has a lower
                // approximate degree than the true degree of the two mass
                // eliminated nodes. There is no edge between the two mass
                // eliminated nodes. They are merged with the current pivot
                // anyway.
                //
                // No fillin occurs in the Schur complement, in any case,
                // and this effect does not decrease the quality of the
                // ordering itself, just the quality of the nonzero and
                // flop count analysis. It also means that the post-ordering
                // is not an exact elimination tree post-ordering.

                if DEBUG_LEVEL >= 1 {
                    print!("  MASS i {} => parent e {}\n", i, me);
                }
                Pe[i] = flip(me);
                let nvi = -Nv[i];
                degme -= nvi;
                nvpiv += nvi;
                nel += nvi;
                Nv[i] = 0;
                Elen[i] = EMPTY;
            } else {
                // Update the upper-bound degree of i.

                // The following degree does not yet include the size
                // of the current element, which is added later:

                degree[i] = min(degree[i], deg);

                // Add me to the list for i.

                // Move first supervariable to end of list.
                Iw[pn] = Iw[p3];
                // Move first element to end of element part of list.
                Iw[p3] = Iw[p1];
                // Add new element, me, to front of list.
                Iw[p1] = me;
                // Store the new length of the list in Len[i].
                Len[i] = pn - p1 + 1;

                // Place in hash bucket. Save hash key of i in Last[i].

                // FIXME: this can fail if hash is negative, because the ANSI C
                // standard does not define a % b when a and/or b are negative.
                // That's why hash is defined as an unsigned int, to avoid this
                // problem.
                hash = hash % n as u32;
                if DEBUG_LEVEL != 0 {
                    debug_assert!(hash >= 0 && hash < n as u32);
                }

                // If the Hhead array is not used:
                let j = Head[hash];
                if j <= EMPTY {
                    // Degree list is empty, hash head is flip(j).
                    Next[i] = flip(j);
                    Head[hash] = flip(i);
                } else {
                    // Degree list is not empty, use Last [Head[hash]] as hash head.
                    Next[i] = Last[j];
                    Last[j] = i;
                }

                // If a separate Hhead array is used:
                // 	Next [i] = Hhead[hash]
                // 	Hhead [hash] = i

                Last[i] = hash as i32;
            }
        }

        degree[me] = degme;

        // Clear the counter array, W [...], by incrementing wflg.

        // Make sure that wflg+n does not cause integer overflow.
        lemax = max(lemax, degme);
        wflg += lemax;
        wflg = clear_flag(wflg, wbig, W, n);
        // at this point, W[0..n-1] < wflg holds

        /* Supervariable Detection */

        if DEBUG_LEVEL >= 1 {
            print!("Detecting supervariables:\n");
        }
        for pme in pme1..=pme2 {
            let i = Iw[pme];
            if DEBUG_LEVEL != 0 {
                debug_assert!(i >= 0 && i < n);
                if DEBUG_LEVEL >= 2 {
                    print!("Consider i {} nv {}\n", i, Nv[i]);
                }
            }

            if Nv[i] < 0 {
                // i is a principal variable in Lme.

                // Examine all hash buckets with 2 or more variables. We do
                // this by examing all unique hash keys for supervariables in
                // the pattern Lme of the current element, me.

                // Let i = head of hash bucket, and empty the hash bucket.
                if DEBUG_LEVEL != 0 {
                    debug_assert!(Last[i] >= 0 && Last[i] < n);
                }
                hash = Last[i] as u32;

                // If Hhead array is not used:
                let j = Head[hash];
                if j == EMPTY {
                    // hash bucket and degree list are both empty.
                    i = EMPTY;
                } else if j < EMPTY {
                    // Degree list is empty.
                    i = flip(j);
                    Head[hash] = EMPTY;
                } else {
                    // Degree list is not empty, restore Last[j] of head j.
                    i = Last[j];
                    Last[j] = EMPTY;
                }

                // If separate Hhead array is used:
                // i = Hhead[hash]
                // Hhead[hash] = empty

                if DEBUG_LEVEL != 0 {
                    debug_assert!(i >= EMPTY && i < n);
                    if DEBUG_LEVEL >= 2 {
                        print!("----i {} hash {}\n", i, hash);
                    }
                }

                while i != EMPTY && Next[i] != EMPTY {
                    // This bucket has one or more variables following i.
                    // scan all of them to see if i can absorb any entries
                    // that follow i in hash bucket. Scatter i into w.

                    let ln = Len[i];
                    let eln = Elen[i];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(ln >= 0 && eln >= 0);
                        debug_assert!(Pe[i] >= 0 && Pe[i] < iwlen);
                    }
                    // Do not flag the first element in the list(me).
                    for p in Pe[i] + 1..=Pe[i] + ln - 1 {
                        if DEBUG_LEVEL != 0 {
                            debug_assert!(Iw[p] >= 0 && Iw[p] < n);
                        }
                        W[Iw[p]] = wflg;
                    }

                    // Scan every other entry j following i in bucket.

                    let jlast = i;
                    j = Next[i];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(j >= EMPTY && j < n);
                    }

                    while j != EMPTY {
                        // Check if j and i have identical nonzero pattern.

                        if DEBUG_LEVEL >= 3 {
                            print!("compare i {} and j {}", i, j);
                        }

                        if DEBUG_LEVEL != 0 {
                            // Check if i and j have the same Len and Elen.
                            debug_assert!(Len[j] >= 0 && Elen[j] >= 0);
                            debug_assert!(Pe[j] >= 0 && Pe[j] < iwlen);
                        }
                        let mut ok = (Len[j] == ln) && (Elen[j] == eln);
                        // Skip the first element in the list(me).
                        // TODO: for p := Pe[j] + 1; ok && p <= Pe[j]+ln-1; p++ {
                        for p in Pe[j] + 1..=Pe[j] + ln - 1 {
                            if DEBUG_LEVEL != 0 {
                                debug_assert!(Iw[p] >= 0 && Iw[p] < n);
                            }
                            if W[Iw[p]] != wflg {
                                ok = false;
                                break;
                            }
                        }
                        if ok {
                            // Found it  j can be absorbed into i.

                            if DEBUG_LEVEL >= 1 {
                                print!("found it! j {} => i {}\n", j, i);
                            }
                            Pe[j] = flip(i);
                            // Both Nv[i] and Nv[j] are negated since they
                            // are in Lme, and the absolute values of each
                            // are the number of variables in i and j:
                            Nv[i] += Nv[j];
                            Nv[j] = 0;
                            Elen[j] = EMPTY;
                            // Delete j from hash bucket.
                            if DEBUG_LEVEL != 0 {
                                debug_assert!(j != Next[j]);
                            }
                            j = Next[j];
                            Next[jlast] = j;
                        } else {
                            // j cannot be absorbed into i.
                            jlast = j;
                            if DEBUG_LEVEL != 0 {
                                debug_assert!(j != Next[j]);
                            }
                            j = Next[j];
                        }
                        if DEBUG_LEVEL != 0 {
                            debug_assert!(j >= EMPTY && j < n);
                        }
                    }

                    // No more variables can be absorbed into
                    // go to next i in bucket and clear flag array.

                    wflg += 1;
                    i = Next[i];
                    if DEBUG_LEVEL != 0 {
                        debug_assert!(i >= EMPTY && i < n);
                    }
                }
            }
        }
        if DEBUG_LEVEL >= 2 {
            println!("detect done");
        }

        // Restore degree lists and remove nonprincipal supervariables from element.

        let mut p = pme1;
        let nleft = n - nel;
        for pme in pme1..=pme2 {
            let i = Iw[pme];
            if DEBUG_LEVEL != 0 {
                debug_assert!(i >= 0 && i < n);
            }
            let nvi = -Nv[i];
            if DEBUG_LEVEL >= 3 {
                print!("Restore i {} {}\n", i, nvi);
            }
            if nvi > 0 {
                // i is a principal variable in Lme.
                // Restore Nv[i] to signify that i is principal.
                Nv[i] = nvi;

                // Compute the external degree (add size of current element).

                deg = degree[i] + degme - nvi;
                deg = min(deg, nleft - nvi);
                if DEBUG_LEVEL != 0 {
                    debug_assert!(implies(aggressive != 0, deg > 0) && deg >= 0 && deg < n);
                }

                // Place the supervariable at the head of the degree list.

                inext = Head[deg];
                if DEBUG_LEVEL != 0 {
                    debug_assert!(inext >= EMPTY && inext < n);
                }
                if inext != EMPTY {
                    Last[inext] = i;
                }
                Next[i] = inext;
                Last[i] = EMPTY;
                Head[deg] = i;

                // Save the new degree, and find the minimum degree.

                mindeg = min(mindeg, deg);
                degree[i] = deg;

                // Place the supervariable in the element pattern.

                Iw[p] = i;
                p += 1;
            }
        }
        if DEBUG_LEVEL >= 2 {
            println!("restore done");
        }

        // Finalize the new element.

        if DEBUG_LEVEL >= 2 {
            print!("ME = {} DONE\n", me);
        }
        Nv[me] = nvpiv;
        // Save the length of the list for the new element me.
        Len[me] = p - pme1;
        if Len[me] == 0 {
            // There is nothing left of the current pivot element.
            // It is a root of the assembly tree.
            Pe[me] = EMPTY;
            W[me] = 0;
        }
        if elenme != 0 {
            // Element was not constructed in place: deallocate part of
            // it since newly nonprincipal variables may have been removed.
            pfree = p;
        }

        // The new element has nvpiv pivots and the size of the contribution
        // block for a multifrontal method is degme-by-degme, not including
        // the "dense" rows/columns. If the "dense" rows/columns are included,
        // the frontal matrix is no larger than
        // (degme+ndense)-by-(degme+ndense).

        {
            let f = nvpiv as f64;
            let r = (degme + ndense) as f64;
            dmax = max(dmax, f + r);

            // Number of nonzeros in L (excluding the diagonal).
            let lnzme = f * r + (f - 1) * f / 2;
            lnz += lnzme;

            // Number of divide operations for LDL' and for LU.
            ndiv += lnzme;

            // Number of multiply-subtract pairs for LU.
            let s = f * r * r + r * (f - 1) * f + (f - 1) * f * (2 * f - 1) / 6;
            nms_lu += s;

            // Number of multiply-subtract pairs for LDL'.
            nms_ldl += (s + lnzme) / 2;
        }

        if DEBUG_LEVEL >= 2 {
            print!("finalize done nel {} n {}\n   ::::\n", nel, n);
            if DEBUG_LEVEL >= 3 {
                for pme in Pe[me]..=Pe[me] + Len[me] - 1 {
                    print!(" {}", Iw[pme]);
                }
                println!();
            }
        }
    }

    // Done selecting pivots.

    {
        // Count the work to factorize the ndense-by-ndense submatrix.
        let f = ndense as f64;
        dmax = max(dmax, ndense as f64);

        // Number of nonzeros in L (excluding the diagonal).
        let lnzme = (f - 1) * f / 2;
        lnz += lnzme;

        // Number of divide operations for LDL' and for LU.
        ndiv += lnzme;

        // Number of multiply-subtract pairs for LU.
        let s = (f - 1) * f * (2 * f - 1) / 6;
        nms_lu += s;

        // Number of multiply-subtract pairs for LDL'.
        nms_ldl += (s + lnzme) / 2;

        // Number of nz's in L (excl. diagonal).
        info.lnz = lnz;

        // Number of divide ops for LU and LDL'.
        info.n_div = ndiv;

        // Number of multiply-subtract pairs for LDL'.
        info.n_mult_subs_ldl = nms_ldl;

        // Number of multiply-subtract pairs for LU.
        info.n_mult_subs_lu = nms_lu;

        // Number of "dense" rows/columns.
        info.n_dense = ndense;

        // Largest front is dmax-by-dmax.
        info.d_max = dmax;

        // Number of garbage collections in AMD.
        info.n_cmp_a = ncmpa;

        // Successful ordering.
        info.status = Status::OK;
    }

    /* Post-ordering */

    // Variables at this point:
    //
    // Pe: holds the elimination tree. The parent of j is flip(Pe[j]),
    // or EMPTY if j is a root. The tree holds both elements and
    // non-principal (unordered) variables absorbed into them.
    // Dense variables are non-principal and unordered.
    //
    // Elen: holds the size of each element, including the diagonal part.
    // flip(Elen[e]) > 0 if e is an element. For unordered
    // variables i, Elen[i] is EMPTY.
    //
    // Nv: Nv[e] > 0 is the number of pivots represented by the element e.
    // For unordered variables i, Nv[i] is zero.
    //
    // Contents no longer needed:
    // W, Iw, Len, Degree, Head, Next, Last.
    //
    // The matrix itself has been destroyed.
    //
    // n: the size of the matrix.
    // No other scalars needed (pfree, iwlen, etc.)

    // Restore Pe.
    for i in 0..n {
        Pe[i] = flip(Pe[i]);
    }

    // Restore Elen, for output information, and for postordering.
    for i in 0..n {
        Elen[i] = flip(Elen[i]);
    }

    // Now the parent of j is Pe[j], or EMPTY if j is a root. Elen[e] > 0
    // is the size of element e. Elen [i] is EMPTY for unordered variable i.

    if DEBUG_LEVEL >= 2 {
        println!("\nTree:");
        for i in 0..n {
            print!(" {} parent: {}   \n", i, Pe[i]);
            debug_assert!(Pe[i] >= EMPTY && Pe[i] < n);
            if Nv[i] > 0 {
                // This is an element.
                let e = i;
                print!(" element, size is {}", Elen[i]);
                debug_assert!(Elen[e] > 0);
            }
            println!();
        }
        println!("\nelements:");
        for e in 0..n {
            if Nv[e] > 0 {
                print!("Element e = {} size {} nv {} \n", e, Elen[e], Nv[e]);
            }
        }
        if DEBUG_LEVEL >= 3 {
            println!("\nvariables:");
            for i in 0..n {
                let cnt: i32;
                if Nv[i] == 0 {
                    print!("i unordered: {}\n", i);
                    let j = Pe[i];
                    cnt = 0;
                    print!("  j: {}\n", j);
                    if j == EMPTY {
                        println!(" i is a dense variable");
                    } else {
                        debug_assert!(j >= 0 && j < n);
                        while Nv[j] == 0 {
                            print!(" j : {}\n", j);
                            j = Pe[j];
                            print!(" j:: {}\n", j);
                            cnt += 1;
                            if cnt > n {
                                break;
                            }
                        }
                        let e = j;
                        print!(" got to e: {}\n", e);
                    }
                }
            }
        }
    }

    // Compress the paths of the variables.

    for i in 0..n {
        if Nv[i] == 0 {
            // i is an un-ordered row. Traverse the tree from i until
            // reaching an element, e. The element, e, was the principal
            // supervariable of i and all nodes in the path from i to when e
            // was selected as pivot.

            if DEBUG_LEVEL >= 1 {
                print!("Path compression, i unordered: {}\n", i);
            }
            let j = Pe[i];
            if DEBUG_LEVEL != 0 {
                debug_assert!(j >= EMPTY && j < n);
                if DEBUG_LEVEL >= 3 {
                    print!(" j: {}\n", j);
                }
            }
            if j == EMPTY {
                // Skip a dense variable. It has no parent.
                if DEBUG_LEVEL >= 3 {
                    print!("      i is a dense variable\n");
                }
                continue;
            }

            // while (j is a variable)
            while Nv[j] == 0 {
                if DEBUG_LEVEL >= 3 {
                    print!("  j : {}\n", j);
                }
                j = Pe[j];
                if DEBUG_LEVEL != 0 {
                    if DEBUG_LEVEL >= 3 {
                        print!("  j:: {}\n", j);
                    }
                    debug_assert!(j >= 0 && j < n);
                }
            }
            // Got to an element e.
            let e = j;
            if DEBUG_LEVEL >= 3 {
                print!("got to e: {}\n", e);
            }

            // Traverse the path again from i to e, and compress the path
            // (all nodes point to e). Path compression allows this code to
            // compute in O(n) time.

            j = i;
            // while (j is a variable)
            while Nv[j] == 0 {
                let jnext = Pe[j];
                if DEBUG_LEVEL >= 3 {
                    print!("j {} jnext {}\n", j, jnext);
                }
                Pe[j] = e;
                j = jnext;
                if DEBUG_LEVEL != 0 {
                    debug_assert!(j >= 0 && j < n);
                }
            }
        }
    }

    // postorder the assembly tree

    postorder(n, Pe, Nv, Elen, W, Head, Next, Last);

    // Compute output permutation and inverse permutation.

    // W[e] = k means that element e is the kth element in the new
    // order. e is in the range 0 to n-1, and k is in the range 0 to
    // the number of elements. Use Head for inverse order.

    for k in 0..n {
        Head[k] = EMPTY;
        Next[k] = EMPTY;
    }
    for e in 0..n {
        let k = W[e];
        if DEBUG_LEVEL != 0 {
            debug_assert!((k == EMPTY) == (Nv[e] == 0));
        }
        if k != EMPTY {
            if DEBUG_LEVEL != 0 {
                debug_assert!(k >= 0 && k < n);
            }
            Head[k] = e;
        }
    }

    // Construct output inverse permutation in Next, and permutation in Last.
    nel = 0;
    for k in 0..n {
        let e = Head[k];
        if e == EMPTY {
            break;
        }
        if DEBUG_LEVEL != 0 {
            debug_assert!(e >= 0 && e < n && Nv[e] > 0);
        }
        Next[e] = nel;
        nel += Nv[e];
    }
    if DEBUG_LEVEL != 0 {
        debug_assert!(nel == n - ndense);
    }

    // Order non-principal variables (dense, & those merged into supervar's).
    for i in 0..n {
        if Nv[i] == 0 {
            let e = Pe[i];
            if DEBUG_LEVEL != 0 {
                debug_assert!(e >= EMPTY && e < n);
            }
            if e != EMPTY {
                // This is an unordered variable that was merged
                // into element e via supernode detection or mass
                // elimination of i when e became the pivot element.
                // Place i in order just before e.
                if DEBUG_LEVEL != 0 {
                    debug_assert!(Next[i] == EMPTY && Nv[e] > 0);
                }
                Next[i] = Next[e];
                Next[e] += 1;
            } else {
                // This is a dense unordered variable, with no parent.
                // Place it last in the output order.
                Next[i] = nel;
                nel += 1;
            }
        }
    }
    if DEBUG_LEVEL != 0 {
        debug_assert!(nel == n);
    }

    if DEBUG_LEVEL >= 2 {
        print!("\n\nPerm:\n");
    }
    for i in 0..n {
        let k = Next[i];
        if DEBUG_LEVEL != 0 {
            debug_assert!(k >= 0 && k < n);
        }
        Last[k] = i;
        if DEBUG_LEVEL >= 2 {
            print!("   perm [{}] = {}\n", k, i);
        }
    }
}
