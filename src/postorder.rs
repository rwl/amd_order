use crate::amd::*;
use crate::internal::*;
use crate::post_tree::post_tree;

pub fn postorder(
    nn: i32,
    parent: &[i32],
    Nv: &[i32],
    Fsize: &[i32],
    order: &[i32],
    child: &[i32],
    sibling: &[i32],
    stack: &[i32],
) {
    let mut nchild = 0;

    for j in 0..nn {
        child[j] = EMPTY;
        sibling[j] = EMPTY;
    }

    // Place the children in link lists - bigger elements tend to be last.
    // for j := nn - 1; j >= 0; j-- {
    for j in (nn - 1..=0).rev() {
        if Nv[j] > 0 {
            // This is an element.
            let p = parent[j];
            if p != EMPTY {
                // Place the element in link list of the children its parent
                // bigger elements will tend to be at the end of the list.
                sibling[j] = child[p];
                child[p] = j;
            }
        }
    }

    if DEBUG_LEVEL != 0 {
        // var nels, ff int
        print!("\n\n================================ AMD_postorder:\n");
        let mut nels = 0;
        for j in 0..nn {
            if Nv[j] > 0 {
                print!(
                    "{} :  nels {} npiv {} size {} parent {} maxfr {}\n",
                    j, nels, Nv[j], Fsize[j], parent[j], Fsize[j]
                );
                // This is an element. Dump the link list of children.
                nchild = 0;
                let mut d1 = "    Children: ";
                let ff = child[j];
                while ff != EMPTY {
                    d1 += format!("{} ", ff);
                    assert!(parent[ff] == j);
                    nchild += 1;
                    assert!(nchild < nn);

                    ff = sibling[ff];
                }
                println!("{}", d1);
                let p = parent[j];
                if p != EMPTY {
                    assert!(Nv[p] > 0);
                }
                nels += 1;
            }
        }
        println!(
            "\n\nGo through the children of each node, and put
the biggest child last in each list:"
        );
    }

    // Place the largest child last in the list of children for each node.
    for i in 0..nn {
        if Nv[i] > 0 && child[i] != EMPTY {
            if DEBUG_LEVEL != 0 {
                print!("Before partial sort, element {}\n", i);
                nchild = 0;

                let mut f = child[i];
                while f != EMPTY {
                    assert!(f >= 0 && f < nn);
                    print!("      f: {}  size: {}\n", f, Fsize[f]);
                    nchild += 1;
                    assert!(nchild <= nn);

                    f = sibling[f];
                }
            }

            // Find the biggest element in the child list.
            let mut fprev = EMPTY;
            let mut maxfrsize = EMPTY;
            let mut bigfprev = EMPTY;
            let mut bigf = EMPTY;

            let mut f = child[i];
            while f != EMPTY {
                if DEBUG_LEVEL != 0 {
                    debug_assert!(f >= 0 && f < nn);
                }
                let frsize = Fsize[f];
                if frsize >= maxfrsize {
                    // This is the biggest seen so far.
                    maxfrsize = frsize;
                    bigfprev = fprev;
                    bigf = f;
                }
                fprev = f;

                f = sibling[f];
            }
            if DEBUG_LEVEL != 0 {
                debug_assert!(bigf != EMPTY);
            }

            let fnext = sibling[bigf];

            if DEBUG_LEVEL >= 1 {
                print!(
                    "bigf {} maxfrsize {} bigfprev {} fnext {} fprev {}\n",
                    bigf, maxfrsize, bigfprev, fnext, fprev
                );
            }

            if fnext != EMPTY {
                // If fnext is EMPTY then bigf is already at the end of list.

                if bigfprev == EMPTY {
                    // Delete bigf from the element of the list.
                    child[i] = fnext;
                } else {
                    // Delete bigf from the middle of the list.
                    sibling[bigfprev] = fnext;
                }

                // Put bigf at the end of the list.
                sibling[bigf] = EMPTY;
                if DEBUG_LEVEL != 0 {
                    assert!(child[i] != EMPTY);
                    assert!(fprev != bigf);
                    assert!(fprev != EMPTY);
                }
                sibling[fprev] = bigf;
            }

            if DEBUG_LEVEL != 0 {
                print!("After partial sort, element {}\n", i);
                let mut f = child[i];
                while f != EMPTY {
                    assert!(f >= 0 && f < nn);
                    print!("        {}  {}\n", f, Fsize[f]);
                    assert!(Nv[f] > 0);
                    nchild -= 1;

                    f = sibling[f];
                }
                assert!(nchild == 0);
            }
        }
    }

    // Postorder the assembly tree.
    for i in 0..nn {
        order[i] = EMPTY;
    }

    let mut k = 0;
    for i in 0..nn {
        if parent[i] == EMPTY && Nv[i] > 0 {
            if DEBUG_LEVEL != 0 {
                print!("Root of assembly tree {}\n", i);
            }
            k = post_tree(i, k, child, sibling, order, stack, nn);
        }
    }
}
