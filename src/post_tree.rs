use crate::amd::*;
use crate::internal::*;

pub fn post_tree(
    root: i32,
    mut k: i32,
    child: &[i32],
    sibling: &[i32],
    order: &[i32],
    stack: &[i32],
    nn: i32,
) -> i32 {
    /*if false {
        // Recursive version (stack[] is not used):
        // this is simple, but can can cause stack overflow if nn is large
        i = root;
        f = child[i];
        while f != EMPTY {
            k = post_tree(f, k, child, sibling, order, stack, nn);
            f = sibling[f];
        }
        order[i] = k;
        k += 1;
        return k;
    }*/

    // Non-recursive version, using an explicit stack.

    // Push root on the stack.
    let mut head = 0;
    stack[0] = root;

    while head >= 0 {
        // Get head of stack.
        if DEBUG_LEVEL != 0 {
            debug_assert!(head < nn);
        }
        let i = stack[head];
        if DEBUG_LEVEL != 0 {
            print!("head of stack {} \n", i);
            debug_assert!(i >= 0 && i < nn);
        }

        if child[i] != EMPTY {
            // The children of i are not yet ordered
            // push each child onto the stack in reverse order
            // so that small ones at the head of the list get popped first
            // and the biggest one at the end of the list gets popped last.
            let mut f = child[i];
            while f != EMPTY {
                head += 1;
                if DEBUG_LEVEL != 0 {
                    debug_assert!(head < nn);
                    debug_assert!(f >= 0 && f < nn);
                }
                f = sibling[f];
            }
            let mut h = head;
            if DEBUG_LEVEL != 0 {
                debug_assert!(head < nn);
            }
            let mut f = child[i];
            while f != EMPTY {
                if DEBUG_LEVEL != 0 {
                    debug_assert!(h > 0);
                }
                stack[h] = f;
                h -= 1;
                if DEBUG_LEVEL != 0 {
                    println!("push {} on stack", f);
                    debug_assert!(f >= 0 && f < nn);
                }
                f = sibling[f];
            }
            if DEBUG_LEVEL != 0 {
                debug_assert!(stack[h] == i);
            }

            // Delete child list so that i gets ordered next time we see it.
            child[i] = EMPTY;
        } else {
            // The children of i (if there were any) are already ordered
            // remove i from the stack and order it. Front i is kth front.
            head -= 1;
            if DEBUG_LEVEL >= 1 {
                print!("pop {} order {}\n", i, k);
            }
            order[i] = k;
            k += 1;
            if DEBUG_LEVEL != 0 {
                debug_assert!(k <= nn);
            }
        }

        if DEBUG_LEVEL != 0 {
            let mut d1 = "\nStack:";
            // for h := head; h >= 0; h-- {
            for h in (head..=0).rev() {
                let j = stack[h];
                d1 += format!(" {}", j);
                assert!(j >= 0 && j < nn);
            }
            print!("{}\n\n", d1);
            assert!(head < nn);
        }
    }

    return k;
}
