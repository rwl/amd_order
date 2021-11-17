use crate::amd::*;
use crate::internal::EMPTY;

pub fn valid(n_row: i32, n_col: i32, a_p: &[i32], a_i: &[i32]) -> Status {
    if n_row < 0 || n_col < 0 {
        return Status::Invalid;
    }

    let nz = a_p[n_col as usize];
    if a_p[0] != 0 || nz < 0 {
        // Column pointers must start at `Ap[0] = 0`, and `Ap[n]` must be `>= 0`.
        if DEBUG_LEVEL > 0 {
            println!("column 0 pointer bad or nz < 0");
        }
        return Status::Invalid;
    }

    let mut status = Status::OK;

    for j in 0..n_col {
        let p1 = a_p[j as usize];
        let p2 = a_p[j as usize + 1];
        if DEBUG_LEVEL >= 2 {
            print!("\nColumn: {} p1: {} p2: {}\n", j, p1, p2)
        }

        if p1 > p2 {
            // Column pointers must be ascending.
            if DEBUG_LEVEL >= 0 {
                print!("column {} pointer bad\n", j)
            }
            return Status::Invalid;
        }

        let mut ilast: i32 = EMPTY;

        for p in p1..p2 {
            let i = a_i[p as usize];

            if DEBUG_LEVEL >= 3 {
                print!("row: {}\n", i)
            }
            if i < 0 || i >= n_row {
                // Row index out of range.
                if DEBUG_LEVEL >= 0 {
                    print!("index out of range, col {} row {}\n", j, i)
                }
                return Status::Invalid;
            }
            if i <= ilast {
                // Row index unsorted, or duplicate entry present.
                if DEBUG_LEVEL >= 1 {
                    print!("index unsorted/dupl col {} row {}\n", j, i)
                }
                status = Status::OkButJumbled;
            }
            ilast = i
        }
    }

    return status;
}
