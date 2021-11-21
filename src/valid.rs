use crate::amd::*;
use crate::internal::*;

pub fn valid(n_row: usize, n_col: usize, a_p: &[usize], a_i: &[usize]) -> Status {
    if a_p[0] != 0 {
        // Column pointers must start at `Ap[0] = 0`, and `Ap[n]` must be `>= 0`.
        debug1_println!("column 0 pointer bad or nz < 0");
        return Status::Invalid;
    }

    let mut status = Status::OK;

    for j in 0..n_col {
        let p1 = a_p[j];
        let p2 = a_p[j + 1];
        debug2_print!("\nColumn: {} p1: {} p2: {}\n", j, p1, p2);

        if p1 > p2 {
            // Column pointers must be ascending.
            debug1_print!("column {} pointer bad\n", j);
            return Status::Invalid;
        }

        let mut ilast: isize = EMPTY;

        for p in p1..p2 {
            let i = a_i[p] as isize;
            debug3_print!("row: {}\n", i);

            if i < 0 || i >= n_row as isize {
                // Row index out of range.
                debug1_print!("index out of range, col {} row {}\n", j, i);
                return Status::Invalid;
            }
            if i <= ilast {
                // Row index unsorted, or duplicate entry present.
                debug1_print!("index unsorted/dupl col {} row {}\n", j, i);
                status = Status::OkButJumbled;
            }
            ilast = i
        }
    }

    status
}
