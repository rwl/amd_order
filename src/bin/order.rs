extern crate amd_order;

use amd_order::amd::{default_control_settings, Status};
use amd_order::order;

// A simple test that illustrates the use of the interface to AMD.
fn main() {
    // The symmetric can_24 Harwell/Boeing matrix, including upper and lower
    // triangular parts, and the diagonal entries.  Note that this matrix is
    // 0-based, with row and column indices in the range 0 to n-1.
    let n: usize = 24;

    let a_p = vec![
        0, 9, 15, 21, 27, 33, 39, 48, 57, 61, 70, 76, 82, 88, 94, 100, 106, 110, 119, 128, 137,
        143, 152, 156, 160,
    ];

    let a_i = vec![
        0, 5, 6, 12, 13, 17, 18, 19, 21, // column 0
        1, 8, 9, 13, 14, 17, // column 1
        2, 6, 11, 20, 21, 22, // column 2
        3, 7, 10, 15, 18, 19, // column 3
        4, 7, 9, 14, 15, 16, // column 4
        0, 5, 6, 12, 13, 17, // column 5
        0, 2, 5, 6, 11, 12, 19, 21, 23, // column 6
        3, 4, 7, 9, 14, 15, 16, 17, 18, // column 7
        1, 8, 9, 14, // column 8
        1, 4, 7, 8, 9, 13, 14, 17, 18, // column 9
        3, 10, 18, 19, 20, 21, // column 10
        2, 6, 11, 12, 21, 23, // column 11
        0, 5, 6, 11, 12, 23, // column 12
        0, 1, 5, 9, 13, 17, // column 13
        1, 4, 7, 8, 9, 14, // column 14
        3, 4, 7, 15, 16, 18, // column 15
        4, 7, 15, 16, // column 16
        0, 1, 5, 7, 9, 13, 17, 18, 19, // column 17
        0, 3, 7, 9, 10, 15, 17, 18, 19, // column 18
        0, 3, 6, 10, 17, 18, 19, 20, 21, // column 19
        2, 10, 19, 20, 21, 22, // column 20
        0, 2, 6, 10, 11, 19, 20, 21, 22, // column 21
        2, 20, 21, 22, // column 22
        6, 11, 12, 23, // column 23
    ];

    let mut p_inv = vec![0; 24];
    let control = default_control_settings();
    let mut a = [[""; 24]; 24];

    print!("AMD demo, with the 24-by-24 Harwell/Boeing matrix, can_24:\n");

    amd_order::control::control(&control);

    // Print the input matrix.
    let nz = a_p[n];
    print!(
        "\nInput matrix:  {}-by-{}, with {} entries.
    Note that for a symmetric matrix such as this one, only the
    strictly lower or upper triangular parts would need to be
    passed to AMD, since AMD computes the ordering of A+A'. The
    diagonal entries are also not needed, since AMD ignores them.\n",
        n, n, nz
    );

    for j in 0..n {
        print!(
            "\nColumn: {}, number of entries: {}, with row indices in
 Ai [{} ... {}]:
    row indices:",
            j,
            a_p[j + 1] - a_p[j],
            a_p[j],
            a_p[j + 1] - 1
        );
        for pj in a_p[j]..a_p[j + 1] {
            let i = a_i[pj as usize];
            print!(" {}", i);
        }
        print!("\n");
    }

    // Print a character plot of the input matrix. This is only reasonable
    // because the matrix is small.
    print!("\nPlot of input matrix pattern:\n");
    for j in 0..n {
        for i in 0..n {
            a[i][j] = "."
        }
        for pj in a_p[j]..a_p[j + 1] {
            let i = a_i[pj as usize] as usize;
            a[i][j] = "X";
        }
    }
    print!("    ");
    for j in 0..n {
        print!(" {}", j % 10);
    }
    print!("\n");
    for i in 0..n {
        print!("{}: ", i);
        for j in 0..n {
            print!(" {}", a[i][j]);
        }
        print!("\n");
    }

    // Order the matrix.
    let (p, info) = order(n as i32, &a_p, &a_i, &control).unwrap();
    print!(
        "return value from amd_order: {:?} (should be {:?})\n",
        info.status,
        Status::OK
    );

    // Print the statistics.
    amd_order::info::info(&info);

    if info.status != Status::OK {
        print!("AMD failed\n");
        return;
    }

    // Print the permutation vector, P, and compute the inverse permutation.
    print!("Permutation vector:\n");
    for k in 0..n {
        // Row/column j is the kth row/column in the permuted matrix.
        let j = p[k];
        p_inv[j as usize] = k;
        print!(" {}", j);
    }
    print!("\n\n");

    print!("Inverse permutation vector:\n");
    for j in 0..n {
        let k = p_inv[j];
        print!(" {}", k);
    }
    print!("\n\n");

    // Print a character plot of the permuted matrix.
    print!("\nPlot of permuted matrix pattern:\n");
    for jnew in 0..n {
        let j = p[jnew] as usize;
        for inew in 0..n {
            a[inew][jnew] = ".";
        }
        for pj in a_p[j]..a_p[j + 1] {
            let i = a_i[pj as usize];
            let inew = p_inv[i as usize] as usize;
            a[inew][jnew] = "X";
        }
    }
    print!("    ");
    for j in 0..n {
        print!(" {}", j % 10);
    }
    print!("\n");
    for i in 0..n {
        print!("{}: ", i);
        for j in 0..n {
            print!(" {}", a[i][j]);
        }
        print!("\n");
    }
}
