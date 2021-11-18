extern crate amd_order;

use amd_order::amd::default_control_settings;
use amd_order::order;

fn main() {
    let n: i32 = 5;
    let a_p: Vec<i32> = vec![0, 2, 6, 10, 12, 14];
    let a_i: Vec<i32> = vec![
        0, 1, // 1st column
        0, 1, 2, 4, // 2nd column
        1, 2, 3, 4, // 3rd column
        2, 3, // 4th column
        1, 4, // 5th column
    ];
    let control = default_control_settings();

    let (permutation, info) = order(n, &a_p, &a_i, control);

    println!("P = {:?}", permutation);
    // Output:
    //   [0 3 2 4 1]

    amd_order::info::info(&info);
}
