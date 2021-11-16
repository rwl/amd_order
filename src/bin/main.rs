extern crate amd_order;

use amd_order::*;

fn main() {
    let n: i32 = 5;
    let a_p: Vec<i32> = vec![0, 2, 6, 10, 12, 14];
    let a_i: Vec<i32> = vec![
        0, 1, // 1st row
        0, 1, 2, 4, // 2nd row
        1, 2, 3, 4, // 3rd row
        2, 3, // 4th row
        1, 4, // 5th row
    ];
    let control = amd_order::amd::default_control_settings();

    let (permutation, info) = order(n, &a_p, &a_i, control);

    println!("P = {:?}", permutation);
    println!("Info = {:?}", info);
}
