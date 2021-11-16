// Flip is a "negation about -1", and is used to mark an integer i that is
// normally non-negative. flip(EMPTY) is empty. Flip of a number > EMPTY
// is negative, and flip of a number < EMTPY is positive. flip(flip(i)) = i
// for all integers i. unflip(i) is >= EMPTY.
pub const EMPTY: i32 = -1;

pub fn flip(i: i32) -> i32 {
    -i - 2
}
