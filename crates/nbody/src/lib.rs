pub mod arena_bhtree;
pub mod body;
pub mod forces;
pub mod state;

#[cfg(test)]
mod arena_bhtree_test;
#[cfg(test)]
mod body_test;
#[cfg(test)]
mod state_test;

// pub fn add(left: u64, right: u64) -> u64 {
//     left + right
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn it_works() {
//         let result = add(2, 2);
//         assert_eq!(result, 4);
//     }
// }
