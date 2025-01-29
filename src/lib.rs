//#![no_std]
#![allow(improper_ctypes_definitions)]

//! Machine-factor is a relatively fast library for factoring integers up to 2^128. Most integers can be factored in less than  
//! 60 seconds (on i5-5300U). Machine-factor can be used in const contexts, however this is will often exceed the allowed time. 

mod mfactor;
mod primes;
mod wfactor;

pub use mfactor::{get_factor,factorize,Factorization};

pub use wfactor::{get_factor_128,factorize_128,Factorization};
#[cfg(feature="internal")]
pub use {primes::PRIMES_128,wfactor::gcd_128,mfactor::{drbg,gcd}};
/*
#[panic_handler]
fn panic(_info: &core::panic::PanicInfo) -> ! {
    loop {}
}
*/
