use crate::mfactor::drbg;
use crate::primes::PRIMES_130;

/// GCD(A,B) for 128-bit integers
#[no_mangle]
pub const extern "C" fn gcd_128(mut a: u128, mut b: u128) -> u128{
        if b == 0 {
            return a;
        } else if a == 0 {
            return b;
        }

        let self_two_factor = a.trailing_zeros();
        let other_two_factor = b.trailing_zeros();
        let mut min_two_factor = self_two_factor;
        
        if other_two_factor < self_two_factor{
           min_two_factor=other_two_factor;
        }
        
        a >>= self_two_factor;
        b >>= other_two_factor;
        loop {
            if b > a {
            let interim = b;
                b = a;
                a = interim;

            }
            a -= b;

            if a == 0 {
                return b << min_two_factor;
            }
            a >>= a.trailing_zeros();
        }
}

const fn poly_eval_128(x: u128, subtrahend: u128, n: u128, npi: u128) -> u128{
     machine_prime::mont_sub_128(machine_prime::mont_sqr_128(x,n,npi),subtrahend,n)
}


const fn pollard_brent_128(base: u128,inv:u128,subtrahend: u128, n: u128) -> Option<u128>{
    let m = 512;
    let mut r = 1;
    let mut q = 1;
    let mut g = 1;
    let mut ys = 1;
    let mut y = base;
    let mut x = y;
    let mut cycle = 0;
    
    while cycle < 33{
    
     cycle+=1;    
      x = y;
      
      let mut yloop = 0;
      
      while yloop < r{
      yloop+=1; 
        y = poly_eval_128(y,subtrahend,n,inv);      
      }
      
      let mut k = 0;
      
      loop{
      
      let mut i = 0;
      
       while i < m*cycle{
           if i >= r-k{
             break;
           }
         
         y=poly_eval_128(y,subtrahend,n,inv);
         q=machine_prime::mont_prod_128(q,x.abs_diff(y),n,inv);
         i+=1;
         } // end loop

         ys=y;
         g = gcd_128(q,n);
         k+=m;
         if k >= r || g !=1{
            break;
         }
      }
      
      r<<=1;
      if g != 1{
         break;
      }
      
    }
    
    if g ==n{
       while g==1{
         ys=poly_eval_128(ys,subtrahend,n,inv);
         g=gcd_128(x.abs_diff(ys),n);
      
       }
    }
    if g !=1 && g !=n && machine_prime::is_prime_wc_128(g){
       return Some(g);
    }
    None
}

/// Returns some prime factor of an 128-bit integer 
#[no_mangle]
pub const extern "C" fn get_factor_128(n: u128) -> u128{
    let inv = machine_prime::mul_inv2_128(n);
   let one = machine_prime::one_mont_128(n);
   let base = machine_prime::two_mont_128(one,n);
   
   match pollard_brent_128(base,inv,one,n){
      Some(factor) => return factor,
      None => {
      // if x^2 -1 failed try x^2+1
      // No particular reason except to reuse some values 
        let coef = machine_prime::to_mont_128(n-1,n);
        match pollard_brent_128(base,inv,coef,n){
           Some(factor) => return factor,
           None => {
             // Loop that has a roughly 0.5 probability of factoring each iteration
             // The probability of being unable to factor a composite is 1/2^period and the period is likely very large
             // So the risk of error is vanishingly small
            let mut param = drbg(n as u64);
              loop{
                 
                 let  rand_base= (param as u128)%(n-3)+3;
                match pollard_brent_128(rand_base,inv,one,n){
                   Some(factor) => return factor,
                   None => param=drbg(param),
                 }
              }
           }
        }
      }
   }
}

/// Factorization of a 128-bit integer
#[repr(C)]
pub struct Factorization128{
   pub factors: [u128;26],
   pub powers: [u8;26],
   pub len: usize,
}

impl Factorization128{
   const fn new() -> Self{
      Factorization128{factors: [0u128;26],powers: [0u8;26], len: 0}
   }
}


/// Complete factorization of 128-bit N
#[no_mangle]
pub const extern "C" fn factorize_128(mut n: u128) -> Factorization128{

      let mut t = Factorization128::new();
      
      let mut idx = 0usize;
      
      if n == 0{
         return t;
      }
      if n == 1{
         t.factors[0]=1;
         t.powers[1]=1;
         t.len = 1;
         return t;
      }
      
      let twofactor = n.trailing_zeros();
      if twofactor != 0{
        t.factors[0]=2u128;
        t.powers[0]=twofactor as u8;
        n>>=twofactor;
        idx+=1;
      }

      let mut i = 0usize;
      while i < 130 {
         let fctr = PRIMES_130[i] as u128;
            // strips out small primes
            if n % fctr == 0 {
                t.factors[idx]=fctr;
                let mut count = 0u8;
                while n % fctr == 0 {
                    count += 1;
                    n /= fctr;
                }
                
                t.powers[idx]=count;
                idx+=1;
            }
            i+=1;
        }
        
        if n == 1 {
            t.len=idx;
            return  t;
        }

        if machine_prime::is_prime_wc_128(n){
            t.factors[idx]=n;
            t.powers[idx]=1;
            idx+=1;
            t.len=idx;
            return  t;
        }
        while n != 1 {
            let k = get_factor_128(n);
            t.factors[idx]=k;
            let mut count = 0u8;
            while n % k == 0 {
                count += 1;
                n /= k;
            }
            t.powers[idx]=count;
            idx+=1;
            if n == 1{
               t.len=idx;
               return t;
            }
            if machine_prime::is_prime_wc_128(n){
               t.factors[idx]=n;
               t.powers[idx]=1;
               idx+=1;
               t.len=idx;
               return t;
            }
        }
        t
    }
