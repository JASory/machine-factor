use crate::primes::PRIMES_130;

/// Deterministic Random Bit Generator (xorshift)
pub const fn drbg(mut x: u64) -> u64{
    x ^= x.wrapping_shr(12);
    x ^= x.wrapping_shl(25);
    x ^= x.wrapping_shr(27);
    x.wrapping_mul(0x2545F4914F6CDD1D)
}

/// GCD(A,B)
#[no_mangle]
pub const extern "C" fn gcd(mut a: u64, mut b: u64) -> u64{
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

const fn poly_eval(x: u64, subtrahend: u64, n: u64, npi: u64) -> u64{
     machine_prime::mont_sub(machine_prime::mont_prod(x,x,n,npi),subtrahend,n)
}


const fn pollard_brent(base: u64,inv:u64,subtrahend: u64, n: u64) -> Option<u64>{
    let m = 128;
    let mut r = 1;
    let mut q = 1;
    let mut g = 1;
    let mut ys = 1;
    let mut y = base;
    let mut x = y;
    let mut cycle = 0;
    
    while cycle < 17{
    
     cycle+=1;    
      x = y;
      
      let mut yloop = 0;
      
      while yloop < r{
      yloop+=1; 
        y = poly_eval(y,subtrahend,n,inv);      
      }
      
      let mut k = 0;
      
      loop{
      
      let mut i = 0;
      
       while i < m*cycle{
           if i >= r-k{
             break;
           }
         
         y=poly_eval(y,subtrahend,n,inv);
         q=machine_prime::mont_prod(q,x.abs_diff(y),n,inv);
         i+=1;
         } // end loop

         ys=y;
         g = gcd(q,n);
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
         ys=poly_eval(ys,subtrahend,n,inv);
         g=gcd(x.abs_diff(ys),n);
      
       }
    }
    if g !=1 && g !=n && machine_prime::is_prime_wc(g){
       return Some(g);
    }
    None
}

/// Returns some prime factor of an 64-bit integer
/// 
/// This function uses the Pollard-rho algorithm and mostly exists for FFI
#[no_mangle]
pub const extern "C" fn get_factor(n: u64) -> u64{
    let inv = machine_prime::mul_inv2(n);
   let one = machine_prime::one_mont(n);
   let base = machine_prime::two_mont(one,n);
   
   match pollard_brent(base,inv,one,n){
      Some(factor) => return factor,
      None => {
      // if x^2 -1 failed try x^2+1
      // No particular reason except to reuse some values 
        let coef = machine_prime::to_mont(n-1,n);
        match pollard_brent(base,inv,coef,n){
           Some(factor) => return factor,
           None => {
             // Loop that has a roughly 0.5 probability of factoring each iteration
            let mut param = drbg(n);
              loop{

                 let  rand_base= param%(n-3)+3;
                match pollard_brent(rand_base,inv,one,n){
                   Some(factor) => return factor,
                   None => param=drbg(param),
                 }
              }
           }
        }
      }
   }
}

/// Factorisation of an integer
///
/// Representation of the factors in the form p^k p2^k
#[repr(C)]
pub struct Factorization{
   pub factors: [u64;15],
   pub powers: [u8;15],
   pub len: usize,
}

impl Factorization{
   const fn new() -> Self{
      Factorization{factors: [0u64;15],powers: [0u8;15], len: 0}
   }
}



/// Complete factorization of N
#[no_mangle]
pub const extern "C" fn factorize(mut n: u64) -> Factorization{

      let mut t = Factorization::new();
      
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
        t.factors[0]=2u64;
        t.powers[0]=twofactor as u8;
        n>>=twofactor;
        idx+=1;
      }

      let mut i = 0usize;
      while i < 130 {
         let fctr = PRIMES_130[i] as u64;
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

        if machine_prime::is_prime_wc(n){
            t.factors[idx]=n;
            t.powers[idx]=1;
            idx+=1;
            t.len=idx;
            return  t;
        }
        while n != 1 {
            let k = get_factor(n);
            t.factors[idx]=k;
            //idx+=1;
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
            if machine_prime::is_prime_wc(n){
               t.factors[idx]=n;
               t.powers[idx]=1;
               idx+=1;
               t.len=idx;
               return t;
            }
          
        }
        t
    }
