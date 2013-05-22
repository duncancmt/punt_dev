import primes
import SecureRandom
import math
import sys
import cProfile
import cPickle
import sieve
random = SecureRandom.SecureRandom()

try:
    from gmpy2 import mpz
    has_gmpy = True
except ImportError:
    try:
        from gmpy import mpz
        has_gmpy = True
    except ImportError:
        has_gmpy = False


def gen_special_prime(bits, certainty=128, random=random):
    # In addition to the well-known constraint that p must be congruent
    # to 3, mod 4, in order to avoid short cycles, we impose the
    # additional constraints described here:
    # http://www.ciphersbyritter.com/NEWS6/BBS.HTM
    # Essentially, p must be a "doubly safe" prime: p1 := (p-1)/2, p2 := (p1-1)/2
    # where p1 and p2 are also prime

    # Only numbers with particular residues can generate a doubly safe prime
    # by only generating candidates with those particular residues, we
    # dramatically reduce our search space.
    (sieve_modulus, sieve_residues) = sieve.generate(15) # chosen to fit in a single limb of a gmpy2 mpz

    if has_gmpy:
        m = mpz(sieve_modulus)
        sieve_residues = [ (mpz(n), map(mpz, s), mpz(e)) for (n,s,e) in sieve_residues ]
    else:
        m = sieve_modulus

    def choice(s):
        # random.choice is implemented wrong, so we do it ourselves
        l = int(math.ceil(math.log(len(s),2)))
        i = len(s)
        while i >= len(s):
            i = random.getrandbits(l)
        return s[i]
        
    def choose_residue():
        residue = 0
        for (n,s,e) in sieve_residues:
            residue += e*choice(s) % sieve_modulus
        return int(residue)
    
    bits -= 2 # we'll get the low two bits by left shifting and adding one, twice
    rand_bits = bits
    rand_bits -= int(math.floor(math.log(sieve_modulus,2))) # we get these bits by multiplying by the modulus and adding the residue
    rand_bits -= 1 # we always set the high bit
    rand_bits = int(math.ceil(rand_bits))
    bits = int(math.ceil(bits))

    loops = 0
    p2_pass = 0
    p1_pass = 0
    while True:
        loops += 1
        if loops % 10**6 == 0:
            print "%s loops" % loops
            print "%s p2_pass" % p2_pass
            print "%s p1_pass" % p1_pass
            
        # choose a random p2 that has an allowed residue
        p2 = random.getrandbits(rand_bits)
        p2 |= (1 << (rand_bits))
        p2 *= sieve_modulus
        p2 += choose_residue()
        if p2 >> bits:
            continue

        p1 = p2 * 2 + 1
        p = p1 * 2 + 1
        if primes.mr_test(p2, certainty=certainty):
            p2_pass += 1
            if primes.mr_test(p1, certainty=certainty):
                p1_pass += 1
                if primes.mr_test(p, certainty=certainty):                    
                    break
    print "Found prime after %s iterations." % loops
    print "We found %s primes and %s singly safe primes." % (p2_pass, p1_pass)
    return p
    
if __name__ == "__main__":
    exec("print gen_special_prime(256, 256)")
    # cProfile.run("print gen_special_prime(8192, 256)")
    
