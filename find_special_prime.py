import SecureRandom
import cPickle
import cProfile
import math
import operator
import primes
import sieve
import sys
import time
from random import SystemRandom
random = SystemRandom()

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
        if len(s) == 1:
            return s[0]
        l = int(math.ceil(math.log(len(s),2)))
        i = len(s)
        while i >= len(s):
            i = random.getrandbits(l)
        return s[i]
        
    def choose_residue():
        residue = 0
        for (n,s,e) in sieve_residues:
            residue = (residue + e*choice(s)) % sieve_modulus
        assert sieve_modulus > residue
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
    starttime = time.time()
    lasttime = starttime
    advantage = float(reduce(operator.mul, map(len, map(operator.itemgetter(1), sieve_residues)))) / sieve_modulus
    candidate_probability = 1/(((math.log(2)*bits)**3)*advantage)

    while True:
        loops += 1
        if loops % (2**15) == 0: # output roughly once an hour
            now = time.time()
            print >>sys.stderr, "%s candidates tested" % loops
            print >>sys.stderr, "%s primes found" % p2_pass
            print >>sys.stderr, "%s safe primes found" % p1_pass
            # binomial expansion with 4 terms
            prob = (loops * candidate_probability)
            prob -= float(loops*(loops-1))/2 * candidate_probability**2
            prob += float(loops*(loops-1)*(loops-2))/6 * candidate_probability**3
            prob -= float(loops*(loops-1)*(loops-2)*(loops-3))/24 * candidate_probability**4
            print >>sys.stderr, "%s expected probability of finding prime already" % prob
            print >>sys.stderr, "%s iterations per second" % (loops / (now - starttime))
            print >>sys.stderr, "%s seconds since last output" % (now - lasttime)
            print >>sys.stderr
            lasttime = now
            
        # choose a random p2 that has an allowed residue
        p2 = random.getrandbits(rand_bits)
        p2 |= (1 << (rand_bits))
        p2 *= sieve_modulus
        p2 += choose_residue()
        if p2 >> bits:
            continue

        p1 = p2 * 2 + 1
        p = p1 * 2 + 1
        # first run a few abbreviated Miller-Rabin tests to fail quickly
        if primes.mr_test(p2, rounds=1) \
           and primes.mr_test(p1, rounds=1) \
           and primes.mr_test(p, rounds=1):
            if primes.mr_test(p2, certainty=certainty):
                p2_pass += 1
                if primes.mr_test(p1, certainty=certainty):
                    p1_pass += 1
                    if primes.mr_test(p, certainty=certainty):
                        break
    print >>sys.stderr, "Found doubly safe prime after %s iterations." % loops
    print >>sys.stderr, "We found %s primes and %s singly safe primes." % (p2_pass, p1_pass)
    return p
    
if __name__ == "__main__":
    exec("print gen_special_prime(8192, 256)")
    # cProfile.run("print gen_special_prime(8192, 256)")
    
