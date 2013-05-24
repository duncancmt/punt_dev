import SecureRandom
import cPickle
import cProfile
import math
import operator
import primes
import sys
import time
from random import SystemRandom
from sieve import Sieve
random = SystemRandom()

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
    sieve_index = 15 # chosen to fit in a python int
    sieve = Sieve(sieve_index)
        
    loops = 0
    p2_pass = 0
    p1_pass = 0
    starttime = time.time()
    lasttime = starttime

    print >>sys.stderr, "Trying to find a special prime with %s bits using a sieve of index %s, size %s, advantage: %s" % (bits, sieve_index, math.log(sieve.modulus,2), sieve.advantage)
    while True:
        loops += 1
        if loops % (2**14) == 0: # output roughly once an hour
            now = time.time()
            print >>sys.stderr, "%s candidates tested" % loops
            print >>sys.stderr, "%s primes found" % p2_pass
            print >>sys.stderr, "%s safe primes found" % p1_pass
            # binomial expansion with 4 terms
            candidate_probability = 1/(((math.log(2)*bits)**3)*sieve.advantage)
            prob = (loops * candidate_probability)
            prob -= float(loops*(loops-1))/2 * candidate_probability**2
            prob += float(loops*(loops-1)*(loops-2))/6 * candidate_probability**3
            prob -= float(loops*(loops-1)*(loops-2)*(loops-3))/24 * candidate_probability**4
            print >>sys.stderr, "%s expected probability of finding prime already" % prob
            print >>sys.stderr, "%s iterations per second" % (loops / (now - starttime))
            print >>sys.stderr, "%s effective guesses per second" % ((loops / (now - starttime)) / sieve.advantage)
            print >>sys.stderr, "%s seconds since last output" % (now - lasttime)
            print >>sys.stderr
            lasttime = now
            
        # choose a random p2
        p2 = sieve.make_candidate(bits, random)
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
    nbits = 8192
    command = "print gen_special_prime(%s, 256)" % nbits
    exec(command)
    # cProfile.run(command)
    
