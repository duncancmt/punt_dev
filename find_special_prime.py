import cProfile
import time
import sys
import math
from blumblumshub import BlumBlumShubRandom
from sieve import Sieve
from random import SystemRandom
random = SystemRandom()

starttime = time.time()
lasttime = starttime
def callback(loops, p2_pass, p1_pass, done):
    global lasttime
    if done:
        print >>sys.stderr, "Found doubly safe prime after %s iterations." % loops
        print >>sys.stderr, "We found %s primes and %s singly safe primes." % (p2_pass, p1_pass)        
    else:
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

    
if __name__ == "__main__":
    bits = 8192
    certainty = 256
    command = "print BlumBlumShubRandom.gen_special_prime(%s, %s, random=random, callback=callback, callback_period=2**10)" % (bits, certainty)
    sieve = Sieve(max(bits-certainty,64))

    print >>sys.stderr, "Trying to find a special prime with %s bits using a sieve of index %s, size %s, advantage %s" % (bits, sieve.index, math.log(sieve.modulus,2), sieve.advantage)
    # exec(command)
    cProfile.run(command)

