import operator
import cPickle
import gmpy2
from itertools import *

primes = cPickle.Unpickler(open("10000_primes.pkl")).load()

def successors(n, s):
    return [ (2*i + 1) % n for i in s ]

def allowed_residues_prime(n):
    possibilities = range(n)
    first_residues = successors(n, possibilities)
    second_residues = successors(n, first_residues)
    return [ r for (i, r) in enumerate(possibilities)
             if r != 0
             if first_residues[i] != 0
             if second_residues[i] != 0 ]

def allowed_residues_primorial(i):
    pairs = [ (n,set(allowed_residues_prime(n))) for n in primes[:i]]
    m = reduce(operator.mul, map(operator.itemgetter(0), pairs))
    residues = [ r for r in xrange(m)
                 if all(map(lambda (n,rs): (r % n) in rs, pairs)) ]
    return (m,residues)

def fraction_passed(i):
    (m, rs) = allowed_residues_primorial(i)
    return float(len(rs))/(m-1)

def generate_sieve(i):
    m = reduce(operator.mul, primes[:i])
    return (m, [ (n, allowed_residues_prime(n), (m*int(gmpy2.invert(m/n,n))/n) % m) for n in primes[:i] ] )

def check_sieve(nth_primorial):
    sieve = generate_sieve(nth_primorial)
    residues = list()
    for i in product(*map(operator.itemgetter(1), sieve[1])):
        residue = 0
        for (j, e) in izip(i, imap(operator.itemgetter(2), sieve[1])):
            residue += j*e
        residues.append(residue % sieve[0])
    residues.sort()

    (m, other_residues) = allowed_residues_primorial(nth_primorial)
    return (m == sieve[0] and residues == other_residues)
