import operator
import cPickle
import gmpy2
from itertools import *

primes = cPickle.Unpickler(open("10000_primes.pkl")).load()

def _successors(n, s):
    return [ (2*i + 1) % n for i in s ]

def _allowed_residues(n):
    """For a given number n, find the allowed residues mod n that will not obviously invalidate primality"""
    possibilities = range(n)
    first_residues = _successors(n, possibilities)
    second_residues = _successors(n, first_residues)
    return [ r for (i, r) in enumerate(possibilities)
             if r != 0
             if first_residues[i] != 0
             if second_residues[i] != 0 ]

def generate_naive(i):
    """
    Generate the exhaustive list of all residues mod the i'th primorial that do not obviously violate primality
    returns a tuple (i'th primorial, list of residues)
    """
    pairs = [ (n,set(_allowed_residues(n))) for n in primes[:i]]
    m = reduce(operator.mul, map(operator.itemgetter(0), pairs))
    residues = [ r for r in xrange(m)
                 if all(map(lambda (n,rs): (r % n) in rs, pairs)) ]
    return (m,residues)

def generate(i):
    """
    Generate a more efficient way of representing the allowed residues for the i'th primorial
    returns a tuple of (m, sieve)
    where sieve is a list of tuples (prime, allowed residues mod prime, expander)
    where expander is 1 mod the prime and 0 mod all the other primes
    """
    m = reduce(operator.mul, primes[:i])
    return (m, [ (n, _allowed_residues(n), ((m*int(gmpy2.invert(m/n,n)))/n) % m) for n in primes[:i] ] )

def _check(i):
    """Check the naive and efficient sieves against each other for a given i'th primorial"""
    sieve = generate_sieve(i)
    residues = list()
    for j in product(*map(operator.itemgetter(1), sieve[1])):
        residue = 0
        for (k, e) in izip(j, imap(operator.itemgetter(2), sieve[1])):
            residue += k*e
        residues.append(residue % sieve[0])
    residues.sort()

    (m, other_residues) = generate_naive(i)
    return (m == sieve[0] and residues == other_residues)
