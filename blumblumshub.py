from __future__ import division

import decimal
from decimal import Decimal

import math
import random
import struct
import cPickle
from fractions import gcd

import primes
from memoize import memoize
from sieve import Sieve

try:
    from gmpy2 import mpz
    has_gmpy = True
except ImportError:
    try:
        from gmpy import mpz
        has_gmpy = True
    except ImportError:
        has_gmpy = False

def lcm(a,b):
    return (a * b) // gcd(a,b)    
    
class BlumBlumShubRandom(random.Random):
    def __init__(self, security, bits_to_supply=2**20, tolerance=0.5,
                 seed=None, paranoid=True, _state=None):
        """
        Arguments:
            security: The log (base 2) of the number of operations required to distinguish BBS output from random noise.
            tolerance: The difference between the probability that the sequence produced by BBS be declared non-random and the probability that a truly random sequence will be declared non-random.
            bits_to_supply: The number of bits that this RNG will be expected to produce securely. Default: 2**20
            seed: Deterministically initializes the RNG. Must be a long or integer with length 2*security. (Optional)
            paranoid: If paranoid is true, this RNG will raise a RuntimeError if asked to supply >bits_to_supply bits
            _state: If supplied, we skip normal initialization and use the supplied state instead. (Optional)
        """
        if _state is not None:
            self.setstate(_state)
        else:
            self.security = security
            self.paranoid = paranoid
            self._bit_count = 0L
            self._bit_count_limit = bits_to_supply
            self._cache = 0L
            self._cache_len = 0L
            
            (self._modulus_length, self._bits_per_iteration) = self.optimal_parameters(self.security, self._bit_count_limit, tolerance)

            # python random.SystemRandom() is unsuitable because it uses /dev/urandom
            if seed is None:
                seed = 0L
                with open('/dev/random','r') as randfile:
                    for i in xrange(int(math.ceil(2*self.security / 8.0))):
                        seed |= struct.unpack('B', randfile.read(1))[0] << i*8
            self.seed(seed)

    @classmethod
    def from_state(cls, state):
        """Create a BlumBlumShubRandom instance from a state tuple"""
        return cls(security=None, bits_to_supply=None, tolerance=None,
                   seed=None, paranoid=None, _state=state)

    def random(self):
        return float(self.getrandbits(53)) / 2**53

    def getrandbits(self,n):
        if (not self.paranoid) or self._bit_count + n < self._bit_count_limit:
            self._bit_count += n
        else:
            raise RuntimeError("Too many bits supplied. Security guarantees invalidated.")
        
        iterations = max(int(math.ceil((1.0*n-self._cache_len) / self._bits_per_iteration)),0)

        mask = 2**self._bits_per_iteration - 1
        for i in xrange(iterations):
            self._cache <<= self._bits_per_iteration
            self._cache |= self.state & mask
            self._cache_len += self._bits_per_iteration
            self.state = pow(self.state, 2, self.modulus)

        result = self._cache & ((1<<n) - 1)
        self._cache >>= n
        self._cache_len -= n
        return result

    def seed(self, seed):
        assert isinstance(seed, (int, long))
        self.state = x_0 = seed >> self.security
        random.seed(seed & (2**self.security - 1))

        p = self.gen_special_prime(self._modulus_length / 2, certainty=self.security, random=random)
        q = self.gen_special_prime(self._modulus_length / 2, certainty=self.security, random=random)
        random.seed()

        assert p != q
        self.modulus = p*q
        self.skip_modulus = lcm(p-1, q-1) # TODO: it probably isn't safe to keep this around

        if has_gmpy:
            self.modulus = mpz(self.modulus)
            self.skip_modulus = mpz(self.skip_modulus)
            self.state = mpz(self.state)

        # degeneracy test as described by
        # http://www.ciphersbyritter.com/NEWS6/BBS.HTM
        self.jumpahead(1)
        assert self.state > 0
        x_1 = self.state
        self.jumpahead(1)
        x_2 = self.state
        assert x_0 != x_1 and x_1 != x_2 and x_2 != x_0, \
               "Given parameters resulted in a degenerate cycle. Choose another initial state."

    def getstate(self):
        return (self.security, self.paranoid,
                self._bit_count, self._bit_count_limit,
                self._cache, self._cache_len,
                self._modulus_length, self._bits_per_iteration,
                int(self.modulus), int(self.state),
                int(self.skip_modulus) if self.skip_modulus is not None else None)

    def setstate(self, state):
        (self.security, self.paranoid,
         self._bit_count, self._bit_count_limit,
         self._cache, self._cache_len,
         self._modulus_length, self._bits_per_iteration,
         self.modulus, self.state, self.skip_modulus) = state

        if has_gmpy:
            self.modulus = mpz(self.modulus)
            self.state = mpz(self.state)
            if self.skip_modulus is not None:
                self.skip_modulus = mpz(self.skip_modulus)

    def jumpahead(self, n):
        if self.skip_modulus is not None:
            self.state = pow(self.state, pow(2, n, self.skip_modulus), self.modulus)
        else:
            for _ in xrange(n):
                self.state = pow(self.state, 2, self.modulus)

    @staticmethod
    def gen_special_prime(bits, certainty=128, random=random,
                          callback=lambda loops, p2_pass, p1_pass, done: None,
                          callback_period=2**13): # run callback roughly once an hour
        # In addition to the well-known constraint that p must be congruent
        # to 3, mod 4, in order to avoid short cycles, we impose the
        # additional constraints described here:
        # http://www.ciphersbyritter.com/NEWS6/BBS.HTM
        # Essentially, p must be a "doubly safe" prime: p1 := (p-1)/2, p2 := (p1-1)/2
        # where p1 and p2 are also prime

        # Only numbers with particular residues can generate a doubly safe prime.
        # By only generating candidates with those particular residues, we
        # dramatically reduce our search space.
        sieve = Sieve(max(bits-certainty,64))

        loops = 0
        p2_pass = 0
        p1_pass = 0

        # TODO: adjust certainty to account for the number of tests that we run
        while True:
            loops += 1
            if loops % callback_period == 0:
                callback(loops, p2_pass, p1_pass, False)

            # choose a random p2 that has an allowed residue
            p2 = sieve.make_candidate(bits - 2, random) # we get the low two bits by left shifting and adding one, twice
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
        callback(loops, p2_pass, p1_pass, True)
        return p
    
    @staticmethod
    def optimal_parameters(log_Ta, M, epsilon, maxn=2**32):
        """
        Arguments:
            log_Ta: Log (base 2) of the number of operations required to distinguish BBS output from random noise.
            M: Number of bits that this instance of BBS will output. If more than M bits are extracted from this instance, security guarantees are invalidated.
            epsilon: Also called the tolerance. The difference between the probability that the sequence produced by BBS be declared non-random and the probability that a truly random sequence will be declared non-random.
            maxn: Maximum bit-length of the modulus. Default: 2**32
        Returns:
            A tuple of (bit-length of modulus (n), bits to extract per cycle (j))
        """
    
        log_Ta = Decimal(log_Ta)
        M = Decimal(M)
        epsilon = Decimal(epsilon)
        maxn = Decimal(maxn)
    
        @memoize
        def attack_difficulty(n, j):
            n = Decimal(n)
            j = Decimal(j)
            ln2 = Decimal(2).ln()
            
            # number of operations required by GNFS to factor an integer of bit length n
            # Equation (8) in http://www.win.tue.nl/~berry/papers/ima05bbs.pdf recommends
            # gamma = 2.8e-3. However, ECRYPT is a better authority on the subject.
            # Section 6.2.1 of ECRYPT 2012 suggests gamma = 2**-14 to 2**-12
            # See: http://www.ecrypt.eu.org/documents/D.SPA.20.pdf
            gamma = Decimal(1)/2**14
            tmp1 = (Decimal(64)/9)**(Decimal(1)/3)
            tmp2 = (n*ln2)**(Decimal(1)/3)
            tmp3 = ((n*ln2).ln())**(Decimal(2)/3)
            L = gamma*(tmp1*tmp2*tmp3).exp()

            # difficulty of distinguishing BBS from a random sequence
            # See equation (9) in http://www.win.tue.nl/~berry/papers/ima05bbs.pdf
            delta = epsilon/((2**j - 1) * M)
            left = (L * delta**2)/(36 * n * (n.ln()/ln2))
            right = 2**(2*j+9) * n / (delta**4)
            return left - right

        @memoize
        def find_n(j):
            upper = maxn
            lower = 2**j
            n = None
            while upper > lower:
                n = (upper + lower)/2
                n = n.to_integral_value(decimal.ROUND_CEILING)
                if attack_difficulty(n, j) > 2**log_Ta:
                    if n == upper:
                        # because of the ceiling, we can get stuck here
                        break
                    upper = n
                else:
                    lower = n
                    
            n += n % 2 # force n to be even
            assert 2**log_Ta <= attack_difficulty(n, j)
            return n

        @memoize
        def work_per_bit(j):
            # See section 7.3 in http://www.win.tue.nl/~berry/papers/ima05bbs.pdf
            n = find_n(j)
            return M*n**2/j
    
        with decimal.localcontext() as ctx:
            ctx.prec = 200
            maxj = (maxn.ln()/Decimal(2).ln()).to_integral_value(decimal.ROUND_FLOOR)
            upper = maxj+1 # we can never actually select this value because of the floor
            lower = Decimal(1)
            j = None
            # valley finding
            while upper > lower:
                j = (upper + lower)/2
                j = j.to_integral_value(decimal.ROUND_FLOOR)
                work = work_per_bit(j)
                if work_per_bit(j-1) < work:
                    upper = j
                elif work_per_bit(j+1) < work:
                    if lower == j:
                        # because of the floor, we can get stuck here
                        break
                    lower = j
                else:
                    break

            n = find_n(j)
            assert 2**log_Ta <= attack_difficulty(n, j)
            assert work_per_bit(j-1) >= work_per_bit(j) or j == maxj
            assert work_per_bit(j+1) >= work_per_bit(j) or j == maxj
            n = int(n)
            j = int(j)
            return (n, j)
