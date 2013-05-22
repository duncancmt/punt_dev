import primes
import SecureRandom
import math
import sys
import cProfile
import cPickle
random = SecureRandom.SecureRandom()

def gen_special_prime(bits, certainty=128, random=random):
    # In addition to the well-known constraint that p must be congruent
    # to 3, mod 4, in order to avoid short cycles, we impose the
    # additional constraints described here:
    # http://www.ciphersbyritter.com/NEWS6/BBS.HTM
    # Essentially, p must be a "doubly safe" prime: p1 := (p-1)/2, p2 := (p1-1)/2
    # where p1 and p2 are also prime

    # Only numbers that have a residue in the following set, mod the sieve modulus
    # can generate a doubly safe prime. This dramatically reduces our search space.
    # See sieve_generator.py for how these numbers were found
    # For example, two other possible choices of values are:
    # sieve_modulus = 210 # 4th primorial, advantage 0.03827751196172249
    # allowed_residues = [11, 41, 89, 131, 149, 179, 191, 209]
    # === OR ===
    # sieve_modulus = 2310 # 5th primorial, advantage 0.027717626678215677
    # allowed_residues = [41, 89, 131, 149, 179, 191, 221, 251, 359, 389, 419,
    #                     461, 509, 551, 569, 611, 641, 719, 779, 809, 821, 839,
    #                     851, 881, 971, 989, 1019, 1031, 1049, 1139, 1181, 1229,
    #                     1241, 1271, 1301, 1349, 1409, 1439, 1451, 1469, 1481,
    #                     1511, 1559, 1601, 1649, 1679, 1691, 1769, 1811, 1829,
    #                     1871, 1889, 1901, 1931, 1979, 2021, 2039, 2069, 2099,
    #                     2111, 2141, 2231, 2291, 2309]
    sieve_modulus = 223092870 # 9th primorial, advantage 0.012852046773166708
    allowed_residues = cPickle.Unpickler(open("allowed_residues.pkl")).load()
    assert len(allowed_residues) == 2867200 and sum(allowed_residues) == 319826631304040

    bits -= 2 # we'll get the low two bits by left shifting and adding one, twice
    rand_bits = bits
    rand_bits -= int(math.floor(math.log(sieve_modulus,2))) # we get these bits by multiplying by the modulus and adding the residue
    rand_bits -= 1 # we always set the high bit
    rand_bits = int(math.ceil(rand_bits))
    bits = int(math.ceil(bits))

    l = int(math.ceil(math.log(len(allowed_residues),2)))
    loops = 0
    p2_pass = 0
    p1_pass = 0
    while True:
        loops += 1
        if loops % 10**6 == 0:
            print "%s loops" % loops
            print "%s p2_pass" % p2_pass
            print "%s p1_pass" % p1_pass
        # randomly choose a residue
        # random.choice is implemented wrong, so we do it ourselves
        i = len(allowed_residues)
        while i >= len(allowed_residues):
            i = random.getrandbits(l)

        # choose a random p2 that has the chosen residue
        p2 = random.getrandbits(rand_bits)
        p2 |= (1 << (rand_bits))
        p2 *= sieve_modulus
        p2 += allowed_residues[i]            
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
    exec("print gen_special_prime(8192, 256)")
    # cProfile.run("print gen_special_prime(8192, 256)")
    
