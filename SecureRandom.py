import random, struct, math

# TODO: re-implement random methods for correctness
class SecureRandom(random.Random):
    """
    A random number generator that gets its randomness from /dev/random.
    If you want to get randomness from /dev/urandom, use random.SystemRandom instead.
    """
    def __init__(self):
        self._file = None
        self._cache = 0L
        self._cache_len = 0L
        self.seed(None)

    def seed(self, ignore):
        if self._file is None:
            self._file = open('/dev/random', 'r')
        else:
            try:
                self._file.close()
            except:
                pass
            self._file = None
            self.seed(None)

    def getstate(self):
        return None
    def setstate(self, ignore):
        pass
    def jumpahead(self, ignore):
        pass

    def random(self):
        bits_per_float = 53
        randomness = self.getrandbits(bits_per_float)
        return float(randomness) / (2**bits_per_float)

    def getrandbits(self, n):
        bytes = max(int(math.ceil((n-self._cache_len) / 8.0)),0)

        if bytes >= 8:
            for _ in xrange(0,bytes,8):
                self._cache <<= 8*8
                self._cache |= struct.unpack('Q', self._file.read(8))[0]
                self._cache_len += 8*8
                bytes -= 8
        if bytes >= 4:
            for _ in xrange(0,bytes,4):
                self._cache <<= 8*4
                self._cache |= struct.unpack('L', self._file.read(4))[0]
                self._cache_len += 8*4
                bytes -= 4
        if bytes >= 2:
            for _ in xrange(0,bytes,2):
                self._cache <<= 8*2
                self._cache |= struct.unpack('H', self._file.read(2))[0]
                self._cache_len += 8*2
                bytes -= 2
        if bytes >= 1:
            for _ in xrange(0,bytes,1):
                self._cache <<= 8*1
                self._cache |= struct.unpack('B', self._file.read(1))[0]
                self._cache_len += 8*1
                bytes -= 1

        assert bytes == 0
        assert self._cache_len >= n
        retval = ((1<<n) - 1) & self._cache
        self._cache >>= n
        self._cache_len -= n
        return retval
