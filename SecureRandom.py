import random, struct, math

class SecureRandom(random.Random):
    def __init__(self, paranoid=False):
        self._file = None
        self._paranoid = paranoid
        self._cache = 0L
        self._cache_len = 0L
        self.seed(None)

    def seed(self, ignore):
        if self._file is None:
            try:
                close(self._file)
            except:
                pass
            if self._paranoid:
                fname = '/dev/random'
            else:
                fname = '/dev/urandom'
            self._file = open(fname, 'r')

    def getstate(self):
        return None
    def setstate(self, ignore):
        pass
    def jumpahead(self, ignore):
        pass
    def random(self):
        return abs(struct.unpack('l', self._file.read(4))[0])/(0.+(~(1<<31)))
    def getrandbits(self, n):
        bytes = max(int(math.ceil((1.0*n-self._cache_len) / 8.0)),0)

        for i in xrange(bytes):
            self._cache <<= 8
            self._cache |= struct.unpack('B', self._file.read(1))[0]
            self._cache_len += 8

        assert self._cache_len >= n
        retval = ((1<<n) - 1) & self._cache
        self._cache >>= n
        self._cache_len -= n
        return retval
