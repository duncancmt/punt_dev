import random, struct, math

class SecureRandom(random.Random):
    def __init__(self, paranoid=False):
        self._file = None
        self._paranoid = paranoid
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
        retval = 0L
        bytes = int(math.ceil(n / 8.0))
        for i in xrange(bytes):
            retval |= struct.unpack('B', self._file.read(1))[0] << i*8
        return retval >> (bytes * 8 - n)
