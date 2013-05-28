import cProfile
import time
import sys
import math
import multiprocessing
import select
import operator
import signal
from blumblumshub import BlumBlumShubRandom
from sieve import Sieve
from random import SystemRandom
random = SystemRandom()

class AllDone(Exception):
    pass

class TimeoutException(Exception):
    pass

def print_status(statuses):
    global lasttime
    now = time.time()
    intermediate_statuses = filter(lambda x: isinstance(x, tuple), statuses)
    loops = sum(map(operator.itemgetter(0), intermediate_statuses))
    p2_pass = sum(map(operator.itemgetter(1), intermediate_statuses))
    p1_pass = sum(map(operator.itemgetter(2), intermediate_statuses))
    done = reduce(operator.or_, map(operator.itemgetter(3), intermediate_statuses), False)
    result = filter(lambda x: isinstance(x, (int, long)), statuses)

    # binomial expansion with 4 terms
    candidate_probability = 1/(((math.log(2)*bits)**3)*sieve.advantage)
    prob = (loops * candidate_probability)
    prob -= float(loops*(loops-1))/2 * candidate_probability**2
    prob += float(loops*(loops-1)*(loops-2))/6 * candidate_probability**3
    prob -= float(loops*(loops-1)*(loops-2)*(loops-3))/24 * candidate_probability**4

    if done:
        print >>sys.stderr, "Found doubly safe prime after %s iterations (probability %s)." % (loops, prob)
        print >>sys.stderr, "We found %s primes and %s singly safe primes." % (p2_pass, p1_pass)
    if len(result) > 0:
        print result[0]
        raise AllDone
    elif now - lasttime > 3600:
        print >>sys.stderr, "%s candidates tested" % loops
        print >>sys.stderr, "%s primes found" % p2_pass
        print >>sys.stderr, "%s safe primes found" % p1_pass
        print >>sys.stderr, "%s expected probability of finding prime already" % prob
        print >>sys.stderr, "%s iterations per second" % (loops / (now - starttime))
        print >>sys.stderr, "%s effective guesses per second" % ((loops / (now - starttime)) / sieve.advantage)
        print >>sys.stderr, "%s seconds since last output" % (now - lasttime)
        print >>sys.stderr
        lasttime = now

    
if __name__ == "__main__":
    bits = 8192
    certainty = 256
    sieve = Sieve(max(bits-certainty,64))
    gsp = BlumBlumShubRandom.gen_special_prime

    print >>sys.stderr, "Trying to find a special prime with %s bits using a sieve of index %s, size %s, advantage %s" % (bits, sieve.index, math.log(sieve.modulus,2), sieve.advantage)

    # figure out how often to check in
    checkin_seconds = 5
    checkin_iterations = None
    def handler(*args):
        raise TimeoutException
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(checkin_seconds)
    while checkin_iterations is None:
        checkin_seconds *= 2
        try:
            def trial_callback(loops, p2_pass, p1_pass, done):
                global checkin_iterations
                checkin_iterations = loops
            print gsp(bits, certainty, random=random, callback=trial_callback, callback_period=1)
            sys.exit(0)
        except TimeoutException:
            pass
    print >>sys.stderr, "Subprocesses will check in every %s iterations (%s seconds)" % (checkin_iterations, checkin_seconds)
    signal.signal(signal.SIGALRM, signal.SIG_DFL)

    starttime = time.time()
    lasttime = starttime

    def callback(pipe):
        return (lambda *args: pipe.send(args))
    def worker(pipe):
        def thunk():
            signal.signal(signal.SIGINT, signal.SIG_IGN) # ignore Ctrl-C
            pipe.send(gsp(bits, certainty, random=random,
                          callback=callback(pipe), callback_period=checkin_iterations))
            sys.exit(0)
        return thunk

    parent_pipes = [None]*multiprocessing.cpu_count()
    processes = [None]*len(parent_pipes)
    epoll = select.epoll()
    for i in xrange(len(parent_pipes)):
        parent_pipes[i], child_pipe = multiprocessing.Pipe()
        epoll.register(parent_pipes[i].fileno(), select.EPOLLIN)
        processes[i] = multiprocessing.Process(target=worker(child_pipe))
        processes[i].daemon = True
        processes[i].start()
    del child_pipe

    statuses = [None]*len(parent_pipes)
    try:
        while True:
            try:
                epoll.poll()
            except IOError as e:
                import errno
                if e.errno == errno.EINTR: # Interrupted syscall
                    continue
                else:
                    raise
            for i,parent_pipe in enumerate(parent_pipes):
                if parent_pipe.poll():
                    statuses[i] = parent_pipe.recv()
            print_status(statuses)
    except (AllDone, KeyboardInterrupt):
        print >>sys.stderr, "Done! Cleaning up..."
        pass
    map(operator.methodcaller("terminate"), processes)

