"""
MOST OF THIS CODE WAS STOLEN FROM SciPy. THE LICENSE FOR SciPy IS REPRODUCED BELOW


Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2012 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
"""

from decimal import Decimal

def decimalify(f):
    def decorated(*args, **kwargs):
        return Decimal(f(*args, **kwargs))
    return decorated
        
def brent(f, brack, args=tuple(), kwargs=dict(), guess=None, tol=Decimal('1.48e-8'), maxiters=500):
    """
    Find the minimum of f(x, *args, **kwargs) on the interval brack[0] <= x <= brack[1]
    brent converts all values to decimal.Decimal objects before working on them
    values for x passed to f will be decimal.Decimal objects
    result is cast from decimal.Decimal to the same type as brack[0]

    Arguments:
        f: The objective function. Must have exactly one minimum in the range brack
        brack: The bounds on which to search. Supplied as a list-like object (supports subscripting).
        args: Additional positional arguments to f (optional)
        kwargs: Additional keyword arguments to f (optional)
        guess: Initial guess of where on the range the minimum lies. Default: (brack[0]+brack[1])/2
        tol: Terminate when the change in the objective function is below this value. Default: 1.48e-8
        maxiters: Maximum number of steps brent will take before terminating. Default: 500

    Returns:
        The x corresponding to the minimum value of f(x, *args, **kwargs)
    """

    # TODO: consider stealing the bracketing code from scipy
    if brack[0] <= brack[1]:
        a = brack[0]
        c = brack[1]
    else:
        a = brack[1]
        c = brack[0]
    a = Decimal(a)
    c = Decimal(c)
    if guess is None:
        b = (a+c)/2
    else:
        b = guess
    b = Decimal(b)
    if not (a <= b and b <= c):
        raise ValueError("Invalid bracketing interval")

    f = decimalify(f)
    cg = (3-Decimal(5).sqrt())/2
    deltax = Decimal(0)
    iters = 0

    x = w = v = b
    fx = fw = fv = f(x, *args, **kwargs)
    while(iters < maxiters):
        tol1 = tol * abs(x+1)
        tol2 = 2 * tol1
        xmid = (a + b)/2
        if abs(x - xmid) < (tol2 - (b - a)/2): # check for convergence
            xmin = x
            break
        if abs(deltax) <= tol1: # do a golden section step
            if x >= xmid:
                deltax = a - x
            else:
                deltax = b - x
            rat = cg * deltax
        else: # do a parabolic step
            tmp1 = (x - w)*(fx - fv)
            tmp2 = (x - v)*(fx - fw)
            p = (x - v)*tmp2 - (x - w)*tmp1;
            tmp2 = 2*(tmp2 - tmp1)
            if tmp2 > 0:
                p = -p
            tmp2 = abs(tmp2)
            dx_temp = deltax
            deltax = rat
            if p > tmp2*(a - x) and p < tmp2*(b - x) and abs(p) < abs((tmp2*dx_temp)/2): # check parabolic fit
                # parabolic step is useful
                rat = p / tmp2
                u = x + rat
                if (u - a) < tol2 or (b - u) < tol2:
                    if xmid - x >= 0:
                        rat = tol1
                    else:
                        rat = -tol1
            else:
                # parabolic step is not useful, do a golden section step
                if (x >= xmid):
                    deltax = a - x
                else:
                    deltax = b - x
                rat = cg*deltax
        if abs(rat) < tol1: # update by at least tol1
            if rat >= 0:
                u = x + tol1
            else:
                u = x - tol1
        else:
            u = x + rat

        fu = f(u, *args, **kwargs) # calculate new output value
        if fu > fx:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v = w
                w = u
                fv = fw
                fw = fu
            elif fu <= fv or v == x or v == w:
                v = u
                fv = fu
        else:
            if (u >= x):
                a = x
            else:
                b = x
            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu

        iters += 1

    return type(brack[0])(xmin)
