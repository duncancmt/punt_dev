import math
from decimal import Decimal

def decimalify(f, x_type):
    def decorated(x, *args, **kwargs):
        return Decimal(f(x_type(x), *args, **kwargs))
    return decorated

def get_bracket_info(f, brack=None, args=tuple(), kwargs=dict()):
    if brack is None:
        xa, xb, xc, fa, fb, fc, funcalls = bracket(f, args=args, kwargs=kwargs)
    elif len(brack) == 2:
        if brack[0] > brack[1]:
            xc = brack[0]
            xa = brack[1]
        else:
            xa = brack[0]
            xc = brack[1]
        print xa
        print xc
        def clamp(x):
            if x < xa:
                return xa
            if x > xc:
                return xc
            return x
        xa, xb, xc, fa, fb, fc, funcalls = bracket(f, xa=xa, xb=xc, clamp=clamp, args=args, kwargs=kwargs)
        print "	bracket: [(%s,%s), (%s,%s), (%s,%s)]" % (xa, fa, xb, fb, xc, fc)
    elif len(brack) == 3:
        xa, xb, xc = brack
        if (xa > xc):  # swap so xa < xc can be assumed
            dum = xa; xa = xc; xc = dum
        if not ((xa < xb) and (xb < xc)):
            raise ValueError("Not a bracketing interval.")
        fa = f(xa, *args, **kwargs)
        fb = f(xb, *args, **kwargs)
        fc = f(xc, *args, **kwargs)
        if not ((fb < fa) and (fb < fc)):
            raise ValueError("Not a bracketing interval.")
        funcalls = 3
    else:
        raise ValueError("Bracketing interval must be " \
                         "length 2 or 3 sequence.")
    return xa, xb, xc, fa, fb, fc, funcalls

def brent(f,brack,x_type=Decimal,rtol=1.48e-8,atol=1.0e-11,maxiters=500,args=tuple(),kwargs=dict(),full_output=False):
    #set up for optimization
    f = decimalify(f, x_type)
    rtol = Decimal(rtol)
    atol = Decimal(atol)
    xa, xb, xc, fa, fb, fc, funcalls = get_bracket_info(f,brack,args=args,kwargs=kwargs)
    xa = Decimal(xa)
    xb = Decimal(xb)
    xc = Decimal(xc)
    fa = Decimal(fa)
    fb = Decimal(fb)
    fc = Decimal(fc)
    cg = (3-Decimal(5).sqrt())/2
    #################################
    #BEGIN CORE ALGORITHM
    #we are making NO CHANGES in this
    #################################
    x = w = v = xb
    fw = fv = fx = f(x, *args, **kwargs)
    if (xa < xc):
        a = xa
        b = xc
    else:
        a = xc
        b = xa
    deltax = Decimal(0)
    funcalls += 1
    iters = 0
    while (iters < maxiters):
        tol1 = rtol*abs(x) + atol
        tol2 = 2*tol1
        xmid = (a + b)/2
        if abs(x - xmid) < (tol2 - (b - a)/2):  # check for convergence
            break
        if (abs(deltax) <= tol1):
            if (x >= xmid): deltax = a - x       # do a golden section step
            else: deltax = b - x
            rat = cg*deltax
        else:                              # do a parabolic step
            tmp1 = (x - w)*(fx - fv)
            tmp2 = (x - v)*(fx - fw)
            p = (x - v)*tmp2 - (x - w)*tmp1;
            tmp2 = 2*(tmp2 - tmp1)
            if (tmp2 > 0):
                p = -p
            tmp2 = abs(tmp2)
            dx_temp = deltax
            deltax = rat
            # check parabolic fit
            if ((p > tmp2*(a - x)) and (p < tmp2*(b - x)) and
                (abs(p) < abs(tmp2*dx_temp/2))):
                rat = p / tmp2        # if parabolic step is useful.
                u = x + rat
                if ((u - a) < tol2 or (b - u) < tol2):
                    if xmid - x >= 0: rat = tol1
                    else: rat = -tol1
            else:
                if (x >= xmid): deltax = a - x # if it's not do a golden section step
                else: deltax = b - x
                rat = cg*deltax

        if (abs(rat) < tol1):            # update by at least tol1
            if rat >= 0: u = x + tol1
            else: u = x - tol1
        else:
            u = x + rat
        fu = f(u, *args, **kwargs)      # calculate new output value
        funcalls += 1

        if (fu > fx):                 # if it's bigger than current
            if (u < x): a = u
            else: b = u
            if (fu <= fw) or (w == x):
                v = w; w = u; fv = fw; fw = fu
            elif (fu <= fv) or (v == x) or (v == w):
                v = u; fv = fu
        else:
            if (u >= x): a = x
            else: b = x
            v = w; w = x; x = u
            fv = fw; fw = fx; fx = fu

        iters += 1
    #################################
    #END CORE ALGORITHM
    #################################

    if full_output:
        return x, fx, iters, funcalls
    else:
        return x

    

def bracket(f, xa=0.0, xb=1.0, clamp=None, args=tuple(),kwargs=dict(), grow_limit=110.0, maxiters=1000):
    """
    Bracket the minimum of the function.

    Given a function and distinct initial points, search in the
    downhill direction (as defined by the initital points) and return
    new points xa, xb, xc that bracket the minimum of the function
    f(xa) > f(xb) < f(xc). It doesn't always mean that obtained
    solution will satisfy xa<=x<=xb

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function to minimize.
    xa, xb : float, optional
        Bracketing interval. Defaults `xa` to 0.0, and `xb` to 1.0.
    args : tuple, optional
        Additional arguments (if present), passed to `f`.
    grow_limit : float, optional
        Maximum grow limit.  Defaults to 110.0
    maxiters : int, optional
        Maximum number of iterations to perform. Defaults to 1000.

    Returns
    -------
    xa, xb, xc : float
        Bracket.
    fa, fb, fc : float
        Objective function values in bracket.
    funcalls : int
        Number of function evaluations made.

    """
    xa = Decimal(xa)
    xb = Decimal(xb)
    grow_limit = Decimal(grow_limit)
    gold = (Decimal(5).sqrt()+1)/2
    if clamp is None:
        clamp = lambda x: x
    
    fa = f(xa, *args, **kwargs)
    fb = f(xb, *args, **kwargs)
    if (fa < fb):                      # Switch so fa > fb
        dum = xa; xa = xb; xb = dum
        dum = fa; fa = fb; fb = dum
    xc = clamp(xb + gold*(xb - xa))
    fc = f(xc, *args, **kwargs)
    funcalls = 3
    iters = 0
    while (fc < fb):
        tmp1 = (xb - xa)*(fb - fc)
        tmp2 = (xb - xc)*(fb - fa)
        val = tmp2 - tmp1
        denom = 2*val
        w = clamp(xb - ((xb - xc)*tmp2 - (xb - xa)*tmp1) / denom)
        wlim = xb + grow_limit*(xc - xb)
        if iters > maxiters:
            raise RuntimeError("Too many iterations.")
        iters += 1
        if (w - xc)*(xb - w) > 0.0:
            fw = f(w, *args, **kwargs)
            funcalls += 1
            if (fw < fc):
                xa = xb; xb = w; fa = fb; fb = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            elif (fw > fb):
                xc = w; fc = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            w = clamp(xc + gold*(xc - xb))
            fw = f(w, *args, **kwargs)
            funcalls += 1
        elif (w - wlim)*(wlim - xc) >= 0.0:
            w = clamp(wlim)
            fw = f(w, *args, **kwargs)
            funcalls += 1
        elif (w - wlim)*(xc - w) > 0.0:
            fw = f(w, *args, **kwargs)
            funcalls += 1
            if (fw < fc):
                xb = xc; xc = w; w = clamp(xc + gold*(xc - xb))
                fb = fc; fc = fw; fw = f(w, *args, **kwargs)
                funcalls += 1
        else:
            w = clamp(xc + gold*(xc - xb))
            fw = f(w, *args, **kwargs)
            funcalls += 1
        xa = xb; xb = xc; xc = w
        fa = fb; fb = fc; fc = fw
    return xa, xb, xc, fa, fb, fc, funcalls
