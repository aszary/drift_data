import numpy as np
from scipy.optimize import leastsq

def least_sq(x, y, fun, v0, size=200, xmax=None):
    """
    no errors
    """
    x_0 = min(x)
    x_1 = max(x)
    ## Error function
    errfunc = lambda v, x, y: (fun(v, x) - y)
    res = leastsq(errfunc, v0, args=(np.array(x), np.array(y)), maxfev=10000, full_output=True)
    v, conv = res[0], res[1]
    if xmax is None:
        xmax = x_1
    x_new = np.linspace(x_0, xmax, size)
    y_new = fun(v, x_new)
    return x_new, y_new, v


# odr

def p3_rs(edot, nsp=10):
    res = 5.6 / nsp * (edot / 4e31) ** 0.5
    if res < 1:
        res = 1
    return res

def p3_rs_log10(edot_log10, nsp=10):
    return np.log10(p3_rs(10 ** edot_log10, nsp=nsp))

def p3_rs_log10n1(edot_log10, nsp_fun, v):
    nsp = nsp_fun(v, edot_log10)
    return np.log10(p3_rs(10 ** edot_log10, nsp=nsp))



def p3_log10n2(edot_log10, nsp_fun, v):
    nsp = nsp_fun(v, edot_log10)
    if edot_log10 < 32:
        return np.log10(10 ** 16.1955 * (10 ** edot_log10) **  -0.496491)
    else:
        return np.log10(p3_rs(10 ** edot_log10, nsp=nsp))
