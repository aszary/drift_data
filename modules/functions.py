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
