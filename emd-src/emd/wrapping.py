"""Wrapping

Helper module for plotting wrapped data
"""
import numpy as np


def interpolate_for_plot(x, y, low, high):
    """interpolate_for_plot(x, y, low, high)

    Returns two lists of x- and y- coordinates, where vertical lines due to
    wrapping in the y coordinate are removed.

    Parameters
    ----------
    x : array_like
      1d monotonic increasing x-coordinates
    y : array_like
      1d wrapped (low <= y < high) y-coordinates
    low : float
      lower wrapping boundary
    high : float
      higher wrapping boundary

    Returns
    -------
    ix : list
      monotonic increasing new x-coordinates
    iy : list
      interpolated new y-coordinates (with NaNs to prevent vertical lines)

    """
    # Test low and high boundaries
    assert np.isscalar(low) and np.isscalar(high)
    assert low < high
    interval = float(high - low)

    # Test x and y coordinates
    x = np.squeeze(x)
    y = np.squeeze(y)
    assert x.ndim == 1
    assert np.all(np.diff(x) > 0)
    assert y.ndim == 1
    assert np.all(np.logical_and(low <= y, y < high))

    # calculate points at which the data wraps
    dy = np.diff(y)
    wrapidx, = np.where(np.abs(dy) > (interval / 2))
    jumpsign = np.sign(dy[wrapidx])

    # Generate output arrays
    oarrx = []
    oarry = []
    oldidx = None
    for i, s in zip(wrapidx, jumpsign):
        # extend output with data interval which does not wrap
        oarrx.extend(x[oldidx:i + 1])
        oarry.extend(y[oldidx:i + 1])
        # linear interpolation at wrap point
        x0, x1 = x[i], x[i + 1]
        y0, y1 = y[i], y[i + 1] - s * interval
        m = (y0 - y1) / (x0 - x1)
        t = y0 - m * x0
        # new points added to draw line up to boundaries
        if m > 0:
            oarrx.extend([(high - t) / m] * 3)
            oarry.extend([high, np.nan, low])
        else:
            oarrx.extend([(low - t) / m] * 3)
            oarry.extend([low, np.nan, high])
        # next loop
        oldidx = i + 1
    # ensure that all data is returned
    oarrx.extend(x[oldidx:])
    oarry.extend(y[oldidx:])

    return oarrx, oarry


def wrap_to(x, low, high):
    """wrap_to(x, low, high)

    Returns x wrapped to the half open interval [low, high[

    """
    interval = float(high - low)
    return ((x - low) % interval) + low
