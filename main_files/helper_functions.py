import numpy as np

def veclen(vec):
    return np.sqrt(np.dot(vec, vec))

def NormGeoSeries(x, N):
    """
    Returns numbers following the normalized geometric series
    """
    return (1/2)**x/np.sum((1/2)**np.arange(1, N+1))
