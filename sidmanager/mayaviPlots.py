import numpy as np
from mayavi.mlab import *

def test_flow():
    x, y, z = np.mgrid[0:5, 0:5, 0:5]
    r = np.sqrt(x**2 + y**2 + z**4)
    u = y*np.sin(r)/r
    v = -x*np.sin(r)/r
    w = np.zeros_like(z)
    obj = flow(u, v, w)
    #obj.show()

test_flow()
