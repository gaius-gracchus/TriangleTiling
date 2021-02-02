# -*- coding: UTF-8 -*-

"""Generate a triangle tiling
"""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import numpy as np

from TriangleTiling import (
  MultiTriangleTiling, )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

n = 3
m = 5
N = 100
C = 120

tiling = MultiTriangleTiling( C = C, n = n, m = m )

tiling.plot( N = N, filename = 'Figure_15_c.pdf' )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#