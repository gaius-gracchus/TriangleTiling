# -*- coding: UTF-8 -*-

"""Generate a triangle tiling
"""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import os

import numpy as np
import matplotlib.pyplot as plt

from TriangleTiling import (
  RegularMultiTriangleTiling, )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

OUTPUT_DIR = 'pentagon_animation'

m = 5
N = 100

os.makedirs( OUTPUT_DIR, exist_ok = True )

A_values = np.arange( 1, int( 360 / 5 ), 1 )

for A in A_values:

  tiling = RegularMultiTriangleTiling( A = A, m = m )
  tiling.plot(
    N = N,
    filename = os.path.join( OUTPUT_DIR, f'{A:02d}.png' ) )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#