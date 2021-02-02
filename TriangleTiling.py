# -*- coding: UTF-8 -*-

"""Classes for generating logarithmic triangle tilings, based on:

https://www.mathartfun.com/FathauerBridges2021v1.pdf
"""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.cm import get_cmap

from sympy.solvers import solve
from sympy import Symbol

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class MultiTriangleTiling( ):

  """Class to generate a triangle tiling with the specified parameters

  Parameters
  ----------
  C : float
    Angle of lower-left corner of triangle (in degrees)
  n : int
    Number of transformations of the starting triangle before a triangle
    partially shares an edge with the starting triangle
  m : int
    Number of arms in the tiling

  """

  #---------------------------------------------------------------------------#

  def __init__( self, C, n, m ):

    # store n and m as instance variables
    self.n = int( n )
    self.m = int( m )

    # compute all 3 angles in radians
    self.C = np.deg2rad( C )
    self.A = ( 2 * np.pi / self.m - np.pi + self.C ) / ( self.n - 1 )
    self.B = np.pi - ( self.A + self.C )

    # compute side lengths based on the Law of Sines
    self.a = np.sin( self.A ) / np.sin( self.C )
    self.b = np.sin( self.B ) / np.sin( self.C )
    self.c = 1

    # compute scaling factor
    self.get_s( )

  #---------------------------------------------------------------------------#

  def get_s( self ):

    """Compute the scaling factor `s`
    """

    # initialize scaling factor as unknown variable, assuming it's real and
    # greater than zero
    _s = Symbol( 's', real = True, positive = True )

    # solve for scaling factor (first argument is expression set equal to zero)
    s = solve( self.a * _s ** self.n + self.b * _s - 1, _s )

    # save result as float
    self.s = float( s[ 0 ] )

  #---------------------------------------------------------------------------#

  def next_triangle( self, triangle ):

    """Compute the next tesselation of a specified triangle

    Parameters
    ----------
    triangle : numpy.ndarray
      Array of shape ( 3, 2 ) specifying ( x, y ) coordinates of all 3 points
      of the triangle

    Returns
    -------
    numpy.ndarray
      Array of shape ( 3, 2 ) specifying ( x, y ) coordinates of all 3 points
      of the next triangle in the tesselation sequence, which has been
      sufficiently rotated, scaled, and shifted

    """

    # initialize rotation matrix to rotate by angle ``A``
    rot_mat = np.array( [
      [ np.cos( self.A ), np.sin( self.A ) ],
      [ -np.sin( self.A ), np.cos( self.A ) ] ] )

    # initialize next triangle and translate it such that point C is at the
    # origin
    _next_triangle = triangle.copy( )
    rot_shift = _next_triangle[ 0 ]
    _next_triangle -= rot_shift

    # rotate each point in the next triangle by angle ``A`` using the rotation
    # matrix
    for i in range( 3 ):
      _next_triangle[ i ] = np.dot( rot_mat, _next_triangle[ i ] )

    # scale the next triangle by the scaling factor
    _next_triangle *= self.s

    # translate the next triangle such that point C is at point B of the
    # original triangle
    _next_triangle += rot_shift
    _next_triangle += ( triangle[ 1 ] - triangle[ 0 ] )

    return _next_triangle

  #---------------------------------------------------------------------------#

  def get_triangles( self, N ):

    """Generate a sequence of tesselated triangles, for a single arm

    Parameters
    ----------
    N : int
      Number of triangle tiles to generate for each arm

    Returns
    -------
    numpy.ndarray
      Array of shape ( n, 3, 2 ) specifying ( x, y ) coordinates of all 3 points
      of all triangles in the tessellation sequence, for a single arm

    """

    # store N as an instance variable
    self.N = N

    # initialize array to store locations of points for all triangles in the
    # tessellation sequence
    self.triangles = np.zeros( ( self.N, 3, 2 ) )

    # define points of the first triangle in the tessellation sequence
    point_c = np.array( [ 0, 0 ] )
    point_b = self.a * np.array( [ np.cos( self.C ), np.sin( self.C ) ] )
    point_a = np.array( [ self.b, 0 ] )

    # stack the points into a single array of shape (3, 2 )
    triangle = np.vstack( [ point_c, point_b, point_a ] )

    # loop over the number of triangles in the sequence
    for i in range( self.N ):

      # store the points of the i-th triangle in the array
      self.triangles[ i ] = triangle

      # compute the next triangle in the tessellation sequence
      triangle = self.next_triangle( triangle = triangle )

      # shift the next triangle in the tessellation sequence such that its
      # point C is in the same location as point B of the previous triangle
      triangle += ( self.triangles[ i - 1, 1 ] - self.triangles[ 0, 0 ] )

  #---------------------------------------------------------------------------#

  def get_all_triangles( self, N ):

    """Generate a sequence of tesselated triangles, for all arms

    Parameters
    ----------
    N : int
      Number of triangle tiles to generate for each arm

    Returns
    -------
    numpy.ndarray
      Array of shape ( m, n, 3, 2 ) specifying ( x, y ) coordinates of all 3
      points of all triangles in the tessellation sequence, for all arms

    """

    # compute the array of points for all triangles in one arm of the
    # tessellation
    self.get_triangles( N = N )

    # initialize array of points for all triangles in ALL arms of the
    # tessellation, and store the triangles for the first arm in it
    self.all_triangles = np.zeros( ( self.m, self.N, 3, 2 ) )
    self.all_triangles[ 0 ] = self.triangles

    # compute the angle by which each arm is rotated
    self.theta_m = 2 * np.pi / self.m

    # loop over arms
    for i in range( 1, self.m ):

      # compute rotation angle for the given arm
      theta = self.theta_m * i

      # initialize rotation matrix for the arm rotation angle
      rot_mat = np.array( [
        [ np.cos( theta ), np.sin( theta ) ],
        [ -np.sin( theta ), np.cos( theta ) ] ] )

      # translate a copy of the first arm such that point A of the first
      # triangle in the arm is at the origin
      rot_shift = self.triangles[ 0, 2 ]
      new_arm = self.triangles.copy( )
      new_arm -= rot_shift

      # rotate all points in all triangles in the arm, by the
      # rotation angle
      for tri in range( N ):
        for point in range( 3 ):
          new_arm[ tri, point ] = np.dot( rot_mat, new_arm[ tri, point ] )

      # Translate the new arm such that point A of the first triangle in the
      # new arm is at point B in the ``n``-th triangle in the first arm
      new_arm += self.all_triangles[ i - 1, self.n - 1, 1 ]

      # store rotated arm in the i-th element of the ``all_triangles`` array`
      self.all_triangles[ i ] = new_arm

    return self.all_triangles

  #---------------------------------------------------------------------------#

  def plot( self, N, square = False, filename = None ):

    """Generate and save a visualization of the triangle tiling

    Parameters
    ----------
    N : int
      Number of triangle tiles to generate for each arm
    square : bool
      If ``True``, aspect ratio of the resulting figure is 1
    filename : str
      If not ``None``, the figure is saved to this file.
      File extension can be any extension supported by Matplotlib, including:
      {*.png, *.jpg, *.svg, *.pdf, *.ps}.

    """

    # compute all points for all triangles for all arms
    self.get_all_triangles( N = N )

    # compute spatial bounds of the tessellation
    x_min = np.min( self.all_triangles[ :, :, :, 0 ] ) * 1.02
    x_max = np.max( self.all_triangles[ :, :, :, 0 ] ) * 1.02
    y_min = np.min( self.all_triangles[ :, :, :, 1 ] ) * 1.02
    y_max = np.max( self.all_triangles[ :, :, :, 1 ] ) * 1.02

    # compute aspect ratio of visualization
    if square:
      ratio = 1
      y_max = y_min + ( x_max - x_min )
    else:
      ratio = ( y_max - y_min ) / ( x_max - x_min )

    # initialize figure and axis of visualization
    fig, ax = plt.subplots( figsize = ( 8, 8 * ratio ) )

    # get list of colors used to color each arm of the tessellation
    colors = get_cmap( 'tab10' ).colors

    # loop over arms
    for i , color in zip( range( self.m ), colors ):

      # initialize list of matplotlib Polygons
      patches = list( )

      # loop over all triangles in the arm
      for triangle in self.all_triangles[ i ]:

        # append a matplotlib Polygon of the given triangle to the the
        # polygon list
        polygon = Polygon(
          xy = triangle,
          closed = True,
          joinstyle = 'miter' )
        patches.append( polygon )

      # create a matplotlib PatchCollection of the list of polygons
      p = PatchCollection(
        patches = patches,
        color = color,
        edgecolor = 'k',
        linewidth = 1.0, )

      # add the PatchCollection to the figure axis
      ax.add_collection( p )

    # set bounds of visualization, and remove axis box and tickmarks
    ax.set_xlim( x_min, x_max )
    ax.set_ylim( y_min, y_max )
    ax.axis( 'off' )

    # ensure scaling is consistent between the two axes
    plt.subplots_adjust( 0, 0, 1, 1 )

    # if ``filename`` is specified, save to file, otherwise show`
    if filename is None:
      return ax
    else:
      plt.savefig( filename, bbox_inches = 'tight' )
      plt.close( )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class RegularMultiTriangleTiling( MultiTriangleTiling ):

  """Class to generate a triangle tiling with the specified parameters, for the
  special case of ``n`` equal to 1

  Parameters
  ----------
  A : float
    Angle of lower-right corner of triangle (in degrees)
  m : int
    Number of arms in the tiling

  """

  def __init__( self, A, m ):

    # store n and m as instance variables
    self.n = 1
    self.m = int( m )

    # compute all 3 triangle angles in radians
    self.A = np.deg2rad( A )
    self.C = np.pi - 2 * np.pi / self.m
    self.B = np.pi - ( self.A + self.C )

    # compute side lengths based on the Law of Sines
    self.a = np.sin( self.A ) / np.sin( self.C )
    self.b = np.sin( self.B ) / np.sin( self.C )
    self.c = 1

    # compute scaling factor
    self.get_s( )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#