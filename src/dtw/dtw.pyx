## This code is written by Davide Albanese, <albanese@fbk.eu>
## (C) 2011 mlpy Developers.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from cdtw cimport *

np.import_array()


def dtw_std(x, y, dist_only=True, squared=False):
    """Standard DTW as described in [Muller07]_,
    using the Euclidean distance (absolute value 
    of the difference) or squared Euclidean distance
    (as in [Keogh01]_) as local cost measure.

    :Parameters:
       x : 1d array_like object (N)
          first sequence
       y : 1d array_like object (M)
          second sequence
       dist_only : bool
          compute only the distance
       squared : bool
          squared Euclidean distance

    :Returns:
       dist : float
          unnormalized minimum-distance warp path 
          between sequences
       cost : 2d numpy array (N,M) [if dist_only=False]
          accumulated cost matrix
       path : tuple of two 1d numpy array (path_x, path_y) [if dist_only=False]
          warp path
    
    .. [Muller07] M Muller. Information Retrieval for Music and Motion. Springer, 2007.
    .. [Keogh01] E J Keogh, M J Pazzani. Derivative Dynamic Time Warping. In First SIAM International Conference on Data Mining, 2001.
    """

    cdef np.ndarray[np.float_t, ndim=1] x_arr
    cdef np.ndarray[np.float_t, ndim=1] y_arr
    cdef np.ndarray[np.float_t, ndim=2] cost_arr
    cdef np.ndarray[np.int_t, ndim=1] px_arr
    cdef np.ndarray[np.int_t, ndim=1] py_arr
    cdef Path p
    cdef double dist
    cdef int i
    cdef int sq

    x_arr = np.ascontiguousarray(x, dtype=np.float)
    y_arr = np.ascontiguousarray(y, dtype=np.float)
    cost_arr = np.empty((x_arr.shape[0], y_arr.shape[0]), dtype=np.float)

    if squared: sq = 1
    else: sq = 0

    dist = std(<double *> x_arr.data, <double *> y_arr.data, 
                <int> x_arr.shape[0], <int> y_arr.shape[0],
                <double *> cost_arr.data, sq)
    if dist_only:
        return dist
    else:
        path(<double *> cost_arr.data, <int> cost_arr.shape[0], 
              <int> cost_arr.shape[1], -1, -1, &p)
        px_arr = np.empty(p.k, dtype=np.int)
        py_arr = np.empty(p.k, dtype=np.int)
        for i in range(p.k):
            px_arr[i] = p.px[i]
            py_arr[i] = p.py[i] 
        free (p.px)
        free (p.py)
        return dist, cost_arr, (px_arr, py_arr)

    
def dtw_subsequence(x, y):
    """Subsequence DTW as described in [Muller07]_,
    assuming that the length of `y` is much larger 
    than the length of `x` and using the Manhattan 
    distance (absolute value of the difference) as 
    local cost measure.

    Returns the subsequence of `y` that are close to `x` 
    with respect to the minimum DTW distance.
    
    :Parameters:
       x : 1d array_like object (N)
          first sequence
       y : 1d array_like object (M)
          second sequence

    :Returns:
       dist : float
          unnormalized minimum-distance warp path
          between x and the subsequence of y
       cost : 2d numpy array (N,M) [if dist_only=False]
          complete accumulated cost matrix
       path : tuple of two 1d numpy array (path_x, path_y)
          warp path

    """

    cdef np.ndarray[np.float_t, ndim=1] x_arr
    cdef np.ndarray[np.float_t, ndim=1] y_arr
    cdef np.ndarray[np.float_t, ndim=2] cost_arr
    cdef np.ndarray[np.int_t, ndim=1] px_arr
    cdef np.ndarray[np.int_t, ndim=1] py_arr
    cdef Path p
    cdef int i
    
    x_arr = np.ascontiguousarray(x, dtype=np.float)
    y_arr = np.ascontiguousarray(y, dtype=np.float)
    cost_arr = np.empty((x_arr.shape[0], y_arr.shape[0]), dtype=np.float)

    subsequence(<double *> x_arr.data, <double *> y_arr.data, 
                 <int> x_arr.shape[0], <int> y_arr.shape[0],
                 <double *> cost_arr.data)
    
    idx = np.argmin(cost_arr[-1, :])
    dist = cost_arr[-1, idx]

    subsequence_path(<double *> cost_arr.data, <int> x_arr.shape[0],
                      <int> y_arr.shape[0], <int> idx, &p)
        
    px_arr = np.empty(p.k, dtype=np.int)
    py_arr = np.empty(p.k, dtype=np.int)
    
    for i in range(p.k):
        px_arr[i] = p.px[i]
        py_arr[i] = p.py[i]
            
    free (p.px)
    free (p.py)

    return dist, cost_arr, (px_arr, py_arr)
