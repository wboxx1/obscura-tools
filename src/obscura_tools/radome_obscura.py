# -*- coding: utf-8 -*-
"""
@author: Boxx
"""
import numpy as np
from scipy import interpolate


def radome(distance,
           height,
           radius,
           azCorrection,
           filename=None,
           threshold=5,
           resolution=360):
    """Create azimuth and elevations for radome obscura

    Creates azimuth and elevation arrays for a radome obscura defined
    by distance, height, radius, azCorrection.
    Accepts a filename for csv storage.  Does not save file if none is given.
    Optional tuning parameters are threshold, and resolution.

    Args:
        distance (double): the distance from the transmitter to the
            radome center in meters
        height (double): the height difference between center of beam
            and radome center in meters (value is negative
            if radome center is above center of beam)
        radius (double): radius of the obscura radome in meters
        azCorrection (int): azimuth correction in degrees

    Kwargs:
        filename (string): filename to store array
        threshold (int): error correction threshold.  default = 5
        resolution (int): determines the number of points in the meshgrid as
            res^2. default = 360
        plot (bool): set to true to plot radome meshgrid. default = False

    Returns:
     python dictionary containing azimuths array, elevations array and points array

    Example:
        Find the obscura azimuth and elevations vectors for a 20 meter
        radius radome that is 60 meters away from the radiating source
        and has a center 6 meters above it.  The azimuth correction
        is 300 degrees.

    >>> import obscura_tools.radome_obscura as obscura
    >>> import pandas as pd
    >>> ans = obscura.radome(60, -6, 20, 300)
    >>> table = pd.DataFrame({'azimuths': ans['azimuths'], 'elevations': ans['elevations']})
    >>> table.head()

    ==========  ==========  ==========
       <index>    azimuths  elevations
    ==========  ==========  ==========
    0                  281    9.028552
    1                  282   12.343207
    2                  283   14.459659
    3                  284   16.107143
    4                  285   17.472889
    ==========  ==========  ==========

        Use the points key to plot a 3D representation of the obscura.

    >>> import matplotlib.pyplot as plt
    >>> from matplotlib import cm
    >>> from mpl_toolkits.mplot3d.axes3d import get_test_data
    >>> from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    >>> point = ans['points']
    >>> # numpy variables
    >>> pi = np.pi
    >>> sin = np.sin
    >>> cos = np.cos

    >>> # variable simplifications
    >>> R = radius
    >>> N = 360

    >>> theta = np.linspace(0, pi, N)
    >>> phi = np.linspace(-pi / 2, pi / 2, N)
    >>> theta, phi = np.meshgrid(theta, phi)
    >>> X = R * sin(theta) * cos(phi)
    >>> Y = R * sin(theta) * sin(phi)
    >>> Z = R * cos(theta)

    >>> fig = plt.figure()
    >>> ax = fig.gca(projection='3d')
    >>> surf = ax.plot_surface(X, Y, Z) #, rstride=1, cstride=1, cmap=cm.coolwarm,
                       # linewidth=0, antialiased=False)

    >>> points = pd.DataFrame.from_dict(point, orient='index')
    >>> ax.scatter3D(points[0], points[1], points[2], c='red', linewidth=3)

    >>> x0 = points[0][0]
    >>> if x0 > 2 * R:
    >>>    upp = np.ceil(x0 + 10)
    >>>    if upp % 2:
    >>>        upp += 1
    >>>    ax.set_xlim([0, upp])
    >>>    ax.set_ylim([-upp/2, upp/2])
    >>>    ax.set_zlim([-upp/2, upp/2])

    >>> plt.show()

    .. image:: _static/radome_obscura.png
    """
    # thresh = threshold

    # Numpy definitions
    pi = np.pi
    sin = np.sin
    cos = np.cos
    atan2 = np.arctan2

    N = resolution
    x0 = distance  # Distance from transmitter to radome center
    y0 = 0         # Will always be zero for simplicity
    z0 = height    # Height difference between center of beam and radome center
    R = radius     # Radius of obscura radome

    theta = np.linspace(0, pi / 2, N)
    phi = np.linspace(-pi / 2, pi / 2, N)

    el = []
    az = []
    points = {0: [x0, y0, z0]}
    index = 1

    for i in range(N):
        for j in range(N):
            x1 = R * sin(theta[i]) * cos(phi[j])
            y1 = R * sin(theta[i]) * sin(phi[j])
            z1 = R * cos(theta[i])

            # Evaluate gradient at point
            dx = 2 * x1
            dy = 2 * y1
            dz = 2 * z1

            # Evaluate dot product
            dp = x1 * dx + y1 * dy + z1 * dz

            # Evaluate transmit point
            tp = x0 * dx + y0 * dy + z0 * dz

            if abs(dp - tp) < threshold:
                # Create array of points
                points[index] = [x1, y1, z1]

                # Establish point 1 and point2
                dx = x0 - x1
                dy = y0 - y1
                dz = z0 - z1

                # for elevation, take projections on xz plane
                # calculate arctan of z/x
                el.append(-atan2(dz, np.sqrt(dx * dx + dy*dy)) * 180 / pi)

                # for azimuth, take projections on xy plane
                # calculate arctangent of y/x
                # include azimuth correction
                newAz = -atan2(dy, np.sqrt(dx * dx + dz * dz)) * 180 / pi + azCorrection
                az.append(newAz % 360)

                index += 1

    # Find values at integer azimuths using interpolation
    low = int(np.ceil(min(az)))  # Find lowest az
    high = int(np.ceil(max(az)))  # Find highest az
    azimuths = range(low, high)
    toInterp = interpolate.interp1d(az, el)
    elevations = toInterp(azimuths)

    #obscuraFrame = pd.DataFrame({'azimuths': azimuths, 'elevations': elevations})
    if filename:
        #obscuraFrame.to_csv(filename, index=False)
        with open(filename, 'w') as csvFile:
            csvFile.write("azimuths,elevations\n")
            for i in range(len(azimuths)):
                csvFile.write("%s,%s\n" % (azimuths[i], elevations[i]))

    return {'azimuths': azimuths, 'elevations': elevations, 'points': points}
