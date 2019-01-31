# -*- coding: utf-8 -*-
"""
@author: Boxx
"""
import numpy as np
from scipy import interpolate
import pandas as pd


def radome(distance,
           height,
           radius,
           azCorrection,
           filename=None,
           threshold=5,
           resolution=360,
           plot=False):
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
        python dictionary containing azimuths array and elevations array

    Example:
        Find the obscura azimuth and elevations vectors for a 20 meter
        radius radome that is 60 meters away from the radiating source
        and has a center 6 meters above it.  The azimuth correction
        is 300 degrees.

    >>> import obscura_tools.radome_obscura as obscura
    >>> import pandas as pd
    >>> ans = obscura.radome(60, -6, 20, 300)
    >>> table = pd.DataFrame(ans)
    >>> table.head()

    ========= ========= ========= ========= =========
    <index>   Frequency Signal 1  Signal 2  Order
    ========= ========= ========= ========= =========
    0         1000.0    1.0       0.0       1.0
    1         1000.0    -1.0      1.0       2.0
    2         2000.0    0.0       1.0       1.0
    3         2000.0    2.0       0.0       2.0
    4         3000.0    1.0       1.0       2.0
    ========= ========= ========= ========= =========
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

    phi = np.linspace(0, pi / 2, N)
    theta = np.linspace(-pi / 2, pi / 2, N)

    el = []
    az = []
    point = {}
    index = 0

    for i in range(N):
        for j in range(N):
            x1 = R * sin(phi[i]) * cos(theta[j])
            y1 = R * sin(phi[i]) * sin(theta[j])
            z1 = R * cos(phi[i])

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
                point[index] = [x1, y1, z1]

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
                az.append(-atan2(dy, np.sqrt(dx * dx + dz * dz)) * 180 / pi + azCorrection)

                if az[index] > 359:
                    az[index] -= 360

                index += 1

    # Find values at integer azimuths using interpolation
    low = int(np.ceil(min(az)))  # Find lowest az
    high = int(np.ceil(max(az)))  # Find highest az
    azimuths = range(low, high)
    toInterp = interpolate.interp1d(az, el)
    elevations = toInterp(azimuths)

    obscuraFrame = pd.DataFrame({'azimuths': azimuths, 'elevations': elevations})
    if filename:
        obscuraFrame.to_csv(filename, index=False)

    if plot:
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from mpl_toolkits.mplot3d.axes3d import get_test_data
        # This import registers the 3D projection, but is otherwise unused.
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

        phi = np.linspace(0, pi, N)
        phi, theta = np.meshgrid(phi, theta)
        X = R * sin(phi) * cos(theta)
        Y = R * sin(phi) * sin(theta)
        Z = R * cos(phi)

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        point[len(point)] = [x0, y0, z0]
        points = pd.DataFrame.from_dict(point, orient='index')
        #for i in range(len(point)):
        ax.scatter3D(points[0], points[1], points[2], 'x', linewidth=3)

        #ax.plot3D(x0, y0, z0, 'x', linewidth=3)
        if x0 > 2 * R:
            upp = np.ceil(x0 + 10)
            if upp % 2:
                upp += 1
            ax.set_xlim([0, upp])
            ax.set_ylim([-upp/2, upp/2])
            ax.set_zlim([-upp/2, upp/2])

        plt.show()

    return obscuraFrame
    # phi=linspace(0,pi,N);
    # [phi,theta]=meshgrid(phi,theta);
    # X=R*sin(phi).*cos(theta);
    # Y=R*sin(phi).*sin(theta);
    # Z=R*cos(phi);

    # surf(X,Y,Z);hold on;
    # [h,w,l] = size(point);
    # for i = 1:h
    #     plot3(point(i,1),point(i,2),point(i,3),'x','linewidth',3)
    # end

    # plot3(x0,y0,z0,'x','linewidth',3)
    # if (x0>2*R)
    #     upp = ceil(x0+10);
    #     if mod(upp,2) > 0
    #         upp = upp + 1;
    #     end
    #     axis([0 upp -upp/2 upp/2 -upp/2 upp/2]);
    # else

    # en

    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots()
    #ax.plot(ans.azimuths, ans.elevations)
    #ax.grid(True, which="both")
    #ax.minorticks_on()
    #ax.set_title("Radome Obscura")
    #ax.set_xlabel("Azimuth (deg)")
    #ax.set_ylabel("Elevations (deg)")
    #plt.show()
#radome(60, -6, 20, 300, plot=True)
