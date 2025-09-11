from scipy.interpolate import interpn
from matplotlib.colors import Normalize 
import numpy as np
from matplotlib import pyplot as plt

def density_scatter( x , y, ax = None, fig = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    im = ax.scatter( x, y, c=z, **kwargs )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = plt.colorbar(im)
    
   # cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax, cmap = "plasma")
    cbar.ax.set_ylabel('Density')
    cbar.ax.set_ylim(bottom = 0)

    return ax