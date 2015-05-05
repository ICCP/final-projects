# adopted from code given to me from W.J. Ong 
# for research I was doing at the time

import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker

# input file given in command line
ifile = sys.argv[1]

#read in data for contour plot
f = open(ifile,'r')
x, y, z = np.loadtxt(f, usecols=(0,1,2), unpack = True)

# define x and y grid
xset = set(x)
yset = set(y)
xval = list(xset)
yval = list(yset)
xval.sort()
yval.sort()
z2d = np.transpose(np.reshape(z, (len(xval),len(yval))))

#find absolute minimum
zmin = z2d.argmin()
idx = np.unravel_index(zmin, z2d.shape)
xmax = xval[idx[1]]
ymax = yval[idx[0]]

#normalize contour
#pick cmap
mnl = matplotlib.ticker.MaxNLocator(nbins=30)
levels = mnl.bin_boundaries(z2d.min(),z2d.max())
cmap = plt.get_cmap('Spectral')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

# set axis labels and title
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
strx = r'r (fm)' #x axis
stry = r'as (fm)' #y axis
strt = r'r-as grid, step six' #title

#plot graphs and global min (for given range)
plt.contour(xval,yval,z2d, 20, cmap = cmap)
plt.plot(xmax,ymax,'k*',markersize=10)
cbar = plt.colorbar()
cbar.set_label('$\chi^2$')
plt.xlim(min(xval),max(xval))
plt.ylim(min(yval),max(yval))
plt.xlabel(strx)
plt.ylabel(stry)
plt.title(strt)
plt.show()