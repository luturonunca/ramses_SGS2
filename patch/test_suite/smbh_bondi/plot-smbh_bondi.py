from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.0
sizex = 18.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data = loadtxt('data1.dat')
rho  = data[:,2].reshape(128,128)

# Density
density = subplot(211)
density.imshow(rho,origin='lower',extent=[0.0,1.0,0.0,1.0],interpolation='None',cmap='plasma')
density.set_xlabel('Distance x (cm)')
density.set_ylabel('Distance y (cm)')

savefig('smbh_bondi.pdf',bbox_inches='tight')
