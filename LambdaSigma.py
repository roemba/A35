#=================================================================
# AE2220-II - Computational Modelling.
# Analysis program for work session 2
#
# Line 20: Definition of gammas for 4-stage time march
# Line 18: Definition of the lambda-sigma relation
#
#=================================================================
import numpy as np
import scipy.linalg as spl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#------------------------------------------------------
# Input parameters
#------------------------------------------------------
nx    = 100;  # Number of mesh points (must be even)
alpha = 0.5;  # Courant number
kdt   = 0.04; # Artificial viscosity parameter * Delta t
g1    = 1.;   # LSRK : gamma 1
g2    = 0.;   # LSRK : gamma 2
g3    = 0.5;  # LSRK : gamma 3
g4    = 1.;   # LSRK : gamma 4

#------------------------------------------------------
# Function for the lambda-sigma relation
#------------------------------------------------------
def lamSig(ldt):
#  sigma = 1 + ldt;               # Euler explicit time march
#  sigma = 1/(1- ldt);            # Euler implicit time march
#  sigma = (1+ldt/2)/(1- ldt/2)   # Trapezoidal time march
   sigma = 1+g4*ldt+g4*g3*ldt**2+g4*g3*g2*ldt**3+g4*g3*g2*g1*ldt**4
#  sigma  = 1 + g1*ldt + g1*g2*ldt**2 + g1*g2*g3*ldt**3 + g1*g2*g3*g4*ldt**4 
   return sigma

#------------------------------------------------------
# Define the semi-discrete matrix A * Dt 
# for the linear advection operator 
#------------------------------------------------------
AaDt = np.zeros((nx, nx))
for i in xrange(0, nx):
   if i == 0:      # Left periodic boundary
      AaDt[i,nx-1] = -1;
      AaDt[i,i]    =  0;
      AaDt[i,i+1]  =  1;
   elif i == nx-1: # Right periodic boundary
      AaDt[i,i-1]  = -1;
      AaDt[i,i]    =  0;
      AaDt[i,0]    =  1;
   else :          # Interior
      AaDt[i,i-1]  = -1;
      AaDt[i,i]    =  0;
      AaDt[i,i+1]  =  1;
AaDt *= -alpha/2.;

#------------------------------------------------------
# Define the semi-discrete matrix A * Dt 
# for linear diffusion * (Delta x)^2 
# (This provides an artificial viscosity which 
#  vanishes with decreasing Delta x.  
#------------------------------------------------------
AdDt = np.zeros((nx, nx))
for i in xrange(0, nx):
   if i == 0:      # Left periodic boundary
      AdDt[i,nx-1] =  1;
      AdDt[i,i]    = -2;
      AdDt[i,i+1]  =  1;
   elif i == nx-1: # Right periodic boundary
      AdDt[i,i-1]  =  1;
      AdDt[i,i]    = -2;
      AdDt[i,0]    =  1;
   else :          # Interior
      AdDt[i,i-1]  =  1;
      AdDt[i,i]    = -2;
      AdDt[i,i+1]  =  1;
AdDt *= kdt;

#------------------------------------------------------------
# Define the total semi-discrete matrix A*DT, then
# compute the semi-discrete eigenvalues lambda *DT
# from the expression in the notes for circulant matrices
#------------------------------------------------------------
ADt = AaDt + AdDt;
beta=np.zeros(nx);
ldt=np.zeros(nx,'complex')
for m in xrange(nx):
  beta[m] = 2*np.pi*m/nx;
  if beta[m] > np.pi: beta[m] = 2*np.pi - beta[m]; # negative beta modes
  for j in xrange(0, nx):
    ldt[m] = ldt[m] + ADt[0,j]*np.exp(1j*2.*np.pi*j*m/nx);

#------------------------------------------------------------
# Compute the eigenvalues of C using the lambda-sigma relation, 
# then determine the amplitude and relative phase of each mode
#------------------------------------------------------------
sigma  = lamSig(ldt);
magSig = np.abs(sigma);
relPse = np.ones(nx);
for m in xrange(1,nx):
   relPse[m] = -np.angle(sigma[m])/(alpha*beta[m]); 
   if (m>nx/2): relPse[m] = -relPse[m] # negative beta modes


#===================================================================
# Output results
#===================================================================

#------------------------------------------------------------
# Write the results to the screen
#------------------------------------------------------------
print "mode   beta       ldt          sigma      magSig relPse"
print "----------------------------------------------------------"
for m in xrange(nx):
  print "%3i" % m, "| %5.2f" % beta[m],
  print "| %5.2f" % np.real(ldt[m]),   "%5.2f" % np.imag(ldt[m]),
  print "| %5.2f" % np.real(sigma[m]), "%5.2f" % np.imag(sigma[m]),
  print "| %5.2f" % magSig[m]        , "%5.2f" % relPse[m]        

#------------------------------------------------------
# Define a grid of points on the lambda*Dt plane
# then compute |sigma| at these points 
# so we can plot contours of |sigma| < 1 on the
# lambda*Dt plane.
#------------------------------------------------------
prr = np.linspace(-3.0,0.2,50)   # Real range
pir = np.linspace(-3.0,3.0,50)   # Imagainary range
prc,pic = np.meshgrid(prr,pir);  # Grid of points
pldt=prc + 1j*pic;               # lambda dt values for each point
pSigma = lamSig(pldt);           # sigma at each point
pMagSig = np.abs(pSigma);        # |sigma at each point|

#------------------------------------------------------
# Plot values on the ldt plane and performance vs beta
#------------------------------------------------------
fig = plt.figure(figsize=(15,10))

# Plot the lambda-delta t plane
ax1 = fig.add_subplot(121)
ax1.set_title("$|\sigma|$ and $\lambda_m$ on the $\lambda-\Delta t$ plane")
ax1.set_xlabel("$Re(\lambda\Delta t)$")
ax1.set_ylabel("$Im(\lambda\Delta t)$")
ax1.grid()
lev = np.linspace(0.0,1.0,11)
a = ax1.contour(prc, pic, pMagSig, lev)
fig.colorbar(a, ax=ax1)
ax1.plot(np.real(ldt),np.imag(ldt),'ro',markersize=8)

# Plot the amplification factor and relative phase.
ax2 = fig.add_subplot(122)
ax2.set_title("Performance vs wave number")
ax2.set_xlabel("$\\beta$")
ax2.set_xlim((0, 1))
ax2.grid()
ax2.plot(beta, magSig, '-bo', label="$|\sigma|$")
ax2.plot(beta, relPse, '-ro', label="relative phase")
ax2.legend(loc="lower left")

plt.show()

