import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import time as tm
plt.close("all")
def ndgrid(x_v,y_v):
    x2, y2 = np.meshgrid(x_v,y_v, indexing="ij")
    return (x2,y2)

###############################################################################
# Physics
n      = 3;                  # stress exponent for power law rheology
Vbc    = 66.4437;            # boundary velocity difference
Lx     = 0.86038;            # domain size x
Ly     = Lx;                 # domain size y
T0     = 49.3269/n;          # initial temperature
r      = 0.0737;             # radius of initial perturbation
Tamp   = 0.1*T0;             # amplitude of initial perturbation
# Numerics
nx     = 94;                 # number of cells x
ny     = nx;                 # number of cells y
nt     = 54;                 # number time steps
nout   = 100;                # check residual each nout iteration
noutp  = 2;                  # display graphics each nout time step
niter  = 1e5;                # max nonlinear iterations
epsi   = 1e-5;               # non-linear tolerance
tetp   = 0.5;                # reduction of PT steps for pressure
tetv   = 0.5;                # reduction of PT steps for velocity
tetT   = 0.5;                # reduction of physical time step for temperature
rel    = 1e-1;               # relaxation factor for non-linear viscosity
Vdamp  = 4.0;                # velocity damping for momentum equations
eta_b  = 1.0;                # numerical compressibility
# Pre-processing
dampx  = 1*(1-Vdamp/nx);     # velocity damping for x-momentum equation
dampy  = 1*(1-Vdamp/ny);     # velocity damping for y-momentum equation
mpow   = -(1-1/n)/2;         # exponent for strain rate dependent viscosity
dx     = Lx/nx;              # grid step in x
dy     = Ly/ny;              # grid step in y
time   = 0;
# Mesh
xc      = np.linspace(-0*Lx/2+dx/2,Lx-dx/2,nx);  # coordonnées centres  
yc      = np.linspace(-0*Ly/2+dy/2,Ly-dy/2,ny);  # coordonnées centres
xv      = np.linspace(-0*Lx/2     ,Lx   ,nx+1);  # coordonnées vertices
yv      = np.linspace(-0*Ly/2     ,Ly   ,ny+1);  # coordonnées vertice    
xc2,yc2   = ndgrid(xc,yc);      # maillage 2D - centres
xv2,yv2   = ndgrid(xv,yv);      # maillage 2D - centres
xvx2,yvx2 = ndgrid(xv,yc);      # maillage 2D - centres
xvy2,yvy2 = ndgrid(xc,yv);      # maillage 2D - centres
# Intial fields
Vx         =  Vbc*xvx2/Lx;
Vy         = -Vbc*yvy2/Ly;
T          =  np.zeros((nx  ,ny  ));
P          =  np.zeros((nx  ,ny  ));
etac       =   np.ones((nx  ,ny  ));
qx         =  np.zeros((nx+1,ny  ));
qy         =  np.zeros((nx  ,ny+1));
dVxdtauVx  =  np.zeros((nx-1,ny  ));
dVydtauVy  =  np.zeros((nx  ,ny-1));
dVxdtauVx0 =  np.zeros((nx-1,ny  ));
dVydtauVy0 =  np.zeros((nx  ,ny-1));
Vx_exp     = np.zeros((nx+1,ny+2));
Vy_exp     = np.zeros((nx+2,ny+1));
T[(xc2**2+yc2**2)<r**2] = Tamp;                                            # initial temperature pertubation
dtT        = tetT*1/4.1*min(dx,dy)**2;                                      # explicit timestep for 2D diffusion
E = 0; W = 0;
for it in range(nt): # ------ Physical timesteps
    print(it)
    To    = T;                                                             # temperature from previous step (for backward-Euler integration)
    time  = time + dtT;                                                    # update physical time
    iter = 0; resid = epsi*2; resX=1; resY=1; resPt=1; tic=tm.time();
    while  resid > epsi and iter < niter: # ------ Pseudo-Transient cycles
        dVxdtauVx0      = dVxdtauVx + dampx*dVxdtauVx0;                   # used for damping x momentum residuals
        dVydtauVy0      = dVydtauVy + dampy*dVydtauVy0;                   # used for damping y momentum residuals
        # ------ Kinematics
        Vx_exp[:,1:-1] = Vx[:,:]; Vx_exp[:,0:1] = Vx[:,0:1]; Vx_exp[:,-1:] = Vx[:,-1:]
        Vy_exp[1:-1,:] = Vy[:,:]; Vy_exp[0:1,:] = Vy[0:1,:]; Vy_exp[-1:,:] = Vy[-1:,:]
        divV           = np.diff(Vx,n=1,axis=0)/dx + np.diff(Vy,n=1,axis=1)/dy;
        Exxc           = np.diff(Vx,n=1,axis=0)/dx - 1/2*divV;
        Eyyc           = np.diff(Vy,n=1,axis=1)/dy - 1/2*divV;
        Exyv           = 0.5*(np.diff(Vx_exp,n=1,axis=1)/dy + np.diff(Vy_exp,n=1,axis=0)/dx);
        Exyc           = 0.25*(Exyv[0:-1,0:-1] + Exyv[0:-1,1:] + Exyv[1:,0:-1] + Exyv[1:,1:]);
        Eii2           = 0.5*(Exxc**2 + Eyyc**2) + Exyc**2;                      # strain rate invariant
        # ------ Rheology
        etac_phys       = Eii2**mpow*np.exp( -T*(1./(1+T/T0)) );                 # physical viscosity
        etac            = np.exp(rel*np.log(etac_phys) + (1-rel)*np.log(etac));  # numerical shear viscosity
        etav            = np.zeros((nx+1,ny+1));                                 # expand viscosity fom cell centroids to vertices
        etav[1:-1,1:-1] = 0.25*(etac[0:-1,0:-1] + etac[0:-1,1:] + etac[1:,0:-1] + etac[1:,1:]);
        etav[:,0:1] = etav[:,1:2];
        etav[:,-1:] = etav[:,-2:-1];
        etav[0:1,:] = etav[1:2,:];
        etav[-1:,:] = etav[-2:-1,:];
        # ------ Pseudo-Time steps
        dtauP           = tetp*  4.1/min(nx,ny)*etac*(1.0+eta_b);
        dtauVx          = tetv*1/4.1*(min(dx,dy)**2/( 0.5*(etac[1:,:] + etac[0:-1,:]) ))/(1+eta_b);
        dtauVy          = tetv*1/4.1*(min(dx,dy)**2/( 0.5*(etac[:,1:] + etac[:,0:-1]) ))/(1+eta_b);
        dtauT           = tetT*1/4.1*min(dx,dy)**2;
        # ------ Fluxes
        qx[1:-1,:]      = -np.diff(T,n=1,axis=0)/dx;
        qy[:,1:-1]      = -np.diff(T,n=1,axis=1)/dy;
        Sxx             = -P + 2*etac*(Exxc + eta_b*divV);
        Syy             = -P + 2*etac*(Eyyc + eta_b*divV);
        Txy             =      2*etav*Exyv;
        Hs              = 4*etac*Eii2;
        # ------ Residuals
        dVxdtauVx       =               np.diff(Txy[1:-1,:],n=1,axis=1)/dy + np.diff(Sxx,n=1,axis=0)/dx;
        dVydtauVy       =               np.diff(Txy[:,1:-1],n=1,axis=0)/dx + np.diff(Syy,n=1,axis=1)/dy;
        dPdtauP         =            -  divV;
        dTdtauT         = (To-T)/dtT - (np.diff(qx,n=1,axis=0)/dx + np.diff(qy,n=1,axis=1)/dy) + Hs;
        # ------ Updates
        Vx[1:-1,:]      = Vx[1:-1,:] + dtauVx*(dVxdtauVx + dampx*dVxdtauVx0); # update with damping
        Vy[:,1:-1]      = Vy[:,1:-1] + dtauVy*(dVydtauVy + dampy*dVydtauVy0); # update with damping
        P               = P          + dtauP *dPdtauP;
        T               = T          + dtauT *dTdtauT;
        iter  += 1;
        if np.mod(iter,nout)==0:
             fp = dPdtauP.flatten(); fT = dTdtauT.flatten(); fu = np.hstack((dVxdtauVx.flatten(), dVydtauVy.flatten()))
             ru = la.norm(fu)/len(fu);
             rp = la.norm(fp)/len(fp);
             rT = la.norm(fT)/len(fT);
             resid = max(max( ru, rp), rT)
             print( ru)
# Plot
Wtime = tm.time()-tic;
#print(iter)
#print(Wtime)
#Vxc = 0.5*(Vx[0:-1,:] + Vx[1:,:])
#Vyc = 0.5*(Vy[:,0:-1] + Vy[:,1:])
#V   = np.sqrt( Vxc**2 + Vyc**2 );
#plt.figure(1); plt.clf()
#plt.subplot(221);plt.contourf(xc2, yc2,  Eii2, 256, cmap=plt.cm.jet); plt.axis('image'); plt.colorbar()
#plt.subplot(222);plt.contourf(xc2, yc2,     P, 256, cmap=plt.cm.jet); plt.axis('image'); plt.colorbar()
#plt.subplot(223);plt.contourf(xc2, yc2,     T, 256, cmap=plt.cm.jet); plt.axis('image'); plt.colorbar()
#plt.subplot(224);plt.contourf(xc2, yc2,     V, 256, cmap=plt.cm.jet); plt.axis('image'); plt.colorbar()
#plt.show()
