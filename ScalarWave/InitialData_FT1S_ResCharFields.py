# Generating C code for plane wave initial
#  data for the scalar wave equation in
#  ***Cartesian*** coordinates, in up to
#  *three* spatial dimensions
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Thiago Assumpcao
#          assumpcaothiago **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-Scalarwave.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import sympy as sp                            # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm                # NRPy+: Reference metric support
import indexedexp as ixp                      # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import sys                                    # Standard Python module for multiplatform OS-level functions
from ScalarWave.CommonParams import wavespeed # NRPy+: Common parameters for all ScalarWave modules (defines wavespeed)

# The name of this module ("InitialData") is given by __name__:
thismodule = __name__

# Set up spherically-symmetric Gaussian initial data
def SphericalGaussian(CoordSystem="Cartesian",default_time=0,default_sigma=3,default_center=1/4):
    # Step 1: Set parameters for the wave
    DIM = par.parval_from_str("grid::DIM")

    # Step 2: Set up Cartesian coordinates in terms of the native CoordSystem we have chosen.
    #         E.g., if CoordSystem="Cartesian", then xx_to_Cart = [xx[0],xx[1],xx[2]]
    #         or if CoordSystem="Spherical", then xx_to_Cart = [xx[0]*sp.sin(xx[1])*sp.cos(xx[2]),
    #                                                       xx[0]*sp.sin(xx[1])*sp.sin(xx[2]),
    #                                                       xx[0]*sp.cos(xx[1])]
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric() # Must call this function to specify rfm.xx_to_Cart
    xx_to_Cart = rfm.xx_to_Cart

    # Step 3: Declare free parameters intrinsic to these initial data
    time  = par.Cparameters("REAL", thismodule, "time",  default_time)
    sigma = par.Cparameters("REAL", thismodule, "sigma", default_sigma)
    center = par.Cparameters("REAL", thismodule, "center", default_center)

    # Step 4: Compute r
    r = sp.sympify(0)
    for i in range(DIM):
        r += xx_to_Cart[i]**2
    r = sp.sqrt(r)
    costheta = xx_to_Cart[2]/r
    sintheta = sp.sqrt( 1-costheta**2 )
    tanphi = xx_to_Cart[1]/xx_to_Cart[0]
    cosphi = xx_to_Cart[0]/sp.sqrt( xx_to_Cart[0]**2 + xx_to_Cart[1]**2 )
    secphi = sp.sympify(1)/cosphi
    sinphi = tanphi*cosphi

    # Following copied from Tutorial-BSSN-basis_transforms
    # Step 2.a: Construct Jacobian & Inverse Jacobians:
    Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    # Step 5: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    global gg_ID, gp_ID, alpha_ID, chi_ID, chip_ID, chipp_ID, omega_ID, betaU, psi_ID, radius_ID, sinth_ID, costh_ID, rsquared_ID, rcubed_ID, rsgg_ID, rsgp_ID, rspsi_ID #, rspsip_ID
    Kcmc = -sp.sympify(1) # parameter for the hyperboloidal slices # make into a parameter as wavespeed
    radius_ID = r
    rsquared_ID = r*r
    rcubed_ID = r*r*r
    sinth_ID = sintheta
    costh_ID = costheta

    rsgg_ID = sp.sympify(0)
    rsgp_ID = sp.sympify(0)
    rspsi_ID = ixp.zerorank1()
    #rspsip_ID = ixp.zerorank2()



    omega_ID = (-Kcmc)*((1)*(1)-r*r)/(6*(1)) # time-independent full-spacetime conformal factor # how can I set the value of params.RMAX here???
    #gg_ID = (+(r-wavespeed*time)/r*sp.exp(-(r-wavespeed*time)**2/(2*sigma**2)) + (r+wavespeed*time)/r*sp.exp(-(r+wavespeed*time)**2/(2*sigma**2))) #+ sp.sympify(1) # original
    gg_ID = sp.sqrt( 1 + r**2 )*1*r*r*sp.exp(-(r*r-center*center)**2/(4*sigma**4)) # spherically symmetric omega_ID*omega_ID*
    #gg_ID = 1*r*r*sp.exp(-(r*r-center*center)**2/(4*sigma**4))*cosphi
    #gg_ID = 1*r*r*omega_ID*omega_ID*sp.exp(-(r*r-center*center*(1+4*sp.sqrt(1-costheta**2)))**2/(4*sigma**4)) # axially symmetric
    #gg_ID = 1*r*r*sp.exp(-(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))**2/sigma**4) # not spherically symmetric *omega_ID*omega_ID
    #vv_ID = sp.sympify(0) #sp.diff(uu_ID, time)
    psi_ID = ixp.zerorank1()
    gp_ID = ( 1 + r**2 )*r*sp.exp(-(r*r-center*center)**2/(4*sigma**4))*( 2-r*r*(r*r-center*center)/sigma**4 ) + r*r*r*sp.exp(-(r*r-center*center)**2/(4*sigma**4))  #spherically symmetric
    #gp_ID = 2*r*sp.exp( -(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))**2/sigma**4 )*( 1-2*r*r*(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))/sigma**4 ) #not spherically symmetric
    #psi_ID[0] = sp.diff(uu_ID, r)
    psi_ID[0] = -sp.sqrt( 1 + r**2 )*r*sp.exp(-(r*r-center*center)**2/(4*sigma**4))*( 2-r*r*(r*r-center*center)/sigma**4 ) #spherically symmetric
    #psi_ID[0] = (r*(center**2*r**2 - r**4 + 2*sigma**4)*cosphi)/(sp.exp((center**2 - r**2)**2/(4.*sigma**4))*sigma**4) #CRAPPY ATTEMPT
    #psi_ID[0] = r*sp.exp(-(r*r-center*center)**2/(4*sigma**4))*( 2-r*r*(r*r-center*center)/sigma**4 )*cosphi #CRAPPY ATTEMPT
    #psi_ID[0] = -2*r*sp.exp( -(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))**2/sigma**4 )*( 1-2*r*r*(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))/sigma**4 ) #not spherically symmetric
    psi_ID[1] = sp.sympify(0)
    psi_ID[2] = sp.sympify(0)
    for i in range(DIM):
        rspsi_ID[i] = r*r*psi_ID[i]
    rsgg_ID = r*r*gg_ID
    rsgp_ID = r*r*gp_ID
    #psi_ID[2] = -1*r*r*sp.exp(-(r*r-center*center)**2/(4*sigma**4))*sinphi #CRAPPY ATTEMPT
    #psi_ID[2] = -8*r*r*sp.exp( -(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))**2/sigma**4 )*( (r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))*center**2*(1+tanphi**2)**(-3/2)*tanphi*secphi**2/sigma**4 )
    #psi_ID[2] = -8*center**2*sp.exp(-(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))**2/sigma**4)*r**2*(1+tanphi**2)*tanphi*(r*r-center*center*(1+4*1/sp.sqrt(1+tanphi**2)))/(sigma**4*(1+tanphi**2)**(-3/2))
    # generalize to the possibility of having cartesian or any other coords
    alpha_ID = sp.sqrt(Kcmc*Kcmc*r*r/9 + omega_ID*omega_ID) # lapse
    #chi_ID = sp.sympify(1) # spatial conformal factor
    chi_ID = sp.sqrt( 1 + r**2 )    #sp.sympify(1)
    chip_ID = r/sp.sqrt( 1 + r**2 )
    chipp_ID = 1/sp.sqrt( 1 + r**2 )**3
    betaSphU = ixp.zerorank1()
    betaU = ixp.zerorank1()
    betaSphU[0] = Kcmc*r/3
    betaSphU[1] = sp.sympify(0)
    betaSphU[2] = sp.sympify(0)
    if CoordSystem=="Spherical":
        for i in range(DIM):
            betaU[i] = betaSphU[i] # I may want to generalise this for coordinate changes in a more elegant way ...
    elif CoordSystem=="Cartesian":
        for i in range(DIM):
            betaU[i] = Kcmc*xx_to_Cart[i]/3
#    else:

# Set up monochromatic plane-wave initial data
def PlaneWave(CoordSystem="Cartesian",default_time=0,default_k0=1,default_k1=1,default_k2=1):
    # Step 1: Set parameters defined in other modules
    DIM = par.parval_from_str("grid::DIM")

    # Step 2: Set up Cartesian coordinates in terms of the native CoordSystem we have chosen.
    #         E.g., if CoordSystem="Cartesian", then xx_to_Cart = [xx[0],xx[1],xx[2]]
    #         or if CoordSystem="Spherical", then xx_to_Cart = [xx[0]*sp.sin(xx[1])*sp.cos(xx[2]),
    #                                                       xx[0]*sp.sin(xx[1])*sp.sin(xx[2]),
    #                                                       xx[0]*sp.cos(xx[1])]
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric()
    xx_to_Cart = rfm.xx_to_Cart

    # Step 3: Declare free parameters intrinsic to these initial data
    time = par.Cparameters("REAL", thismodule, "time", default_time)
    kk   = par.Cparameters("REAL", thismodule, ["kk0", "kk1", "kk2"], [default_k0,default_k1,default_k2])

    # Step 4: Normalize the k vector
    kk_norm_factor = sp.sqrt(kk[0] ** 2 + kk[1] ** 2 + kk[2] ** 2)

    # Step 5: Compute k_norm.x
    dot_product = sp.sympify(0)
    for i in range(DIM):
        dot_product += kk[i] * xx_to_Cart[i]
    dot_product /= kk_norm_factor

    # Step 6: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    global uu_ID, vv_ID, alpha_ID, chi_ID, omega_ID, betaU, psi_ID
    #betaU = ixp.zerorank1()
    uu_ID = sp.sin(dot_product - wavespeed * time) + 2
    vv_ID = sp.diff(uu_ID, time)
    alpha_ID = sp.sympify(1) #dummy
    chi_ID = sp.sympify(1) #dummy
    omega_ID = sp.sympify(1) #dummy
    betaU[0] = sp.sympify(0) #dummy
    betaU[1] = sp.sympify(0) #dummy
    betaU[2] = sp.sympify(0) #dummy
    psi_ID[0] = sp.simpify(0) #dummy
    psi_ID[1] = sp.simpify(0) #dummy
    psi_ID[2] = sp.simpify(0) #dummy

# Initial data driver routine
def InitialData(Type="PlaneWave",CoordSystem="Cartesian",
                default_time=0,
                default_k0=1,default_k1=1,default_k2=1,
                default_sigma=3,default_center=1/4):
    if Type=="PlaneWave":
        PlaneWave(CoordSystem=CoordSystem, default_time=default_time,
                  default_k0=default_k0,default_k1=default_k1,default_k2=default_k2)
    elif Type=="SphericalGaussian":
        SphericalGaussian(CoordSystem=CoordSystem,default_time=default_time,default_sigma=default_sigma,default_center=default_center)
    else:
        print("Error: ScalarWave initial data Type="+str(Type)+" not supported.")
        sys.exit(1)
