# Generating C code for the right-hand-side
#  of the scalar wave equation, in
#  ***curvilinear*** coordinates, using a
#  reference-metric formalism
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-ScalarwaveCurvilinear.ipynb

# Recall from the documentation, the scalar wave
#   equation in curvilinear coordinates can be
#   written as
# \hat{g}^{mu nu} \partial_{mu} \partial_{nu} u - \hat{\Gamma}^{tau} \partial_{tau} u ,
# where \hat{\Gamma}^{tau} is the *contracted* Christoffel symbol.
#
# For reference metrics supported by NRPy+,
#   \hat{g}_{t nu} = -\delta_{t nu},
#   where \delta_{t nu} is the Kronecker delta, we have
# \partial_t^2 u = \hat{g}^{i j} \partial_{i} \partial_{j} u - \hat{\Gamma}^i \partial_i u ,
#   where Latin indices imply the spatial components *only*.
#
# We break up the above equation into two first-order-in-time PDEs,
#   defining v = \partial_t u:
#
# \partial_t u = v
# \partial_t v = \hat{g}^{i j} \partial_{i} \partial_{j} u - \hat{\Gamma}^i \partial_i u

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import sympy as sp                            # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri                            # NRPy+: Functionality for handling numerical grids
import indexedexp as ixp                      # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm                # NRPy+: Reference metric support
#from ScalarWave.CommonParams import wavespeed # NRPy+: Common parameters for all ScalarWave modules (defines wavespeed)

# The name of this module ("ScalarWaveCurvilinear") is given by __name__:
thismodule = __name__

def ScalarWaveCurvilinear_RHSs_MetricComponents():
    # Step 1: Get the spatial dimension, defined in the
    #         NRPy+ "grid" module. With reference metrics,
    #         this must be set to 3 or fewer.
    DIM = par.parval_from_str("DIM")

    # Step 2: Set up the reference metric and
    #         quantities derived from the
    #         reference metric.
    rfm.reference_metric()

    # Step 3: Compute the contracted Christoffel symbols:
    contractedGammahatU = ixp.zerorank1()
    for k in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                contractedGammahatU[k] += rfm.ghatUU[i][j] * rfm.GammahatUDD[k][i][j]

    # INES : my gridfunctions
    global cp_dD, cpx_dD, cpvD_dD #cpxb_dD, cpth_dD, cpph_dD
    global cm_dD, cmx_dD, cmxb_dD, cmth_dD, cmph_dD
    global cthp_dD, cthpx_dD, cthpxb_dD, cthpth_dD, cthpph_dD
    global cphp_dD, cphpx_dD, cphpxb_dD, cphpth_dD, cphpph_dD
    global cthm_dD, cthmx_dD, cthmxb_dD, cthmth_dD, cthmph_dD
    global cphm_dD, cphmx_dD, cphmxb_dD, cphmth_dD, cphmph_dD
    global phi_dD, phix_dD, phixb_dD, phith_dD, phiph_dD
    global oorc_dD, oorcx_dD, oorxbc_dD, oorcth_dD, oorcph_dD
    global hp_dD, hpx_dD, hpxb_dD, hpth_dD, hpph_dD
    global hc_dD, hcx_dD, hcxb_dD, hcth_dD, hcph_dD
    global oochi, oochip, oochipp
    global Radius, ooRp, rcoord
    global sinth, costh
    global sinph, cosph

    # INES : when working with Christian before
    global gg_dD, gp_dD, gp, chi, chip, chipp, Rprime #, Radius, sinth, costh, Rprime, rcoord, psiD, psiD_dD,

    # Step 4: Register gridfunctions that are needed as input
    #         to the scalar wave RHS expressions.


    cp, cpx  = gri.register_gridfunctions("EVOL",["cp", "cpx"]) # cpxb, cpth, cpph ... "cpxb", "cpth", "cpph"
    cpvD = ixp.register_gridfunctions_for_single_rank1("EVOL", "cpvD") #cpvD[0] stands for cpxb
    cm, cmx, cmxb, cmth, cmph = gri.register_gridfunctions("EVOL",["cm", "cmx", "cmxb", "cmth", "cmph"])
    cthp, cthpx, cthpxb, cthpth, cthpph = gri.register_gridfunctions("EVOL",["cthp", "cthpx", "cthpxb", "cthpth", "cthpph"])
    cphp, cphpx, cphpxb, cphpth, cphpph = gri.register_gridfunctions("EVOL",["cphp", "cphpx", "cphpxb", "cphpth", "cphpph"])
    cthm, cthmx, cthmxb, cthmth, cthmph = gri.register_gridfunctions("EVOL",["cthm", "cthmx", "cthmxb", "cthmth", "cthmph"])
    cphm, cphmx, cphmxb, cphmth, cphmph = gri.register_gridfunctions("EVOL",["cphm", "cphmx", "cphmxb", "cphmth", "cphmph"])
    phi, phix, phixb, phith, phiph = gri.register_gridfunctions("EVOL",["phi", "phix", "phixb", "phith", "phiph"])
    oorc, oorcx, oorcxb, oorcth, oorcph = gri.register_gridfunctions("EVOL",["oorc", "oorcx", "oorcxb", "oorcth", "oorcph"])
    hp, hpx, hpxb, hpth, hpph = gri.register_gridfunctions("EVOL",["hp", "hpx", "hpxb", "hpth", "hpph"])
    hc, hcx, hcxb, hcth, hcph = gri.register_gridfunctions("EVOL",["hc", "hcx", "hcxb", "hcth", "hcph"])

    #psi, psix, psixb, psith, psiph = gri.register_gridfunctions("EVOL",["psi", "psix", "psixb", "psith", "psiph"])
    rcoord = gri.register_gridfunctions("EVOL",["rcoord"])
    Radius = gri.register_gridfunctions("EVOL",["Radius"])
    sinth, costh = gri.register_gridfunctions("EVOL",["sinth","costh"])
    sinph, cosph = gri.register_gridfunctions("EVOL",["sinph","cosph"])
    ooRp = gri.register_gridfunctions("EVOL",["ooRp"])
    oochi = gri.register_gridfunctions("EVOL",["oochi"])
    oochip = gri.register_gridfunctions("EVOL",["oochip"])
    oochipp = gri.register_gridfunctions("EVOL",["oochipp"])
    #chip = gri.register_gridfunctions("EVOL",["chip"]) ?

    # INES : gfs that i used when working with Christian
    #gg, gp = gri.register_gridfunctions("EVOL",["gg","gp"])
    #chi = gri.register_gridfunctions("EVOL",["chi"])
    #chip, chipp = gri.register_gridfunctions("EVOL",["chip","chipp"])
    #psiD = ixp.register_gridfunctions_for_single_rank1("EVOL", "psiD") #psiD[0] stands for gminus
    #Rprime = gri.register_gridfunctions("EVOL",["Rprime"])

    #vetU = ixp.register_gridfunctions_for_single_rank1("EVOL", "vetU") #rescaled shift
    #betaU = ixp.zerorank1() # shift

    #global cD
    #cD = ixp.register_gridfunctions_for_single_rank1("EVOL", "cD")
    #for i in range(DIM):
        #betaU[i] = vetU[i]*rfm.ReU[i]

    # Step 5a: Declare the rank-1 indexed expression \partial_{i} u,
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.

    cp_dD = ixp.declarerank1("cp_dD")
    cpx_dD = ixp.declarerank1("cpx_dD")
    cpvD_dD = ixp.declarerank2("cpvD_dD","nosym") # Adding this quantity
    #cpxb_dD = ixp.declarerank1("cpxb_dD")
    #cpth_dD = ixp.declarerank1("cpth_dD")
    #cpph_dD = ixp.declarerank1("cpph_dD")

    cm_dD = ixp.declarerank1("cm_dD")
    cmx_dD = ixp.declarerank1("cmx_dD")
    cmxb_dD = ixp.declarerank1("cmxb_dD")
    cmth_dD = ixp.declarerank1("cmth_dD")
    cmph_dD = ixp.declarerank1("cmph_dD")

    cthp_dD = ixp.declarerank1("cthp_dD")
    cthpx_dD = ixp.declarerank1("cthpx_dD")
    cthpxb_dD = ixp.declarerank1("cthpxb_dD")
    cthpth_dD = ixp.declarerank1("cthpth_dD")
    cthpph_dD = ixp.declarerank1("cthpph_dD")

    cphp_dD = ixp.declarerank1("cphp_dD")
    cphpx_dD = ixp.declarerank1("cphpx_dD")
    cphpxb_dD = ixp.declarerank1("cphpxb_dD")
    cphpth_dD = ixp.declarerank1("cphpth_dD")
    cphpph_dD = ixp.declarerank1("cphpph_dD")

    cthm_dD = ixp.declarerank1("cthm_dD")
    cthmx_dD = ixp.declarerank1("cthmx_dD")
    cthmxb_dD = ixp.declarerank1("cthmxb_dD")
    cthmth_dD = ixp.declarerank1("cthmth_dD")
    cthmph_dD = ixp.declarerank1("cthmph_dD")

    cphm_dD = ixp.declarerank1("cphm_dD")
    cphmx_dD = ixp.declarerank1("cphmx_dD")
    cphmxb_dD = ixp.declarerank1("cphmxb_dD")
    cphmth_dD = ixp.declarerank1("cphmth_dD")
    cphmph_dD = ixp.declarerank1("cphmph_dD")

    phi_dD = ixp.declarerank1("phi_dD")
    phix_dD = ixp.declarerank1("phix_dD")
    phixb_dD = ixp.declarerank1("phixb_dD")
    phith_dD = ixp.declarerank1("phith_dD")
    phiph_dD = ixp.declarerank1("phiph_dD")

    oorc_dD = ixp.declarerank1("oorc_dD")
    oorcx_dD = ixp.declarerank1("oorcx_dD")
    oorcxb_dD = ixp.declarerank1("oorcxb_dD")
    oorcth_dD = ixp.declarerank1("oorcth_dD")
    oorcph_dD = ixp.declarerank1("oorcph_dD")

    hp_dD = ixp.declarerank1("hp_dD")
    hpx_dD = ixp.declarerank1("hpx_dD")
    hpxb_dD = ixp.declarerank1("hpxb_dD")
    hpth_dD = ixp.declarerank1("hpth_dD")
    hpph_dD = ixp.declarerank1("hpph_dD")

    hc_dD = ixp.declarerank1("hc_dD")
    hcx_dD = ixp.declarerank1("hcx_dD")
    hcxb_dD = ixp.declarerank1("hcxb_dD")
    hcth_dD = ixp.declarerank1("hcth_dD")
    hcph_dD = ixp.declarerank1("hcph_dD")


    # INES : derivatives of gfs that I used with Christian
    #gg_dD = ixp.declarerank1("gg_dD")
    #gp_dD = ixp.declarerank1("gp_dD")
    #omega_dD = ixp.declarerank1("omega_dD")
    #alpha_dD = ixp.declarerank1("alpha_dD")
    #chi_dD = ixp.declarerank1("chi_dD")
    #vetU_dD = ixp.declarerank2("vetU_dD","nosym")
    #vetU_dupD = ixp.declarerank2("vetU_dupD","nosym") # Needed for upwinded \beta^i_{,j}
    #betaU_dD = ixp.zerorank2()
    #betaU_dupD = ixp.zerorank2() # Needed for, e.g., \beta^i RHS
    #psiD_dD = ixp.declarerank2("psiD_dD","nosym") #Adding this quantity

    #for i in range(DIM):
        #for j in range(DIM):
            #betaU_dD[i][j] = vetU_dD[i][j]*rfm.ReU[i] + vetU[i]*rfm.ReUdD[i][j]

    #for i in range(DIM):
    #    cD[i] = uu_dD[i]-psiD[i]

    # Step 5b: Declare the rank-2 indexed expression \partial_{ij} u,
    #          which is symmetric about interchange of indices i and j
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.


    # INES : we had this before, but in FT1S it won't be needed
    #gg_dDD = ixp.declarerank2("gg_dDD","sym01")
    #omega_dDD = ixp.declarerank2("omega_dDD","sym01")

    # Step 7: Specify RHSs as global variables,
    #         to enable access outside this
    #         function (e.g., for C code output)

    global cp_rhs, cpx_rhs, cpv_rhsD #cpxb_rhs, cpth_rhs, cpph_rhs
    global cm_rhs, cmx_rhs, cmxb_rhs, cmth_rhs, cmph_rhs
    global cthp_rhs, cthpx_rhs, cthpxb_rhs, cthpth_rhs, cthpph_rhs
    global cphp_rhs, cphpx_rhs, cphpxb_rhs, cphpth_rhs, cphpph_rhs
    global cthm_rhs, cthmx_rhs, cthmxb_rhs, cthmth_rhs, cthmph_rhs
    global cphm_rhs, cphmx_rhs, cphmxb_rhs, cphmth_rhs, cphmph_rhs
    global phi_rhs, phix_rhs, phixb_rhs, phith_rhs, phiph_rhs
    global oorc_rhs, oorcx_rhs, oorcxb_rhs, oorcth_rhs, oorcph_rhs
    global hp_rhs, hpx_rhs, hpxb_rhs, hpth_rhs, hpph_rhs
    global hc_rhs, hcx_rhs, hcxb_rhs, hcth_rhs, hcph_rhs
    global oochi_rhs, oochip_rhs, oochipp_rhs
    global Radius_rhs, ooRp_rhs, rcoord_rhs
    global sinth_rhs, costh_rhs
    global sinph_rhs, cosph_rhs

    #global psi_rhs, psix_rhs, psixb_rhs, psith_rhs, psiph_rhs, Radius_rhs, sinth_rhs, costh_rhs, ooRp_rhs, oochi_rhs, rcoord_rhs

    # INES : what we had before
    #global gg_rhs,gp_rhs,chip_rhs,chipp_rhs #chi_rhs,psi_rhsD, Radius_rhs #sinth_rhs, costh_rhs, Rprime_rhs, rcoord_rhs #, alpha_rhs, omega_rhs, vet_rhsU
    #vet_rhsU    = ixp.zerorank1()
    #beta_rhsU    = ixp.zerorank1()
    #psi_rhsD    = ixp.zerorank1() #Adding this quantity

    # Step 6: Define right-hand sides for the evolution.
    # Step 6a: uu_rhs = vv:

    rcoord_rhs = sp.sympify(0)
    Radius_rhs = sp.sympify(0)
    ooRp_rhs = sp.sympify(0)
    sinth_rhs = sp.sympify(0)
    costh_rhs = sp.sympify(0)
    sinph_rhs = sp.sympify(0)
    cosph_rhs = sp.sympify(0)
    #alpha_rhs = sp.sympify(0)
    oochi_rhs = sp.sympify(0)
    oochip_rhs = sp.sympify(0)
    oochipp_rhs = sp.sympify(0)
    cpv_rhsD    = ixp.zerorank1()


    cp_rhs = (cpvD[0] + cpx*oochi + cp*oochip/oochi)/sp.sympify(2)
    cpx_rhs = ( cpvD[0] + cpvD[0]/(Radius*oochi) + oochi*cpx - cpx/Radius - (costh/sinth)*cpvD[1]/(Radius**2*oochi) + cp*oochip/oochi - cp*oochip/(Radius*oochi**2) + cpvD[0]*oochip/(oochi**2) - ooRp*cpvD[0]*oochip/(2*oochi**2) - cpx*oochip/oochi - ooRp*cpx*oochip/(2*oochi) + cp*oochip**2/(oochi**3) - ooRp*cp*oochip**2/(2*oochi**3) - cp*oochipp/(oochi**2) - cpvD_dD[2][2]/(sinth**2*oochi*Radius**2) - cpvD_dD[1][1]/(oochi*Radius**2) - ooRp*oochip*cp_dD[0]/(oochi**2) - ooRp*cpx_dD[0] )/(ooRp-2) # missing the source term
    cpv_rhsD[0] = - oochi*cpvD[0]/ooRp - cpvD[0]/(ooRp*Radius) - oochi**2*cpx/ooRp + oochi*cpx/(ooRp*Radius) + (costh/sinth)*cpvD[1]/(ooRp*Radius**2) - cp*oochip/ooRp + cp*oochip/(oochi*ooRp*Radius) - cpvD[0]*oochip/(oochi*ooRp) + cpvD_dD[2][2]/(Radius**2*ooRp*sinth**2) + cpvD_dD[1][1]/(Radius**2*ooRp) - cpvD_dD[0][0] # missing the source term # cpxb_rhs
    cpv_rhsD[1] = oochip*cp_dD[1]/(2*oochi) + cpvD_dD[0][1]/2 + oochi*cpx_dD[1]/2 # cpth_rhs
    cpv_rhsD[2] = oochip*cp_dD[2]/(2*oochi) + cpvD_dD[0][2]/2 + oochi*cpx_dD[2]/2 # cpph_rhs

    cm_rhs = sp.sympify(0)
    cmx_rhs = sp.sympify(0)
    cmxb_rhs = sp.sympify(0)
    cmth_rhs = sp.sympify(0)
    cmph_rhs = sp.sympify(0)

    cthp_rhs = sp.sympify(0)
    cthpx_rhs = sp.sympify(0)
    cthpxb_rhs = sp.sympify(0)
    cthpth_rhs = sp.sympify(0)
    cthpph_rhs = sp.sympify(0)

    cphp_rhs = sp.sympify(0)
    cphpx_rhs = sp.sympify(0)
    cphpxb_rhs = sp.sympify(0)
    cphpth_rhs = sp.sympify(0)
    cphpph_rhs = sp.sympify(0)

    cthm_rhs = sp.sympify(0)
    cthmx_rhs = sp.sympify(0)
    cthmxb_rhs = sp.sympify(0)
    cthmth_rhs = sp.sympify(0)
    cthmph_rhs = sp.sympify(0)

    cphm_rhs = sp.sympify(0)
    cphmx_rhs = sp.sympify(0)
    cphmxb_rhs = sp.sympify(0)
    cphmth_rhs = sp.sympify(0)
    cphmph_rhs = sp.sympify(0)

    phi_rhs = sp.sympify(0)
    phix_rhs = sp.sympify(0)
    phixb_rhs = sp.sympify(0)
    phith_rhs = sp.sympify(0)
    phiph_rhs = sp.sympify(0)

    oorc_rhs = sp.sympify(0)
    oorcx_rhs = sp.sympify(0)
    oorcxb_rhs = sp.sympify(0)
    oorcth_rhs = sp.sympify(0)
    oorcph_rhs = sp.sympify(0)

    hp_rhs = sp.sympify(0)
    hpx_rhs = sp.sympify(0)
    hpxb_rhs = sp.sympify(0)
    hpth_rhs = sp.sympify(0)
    hpph_rhs = sp.sympify(0)

    hc_rhs = sp.sympify(0)
    hcx_rhs = sp.sympify(0)
    hcxb_rhs = sp.sympify(0)
    hcth_rhs = sp.sympify(0)
    hcph_rhs = sp.sympify(0)


    #omega_rhs = sp.sympify(0)
    #beta_rhsU[0] = sp.sympify(0)
    #for i in range(DIM):
        #vet_rhsU[i]    =   sp.sympify(0)

    # INES : My generic expressions but not yet
    #psi_rhs = (psixb + psix*oochi - psi*chip*oochi)/sp.sympify(2)
    #psix_rhs = ( sp.sympify(1)/(2/ooRp-1) )* (  -chip*psix*oochi/sp.sympify(2) -chipp*psi/ooRp -chip*psix*oochi/ooRp + psi*chip**2*oochi/sp.sympify(2) + psi*oochi*chip**2/ooRp -chip*psi*( 2*rcoord/(1-rcoord**2) ) + psix*( 2*rcoord/(1-rcoord**2) ) -(chip*psixb)/2 + chip*psixb/ooRp -( 2*rcoord/(1-rcoord**2) )*psixb/oochi + costh*psith/(ooRp*oochi*sinth*Radius**2) -chip*psi_dD[0] + psix_dD[0] + psith_dD[1]/(oochi*ooRp*Radius**2) + psiph_dD[2][2]/(oochi*ooRp*sinth**2*Radius**2)  ) +( -chip*psi +psix -psixb/oochi  )*( sp.sympify(1)/(2/ooRp-1) )/(Radius*ooRp)

    #gg_rhs = sp.sympify(0)
    #gp_rhs = sp.sympify(0)

    # FOR WRITING RHSs EXPRESSIONS RECALL GAMMA^R=-2/R , GAMMA^theta=-COT(THETA)/R^2

    # INES : previous expressions (from when working with Christian)
    #gg_rhs = ( -chip*gg/chi + gp/chi + psiD[0] )/sp.sympify(2)

    # IN THE FOLLOWING WE ALREADY USE SOME CHRISTOFFELS ARE ZERO IN SPHERICAL COORDINATES
    # EXPRESSIONS WHERE METRIC AND CHRISTOFFELS HAVE BEEN SUBSTITUTED
    #THE FOLLOWING TWO EXPRESSIONS ARE SUPPOSED TO BE THE SAME. THE FIRST ONE IS 'SIMPLER'
    #gp_rhs =( sp.sympify(1)/(2*Rprime-1) )* (  -chip*gp/(2*chi) -chipp*gg*Rprime -chip*gp*Rprime/chi + gg*chip**2/(2*chi) + gg*Rprime*chip**2/chi -chip*gg*( 2*rcoord/(1-rcoord**2) ) + gp*( 2*rcoord/(1-rcoord**2) ) -(chip*psiD[0])/2 + chip*Rprime*psiD[0] -chi*( 2*rcoord/(1-rcoord**2) )*psiD[0] + chi*costh*Rprime*psiD[1]/(sinth*Radius**2) -chip*gg_dD[0] + gp_dD[0] + chi*Rprime*psiD_dD[1][1]/(Radius**2) + chi*Rprime*psiD_dD[2][2]/(sinth**2*Radius**2)  )    +( -chip*gg +gp -chi*psiD[0]  )*( sp.sympify(1)/(2*Rprime-1) )*Rprime/Radius

    #gp_rhs = (chip**2*gg)/(sp.sympify(2)*chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*gp)/(sp.sympify(2)*chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chip**2*gg*Rprime)/(chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chipp*gg*Rprime)/(sp.sympify(-1)+sp.sympify(2)*Rprime) - (chip*gp*Rprime)/( chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*gg*Rprime)/(Radius*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (gp*Rprime)/(Radius*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*psiD[0])/( sp.sympify(2)*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chip*Rprime*psiD[0])/(sp.sympify(-1) + sp.sympify(2)*Rprime) - (chi*Rprime*psiD[0])/(Radius*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chi*Rprime*(costh/sinth)*psiD[1])/(Radius**2*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*gg_dD[0])/(sp.sympify(-1) + sp.sympify(2)*Rprime) + gp_dD[0]/(sp.sympify(-1) + sp.sympify(2)*Rprime) + (chi*Rprime*psiD_dD[1][1])/(Radius**2*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chi*Rprime*(sp.sympify(1)/sinth**2)*psiD_dD[2][2])/(Radius**2*(sp.sympify(-1) + sp.sympify(2)*Rprime))

    #gp_rhs = -sp.Rational(3/2)*chip*gp/chi +sp.Rational(3/2)*chip**2*gg/chi +chip*psiD[0]/sp.sympify(2) +gp_dD[0] -chip*gg_dD[0] -chipp*gg  +(sp.sympify(1)/Radius)*(-chip*gg -chi*psiD[0])  +gp/Radius  +chi*psiD_dD[1][1]/Radius**2  +chi*psiD_dD[2][2]/(Radius**2*sinth**2) +chi*psiD[1]*costh/(sinth*Radius**2)
    #gp_rhs = -sp.Rational(3/2)*chip*gp/chi +sp.Rational(3/2)*chip**2*gg/chi +chip*psiD[0]/sp.sympify(2) +gp_dD[0] -chip*gg_dD[0] -chipp*gg +chi*psiD_dD[1][1]/radius**2 +chi*psiD_dD[2][2]/(radius**2*sinth**2) +(sp.sympify(2)/radius)*(gp -chip*gg)/sp.sympify(2) +chi*psiD[1]*costh/(sinth*rsquared) -sp.Rational(1/2)*chi*(sp.sympify(3)*rspsiD_dD[0][0]/rcubed_dD[0] - psiD_dD[0][0] ) #

    #psi_rhsD[0] = -chip*gg*( 2*rcoord/(1-rcoord**2) )/chi + gp*( 2*rcoord/(1-rcoord**2) )/chi + chip*Rprime*psiD[0]/chi - (2*rcoord/(1-rcoord**2))*psiD[0] + costh*Rprime*psiD[1]/(sinth*Radius**2) -psiD_dD[0][0] + Rprime*psiD_dD[1][1]/Radius**2 + Rprime*psiD_dD[2][2]/(sinth**2*Radius**2)    + ( -chip*gg/chi +gp/chi -psiD[0] )*Rprime/Radius

    #psi_rhsD[0] = -((chip*gg*Rprime)/(chi*Radius)) + (gp*Rprime)/(chi*Radius) + ( chip*Rprime*psiD[0])/chi - (Rprime*psiD[0])/Radius + Rprime*(costh/sinth)*psiD[1]/Radius**2 - psiD_dD[0][0] + (Rprime*psiD_dD[1][1])/Radius**2 + (Rprime*(sp.sympify(1)/sinth**2)*psiD_dD[2][2])/Radius**2

    #psi_rhsD[0] = chip*psiD[0]/chi -psiD_dD[0][0]  +(sp.sympify(1)/Radius)*(gp/chi -chip*gg/chi -psiD[0]) +psiD[1]*costh/(sinth*Radius**2) +psiD_dD[2][2]/(Radius**2*sinth**2) +psiD_dD[1][1]/Radius**2

    # THE IDEA OF THE EVANS TRICK IS TO GET RID OF THE POWERS OF R IN THE DENOMINATORS BY RESCALING SOME GRIDFUNCTIONS
    # LOOKING AT THE PREVIOUS EXPRESSIONS THE NEEDED GFs TO BE RESCALED ARE psiD_dD[i][j], gg, gp, psiD[i]


    # ORIGINAL EXPRESSIONS
    #gp_rhs = gp_dD[0] -contractedGammahatU[0]*(gp-psiD[0])/sp.sympify(2) + rfm.ghatUU[1][1]*psiD_dD[1][1] + rfm.ghatUU[2][2]*psiD_dD[2][2]  -contractedGammahatU[1]*psiD[1]
    #psi_rhsD[0] = -psiD_dD[0][0] -contractedGammahatU[0]*(gp-psiD[0])/sp.sympify(2) + rfm.ghatUU[1][1]*psiD_dD[1][1] + rfm.ghatUU[2][2]*psiD_dD[2][2]  -contractedGammahatU[1]*psiD[1]

    #psi_rhsD[1] = -chip*gg_dD[1]/(sp.sympify(2)*chi) + gp_dD[1]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psiD_dD[0][1]
    #psi_rhsD[2] = -chip*gg_dD[2]/(sp.sympify(2)*chi) + gp_dD[2]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psiD_dD[0][2]
    #psi_rhsD[1] = ( gp_dD[1]/chi -chip*gg_dD[1]/chi + psiD_dD[0][1] )/sp.sympify(2)
    #psi_rhsD[2] = ( gp_dD[2]/chi -chip*gg_dD[2]/chi + psiD_dD[0][2] )/sp.sympify(2)

    #rspsi_rhsD[0] = rsquared*chip*psiD[0]/chi -rsquared*psiD_dD[0][0] +psiD_dD[1][1] +psiD_dD[2][2]/sinth**2 +radius*(gp/chi -chip*gg/chi -psiD[0]) +psiD[1]*costh/sinth
    #rspsi_rhsD[1] = rsquared*( gp_dD[1]/chi -chip*gg_dD[1]/chi + psiD_dD[0][1] )/sp.sympify(2)
    #rspsi_rhsD[2] = rsquared*( gp_dD[2]/chi -chip*gg_dD[2]/chi + psiD_dD[0][2] )/sp.sympify(2)
    #for i in range(DIM):
        #rspsi_rhsD[i] = rsquared*psi_rhsD[i]




    # Step 7: Generate C code for scalarwave evolution equations,
    #         print output to the screen (standard out, or stdout).
    # fin.FD_outputC("stdout",
    #                [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
    #                 lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)])

    # Add Kreiss-Oliger dissipation to the RHSs: (taken from nrpytutorial/Tutorial-BaikalETK.ipynb)
        #thismodule = "KO_Dissipation"
    #diss_strength = sp.Rational(2/100)
    #par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

    #gg_dKOD = ixp.declarerank1("gg_dKOD")
    #gp_dKOD = ixp.declarerank1("gp_dKOD")
    #psiD_dKOD = ixp.declarerank2("psiD_dKOD","nosym")

    #for i in range(DIM): #(DIM):
        #gg_rhs += diss_strength*gg_dKOD[i]*rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        #gp_rhs += diss_strength*gp_dKOD[i]*rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        #for j in range(DIM):
            #psi_rhsD[j] += diss_strength*psiD_dKOD[j][i]*rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]"""

