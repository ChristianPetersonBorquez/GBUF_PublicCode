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

    global gg_dD, gp_dD, psiD, psiD_dD, gp, chi, chip, chipp, radius, sinth, costh, rsquared, rsquared_dD, rcubed, rcubed_dD, rspsiD, rspsiD_dD
    # Step 4: Register gridfunctions that are needed as input
    #         to the scalar wave RHS expressions.
    gg, gp = gri.register_gridfunctions("EVOL",["gg","gp"])
    omega = gri.register_gridfunctions("EVOL",["omega"])
    alpha, chi = gri.register_gridfunctions("EVOL",["alpha","chi"])
    chip, chipp = gri.register_gridfunctions("EVOL",["chip","chipp"])
    vetU = ixp.register_gridfunctions_for_single_rank1("EVOL", "vetU") #rescaled shift
    betaU = ixp.zerorank1() # shift
    psiD = ixp.register_gridfunctions_for_single_rank1("EVOL", "psiD") #psiD[0] stands for gminus
    radius = gri.register_gridfunctions("EVOL",["radius"])
    rsquared = gri.register_gridfunctions("EVOL",["rsquared"])
    rcubed = gri.register_gridfunctions("EVOL",["rcubed"])
    sinth, costh = gri.register_gridfunctions("EVOL",["sinth","costh"])

    # THE FOLLOWING ARE THE GFs REQUIRED FOR EVANS TRICK (comming from psiD_dD[i][j], gg, gp, psiD[i].. look at other comment)
    rsgg = gri.register_gridfunctions("EVOL",["rsgg"])
    rsgp = gri.register_gridfunctions("EVOL",["rsgp"])
    rspsiD = ixp.register_gridfunctions_for_single_rank1("EVOL", "rspsiD")
    #rspsipDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "rspsipDD", "nosym")

    #global cD
    #cD = ixp.register_gridfunctions_for_single_rank1("EVOL", "cD")
    for i in range(DIM):
        betaU[i] = vetU[i]*rfm.ReU[i]

    # Step 5a: Declare the rank-1 indexed expression \partial_{i} u,
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.

    gg_dD = ixp.declarerank1("gg_dD")
    gp_dD = ixp.declarerank1("gp_dD")
    omega_dD = ixp.declarerank1("omega_dD")
    alpha_dD = ixp.declarerank1("alpha_dD")
    chi_dD = ixp.declarerank1("chi_dD")
    vetU_dD = ixp.declarerank2("vetU_dD","nosym")
    vetU_dupD = ixp.declarerank2("vetU_dupD","nosym") # Needed for upwinded \beta^i_{,j}
    betaU_dD = ixp.zerorank2()
    betaU_dupD = ixp.zerorank2() # Needed for, e.g., \beta^i RHS
    psiD_dD = ixp.declarerank2("psiD_dD","nosym") #Adding this quantity

    rsquared_dD = ixp.declarerank1("rsquared_dD")
    rcubed_dD = ixp.declarerank1("rcubed_dD")

    #EVANS TRICK REQUIRES DERIVATIVES OF NEWLY DEFINED RESCALED GFs. HERE THEY GO
    rsgg_dD = ixp.declarerank1("rsgg_dD")
    rsgp_dD = ixp.declarerank1("rsgp_dD")
    rspsiD_dD = ixp.declarerank2("rspsiD_dD","nosym")
    #rspsipDD_dD = ixp.declarerank3("rspsipDD_dD","nosym")

    for i in range(DIM):
        for j in range(DIM):
            betaU_dD[i][j] = vetU_dD[i][j]*rfm.ReU[i] + vetU[i]*rfm.ReUdD[i][j]

    #for i in range(DIM):
    #    cD[i] = uu_dD[i]-psiD[i]

    # Step 5b: Declare the rank-2 indexed expression \partial_{ij} u,
    #          which is symmetric about interchange of indices i and j
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.
    gg_dDD = ixp.declarerank2("gg_dDD","sym01")
    omega_dDD = ixp.declarerank2("omega_dDD","sym01")

    # Step 7: Specify RHSs as global variables,
    #         to enable access outside this
    #         function (e.g., for C code output)
    global gg_rhs,gp_rhs,alpha_rhs,chi_rhs,chip_rhs,chipp_rhs,omega_rhs,vet_rhsU, psi_rhsD, radius_rhs, sinth_rhs, costh_rhs, rsquared_rhs, rcubed_rhs, rsgg_rhs, rsgp_rhs, rspsi_rhsD, rspsip_rhsDD
    vet_rhsU    = ixp.zerorank1()
    beta_rhsU    = ixp.zerorank1()
    psi_rhsD    = ixp.zerorank1() #Adding this quantity

    # Step 6: Define right-hand sides for the evolution.
    # Step 6a: uu_rhs = vv:
    radius_rhs = sp.sympify(0)
    rsquared_rhs = sp.sympify(0)
    rcubed_rhs = sp.sympify(0)
    sinth_rhs = sp.sympify(0)
    costh_rhs = sp.sympify(0)
    alpha_rhs = sp.sympify(0)
    chi_rhs = sp.sympify(0)
    chip_rhs = sp.sympify(0)
    chipp_rhs = sp.sympify(0)
    omega_rhs = sp.sympify(0)
    beta_rhsU[0] = sp.sympify(0)
    for i in range(DIM):
        vet_rhsU[i]    =   sp.sympify(0)

    #RHSs FOR EVANS TRICK GFs
    rsgg_rhs = sp.sympify(0)
    rsgp_rhs = sp.sympify(0)
    rspsi_rhsD = ixp.zerorank1()
    #rspsip_rhsDD = ixp.zerorank2()


    #Introduce r bc we need it in RHSs .. THIS DOESNT WORK
    #r = sp.sqrt(rfm.ghatDD[1][1])
    #xx_to_Cart = rfm.xx_to_Cart
    #r = sp.sympify(0)
    #for i in range(DIM):
        #r += xx_to_Cart[i]**2
    #r = sp.sqrt(r)
    #costheta = xx_to_Cart[2]/r
    #sintheta = sp.sqrt( xx_to_Cart[0]**2 + xx_to_Cart[1]**2 )/r
    #cottheta = costheta/sintheta
    #csctheta = 1/sintheta

    # FOR WRITING RHSs EXPRESSIONS RECALL GAMMA^R=-2/R , GAMMA^theta=-COT(THETA)/R^2

    gg_rhs = ( gp/chi - chip*gg/chi + psiD[0] )/sp.sympify(2)

    # IN THE FOLLOWING WE ALREADY USE SOME CHRISTOFFELS ARE ZERO IN SPHERICAL COORDINATES
    # EXPRESSIONS WHERE METRIC AND CHRISTOFFELS HAVE BEEN SUBSTITUTED
    gp_rhs = -sp.Rational(3/2)*chip*gp/chi +sp.Rational(3/2)*chip**2*gg/chi +chip*psiD[0]/sp.sympify(2) +gp_dD[0] -chip*gg_dD[0] -chipp*gg +chi*psiD_dD[1][1]/radius**2 +chi*psiD_dD[2][2]/(radius**2*sinth**2) +(sp.sympify(2)/radius)*(gp -chip*gg -chi*psiD[0])/sp.sympify(2) +chi*psiD[1]*costh/(sinth*radius**2)
    #gp_rhs = -sp.Rational(3/2)*chip*gp/chi +sp.Rational(3/2)*chip**2*gg/chi +chip*psiD[0]/sp.sympify(2) +gp_dD[0] -chip*gg_dD[0] -chipp*gg +chi*psiD_dD[1][1]/radius**2 +chi*psiD_dD[2][2]/(radius**2*sinth**2) +(sp.sympify(2)/radius)*(gp -chip*gg)/sp.sympify(2) +chi*psiD[1]*costh/(sinth*rsquared) -sp.Rational(1/2)*chi*(sp.sympify(3)*rspsiD_dD[0][0]/rcubed_dD[0] - psiD_dD[0][0] ) #
    psi_rhsD[0] = chip*psiD[0]/chi -psiD_dD[0][0] +psiD_dD[1][1]/radius**2 +psiD_dD[2][2]/(radius**2*sinth**2) +(sp.sympify(2)/radius)*(gp/chi -chip*gg/chi -psiD[0])/sp.sympify(2) +psiD[1]*costh/(sinth*radius**2)

    # THE IDEA OF THE EVANS TRICK IS TO GET RID OF THE POWERS OF R IN THE DENOMINATORS BY RESCALING SOME GRIDFUNCTIONS
    # LOOKING AT THE PREVIOUS EXPRESSIONS THE NEEDED GFs TO BE RESCALED ARE psiD_dD[i][j], gg, gp, psiD[i]


    # ORIGINAL EXPRESSIONS
    #gp_rhs = gp_dD[0] -contractedGammahatU[0]*(gp-psiD[0])/sp.sympify(2) + rfm.ghatUU[1][1]*psiD_dD[1][1] + rfm.ghatUU[2][2]*psiD_dD[2][2]  -contractedGammahatU[1]*psiD[1]
    #psi_rhsD[0] = -psiD_dD[0][0] -contractedGammahatU[0]*(gp-psiD[0])/sp.sympify(2) + rfm.ghatUU[1][1]*psiD_dD[1][1] + rfm.ghatUU[2][2]*psiD_dD[2][2]  -contractedGammahatU[1]*psiD[1]

    psi_rhsD[1] = ( gp_dD[1]/chi -chip*gg_dD[1]/chi + psiD_dD[0][1] )/sp.sympify(2)
    psi_rhsD[2] = ( gp_dD[2]/chi -chip*gg_dD[2]/chi + psiD_dD[0][2] )/sp.sympify(2)

    for i in range(DIM):
        rspsi_rhsD[i] = rsquared*psi_rhsD[i]


    #vv_rhs = sp.sympify(0)
    #for i in range(DIM):
        #vv_rhs += -contractedGammahatU[i]*psiD[i] #uu_dD[i] # ->
        #for j in range(DIM):
            # PART 1:
            #vv_rhs += + rfm.ghatUU[i][j]*psiD_dD[j][i] #uu_dDD[i][j] # ->
            #for k in range(DIM):
                #vv_rhs +=  0

    #for i in range(DIM): # Also, add constraint
        #psi_rhsD[i]    =   vv_dD[i]

    # Step 6b: The right-hand side of the \partial_t v equation
    #          is given by:
    #          \hat{g}^{ij} \partial_i \partial_j u - \hat{\Gamma}^i \partial_i u.
    #          ^^^^^^^^^^^^ PART 1 ^^^^^^^^^^^^^^^^ ^^^^^^^^^^ PART 2 ^^^^^^^^^^^
#    vv_rhs = 0
#    for i in range(DIM):
#        # PART 2:
#        vv_rhs -= alpha*alpha*chi*contractedGammahatU[i]*uu_dD[i]
#        vv_rhs += beta_rhsU[i]*uu_dD[i] -((vv*betaU[i]*alpha_dD[i])/alpha) - (3*vv*betaU[i]*chi_dD[i])/(2*chi) - (2*vv*betaU[i]*omega_dD[i])/(omega) - (betaU[i]*alpha_rhs*uu_dD[i])/alpha - (3*betaU[i]*chi_rhs*uu_dD[i])/(2*chi) + 2*betaU[i]*vv_dD[i] +vv*betaU_dD[i][i]
#        for j in range(DIM):
#            # PART 1:
#            vv_rhs += alpha*alpha*chi*rfm.ghatUU[i][j]*uu_dDD[i][j] + alpha*chi*rfm.ghatUU[i][j]*uu_dD[i]*alpha_dD[j] - 2*alpha*alpha*chi*rfm.ghatUU[i][j]*uu_dD[i]*omega_dD[j]/omega - alpha*alpha*rfm.ghatUU[i][j]*uu_dD[i]*chi_dD[j]/2 # add all missing beta terms - do it in a systematic way
#            vv_rhs += ((3*betaU[i]*betaU[j]*chi_dD[j]*uu_dD[i])/(2*chi)) + (betaU[i]*betaU[j]*alpha_dD[i]*uu_dD[j])/alpha + (2*betaU[i]*betaU[j]*omega_dD[i]*uu_dD[j])/(omega) - betaU[i]*uu_dD[j]*betaU_dD[j][i] - betaU[i]*uu_dD[i]*betaU_dD[j][j] - betaU[i]*betaU[j]*uu_dDD[i][j]
#            for k in range(DIM):
#                vv_rhs += vv*betaU[k]*rfm.GammahatUDD[i][i][k] - betaU[i]*betaU[k]*uu_dD[j]*rfm.GammahatUDD[j][i][k] - betaU[i]*betaU[k]*uu_dD[i]*rfm.GammahatUDD[j][j][k] + betaU[i]*betaU[j]*uu_dD[k]*rfm.GammahatUDD[k][j][i]
#    vv_rhs += vv*alpha_rhs/alpha + 3*vv*chi_rhs/(2*chi) # add omega_rhs as well? and term for the metric's time derivative?
#     Rescaled uu and vv by omega - there may be a mistake somewhere, so I'm making a copy below to try to find it

    # Step 7: Generate C code for scalarwave evolution equations,
    #         print output to the screen (standard out, or stdout).
    # fin.FD_outputC("stdout",
    #                [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
    #                 lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)])

    # Add Kreiss-Oliger dissipation to the RHSs: (taken from nrpytutorial/Tutorial-BaikalETK.ipynb)
        #thismodule = "KO_Dissipation"
    diss_strength = sp.sympify(0)/sp.sympify(1000)  #FIX DISSIPATION
    #par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

    gg_dKOD = ixp.declarerank1("gg_dKOD")
    gp_dKOD = ixp.declarerank1("gp_dKOD")
    psiD_dKOD = ixp.declarerank2("psiD_dKOD","nosym")

    for i in range(1): #(DIM):
        gg_rhs += radius**2*diss_strength*gg_dKOD[i]*rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        gp_rhs += radius**2*diss_strength*gp_dKOD[i]*rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        for j in range(DIM):
            psi_rhsD[j] += radius**2*diss_strength*psiD_dKOD[j][i]*rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]"""

