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

    global gg_dD, gp_dD, bb_dD, bp_dD, uu_dD, up_dD, fgauge_dD, fgaugep_dD, psigD, psigD_dD, psibD, psibD_dD, psiuD, psiuD_dD, psifD, psifD_dD, gg, gp, bb, bp, uu, up, fgauge, fgaugep, chi, chip, chipp, Radius, sinth, costh, Rprime, rcoord
    # Step 4: Register gridfunctions that are needed as input
    #         to the scalar wave RHS expressions.
    gg, gp = gri.register_gridfunctions("EVOL",["gg","gp"])
    bb, bp = gri.register_gridfunctions("EVOL",["bb","bp"])
    uu, up = gri.register_gridfunctions("EVOL",["uu","up"])
    fgauge, fgaugep = gri.register_gridfunctions("EVOL",["fgauge","fgaugep"])
    chi = gri.register_gridfunctions("EVOL",["chi"])
    chip, chipp = gri.register_gridfunctions("EVOL",["chip","chipp"])

    psigD = ixp.register_gridfunctions_for_single_rank1("EVOL", "psigD") #psiD[0] stands for gminus
    psibD = ixp.register_gridfunctions_for_single_rank1("EVOL", "psibD") #psiD[0] stands for gminus
    psiuD = ixp.register_gridfunctions_for_single_rank1("EVOL", "psiuD") #psiD[0] stands for gminus
    psifD = ixp.register_gridfunctions_for_single_rank1("EVOL", "psifD") #psiD[0] stands for gminus
    rcoord = gri.register_gridfunctions("EVOL",["rcoord"])
    Radius = gri.register_gridfunctions("EVOL",["Radius"])
    Rprime = gri.register_gridfunctions("EVOL",["Rprime"])
    sinth, costh = gri.register_gridfunctions("EVOL",["sinth","costh"])

    # THE FOLLOWING ARE THE GFs REQUIRED FOR EVANS TRICK (comming from psiD_dD[i][j], gg, gp, psiD[i].. look at other comment)
    #rsgg = gri.register_gridfunctions("EVOL",["rsgg"])
    #rsgp = gri.register_gridfunctions("EVOL",["rsgp"])
    #rspsiD = ixp.register_gridfunctions_for_single_rank1("EVOL", "rspsiD")
    #rspsipDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "rspsipDD", "nosym")

    #global cD
    #cD = ixp.register_gridfunctions_for_single_rank1("EVOL", "cD")
    #for i in range(DIM):
        #betaU[i] = vetU[i]*rfm.ReU[i]

    # Step 5a: Declare the rank-1 indexed expression \partial_{i} u,
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.

    gg_dD = ixp.declarerank1("gg_dD")
    gp_dD = ixp.declarerank1("gp_dD")
    bb_dD = ixp.declarerank1("bb_dD")
    bp_dD = ixp.declarerank1("bp_dD")
    uu_dD = ixp.declarerank1("uu_dD")
    up_dD = ixp.declarerank1("up_dD")
    fgauge_dD = ixp.declarerank1("fgauge_dD")
    fgaugep_dD = ixp.declarerank1("fgaugep_dD")

    chi_dD = ixp.declarerank1("chi_dD")

    psigD_dD = ixp.declarerank2("psigD_dD","nosym") #Adding this quantity
    psibD_dD = ixp.declarerank2("psibD_dD","nosym") #Adding this quantity
    psiuD_dD = ixp.declarerank2("psiuD_dD","nosym") #Adding this quantity
    psifD_dD = ixp.declarerank2("psifD_dD","nosym") #Adding this quantity


    #EVANS TRICK REQUIRES DERIVATIVES OF NEWLY DEFINED RESCALED GFs. HERE THEY GO
    #rsgg_dD = ixp.declarerank1("rsgg_dD")
    #rsgp_dD = ixp.declarerank1("rsgp_dD")
    #rspsiD_dD = ixp.declarerank2("rspsiD_dD","nosym")
    #rspsipDD_dD = ixp.declarerank3("rspsipDD_dD","nosym")

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
    gg_dDD = ixp.declarerank2("gg_dDD","sym01")
    bb_dDD = ixp.declarerank2("bb_dDD","sym01")
    uu_dDD = ixp.declarerank2("uu_dDD","sym01")
    fgauge_dDD = ixp.declarerank2("fgauge_dDD","sym01")
    #omega_dDD = ixp.declarerank2("omega_dDD","sym01")

    # Step 7: Specify RHSs as global variables,
    #         to enable access outside this
    #         function (e.g., for C code output)
    global gg_rhs, gp_rhs, chi_rhs, chip_rhs, chipp_rhs, psig_rhsD, psib_rhsD, Radius_rhs, sinth_rhs, costh_rhs, Rprime_rhs, rcoord_rhs, bb_rhs, bp_rhs, uu_rhs, up_rhs, psiu_rhsD, fgauge_rhs, fgaugep_rhs, psif_rhsD

    psig_rhsD = ixp.zerorank1() #Adding this quantity
    psib_rhsD = ixp.zerorank1() #Adding this quantity
    psiu_rhsD = ixp.zerorank1() #Adding this quantity
    psif_rhsD = ixp.zerorank1() #Adding this quantity

    # Step 6: Define right-hand sides for the evolution.
    # Step 6a: uu_rhs = vv:
    rcoord_rhs = sp.sympify(0)
    Radius_rhs = sp.sympify(0)
    Rprime_rhs = sp.sympify(0)
    sinth_rhs = sp.sympify(0)
    costh_rhs = sp.sympify(0)
    #alpha_rhs = sp.sympify(0)
    chi_rhs = sp.sympify(0)
    chip_rhs = sp.sympify(0)
    chipp_rhs = sp.sympify(0)


    #RHSs FOR EVANS TRICK GFs
    #rsgg_rhs = sp.sympify(0)
    #rsgp_rhs = sp.sympify(0)
    #rspsi_rhsD = ixp.zerorank1()
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

    gg_rhs = ( -chip*gg/chi + gp/chi + psigD[0] )/sp.sympify(2)

    # IN THE FOLLOWING WE ALREADY USE SOME CHRISTOFFELS ARE ZERO IN SPHERICAL COORDINATES
    # EXPRESSIONS WHERE METRIC AND CHRISTOFFELS HAVE BEEN SUBSTITUTED
    #THE FOLLOWING TWO EXPRESSIONS ARE SUPPOSED TO BE THE SAME. THE FIRST ONE IS 'SIMPLER'
    gp_rhs =( sp.sympify(1)/(2*Rprime-1) )* ( -chip*gg*(2*rcoord/(1-rcoord**2)) + gp*(2*rcoord/(1-rcoord**2)) -chi*(2*rcoord/(1-rcoord**2))*psigD[0]      -chip*gp/(2*chi) -chipp*gg*Rprime -chip*gp*Rprime/chi + gg*chip**2/(2*chi) + gg*Rprime*chip**2/chi -(chip*psigD[0])/2 + chip*Rprime*psigD[0] + chi*costh*Rprime*psigD[1]/(sinth*Radius**2) -chip*gg_dD[0] + gp_dD[0] + chi*Rprime*psigD_dD[1][1]/(Radius**2) + chi*Rprime*psigD_dD[2][2]/(sinth*Radius**2)  )     # """substitute Evans:"""  +( -chip*gg +gp -chi*psigD[0]  )*( sp.sympify(1)/(2*Rprime-1) )*Rprime/Radius

    #gp_rhs = (chip**2*gg)/(sp.sympify(2)*chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*gp)/(sp.sympify(2)*chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chip**2*gg*Rprime)/(chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chipp*gg*Rprime)/(sp.sympify(-1)+sp.sympify(2)*Rprime) - (chip*gp*Rprime)/( chi*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*gg*Rprime)/(Radius*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (gp*Rprime)/(Radius*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*psiD[0])/( sp.sympify(2)*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chip*Rprime*psiD[0])/(sp.sympify(-1) + sp.sympify(2)*Rprime) - (chi*Rprime*psiD[0])/(Radius*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chi*Rprime*(costh/sinth)*psiD[1])/(Radius**2*(sp.sympify(-1) + sp.sympify(2)*Rprime)) - (chip*gg_dD[0])/(sp.sympify(-1) + sp.sympify(2)*Rprime) + gp_dD[0]/(sp.sympify(-1) + sp.sympify(2)*Rprime) + (chi*Rprime*psiD_dD[1][1])/(Radius**2*(sp.sympify(-1) + sp.sympify(2)*Rprime)) + (chi*Rprime*(sp.sympify(1)/sinth**2)*psiD_dD[2][2])/(Radius**2*(sp.sympify(-1) + sp.sympify(2)*Rprime))

    #gp_rhs = -sp.Rational(3/2)*chip*gp/chi +sp.Rational(3/2)*chip**2*gg/chi +chip*psiD[0]/sp.sympify(2) +gp_dD[0] -chip*gg_dD[0] -chipp*gg  +(sp.sympify(1)/Radius)*(-chip*gg -chi*psiD[0])  +gp/Radius  +chi*psiD_dD[1][1]/Radius**2  +chi*psiD_dD[2][2]/(Radius**2*sinth**2) +chi*psiD[1]*costh/(sinth*Radius**2)
    #gp_rhs = -sp.Rational(3/2)*chip*gp/chi +sp.Rational(3/2)*chip**2*gg/chi +chip*psiD[0]/sp.sympify(2) +gp_dD[0] -chip*gg_dD[0] -chipp*gg +chi*psiD_dD[1][1]/radius**2 +chi*psiD_dD[2][2]/(radius**2*sinth**2) +(sp.sympify(2)/radius)*(gp -chip*gg)/sp.sympify(2) +chi*psiD[1]*costh/(sinth*rsquared) -sp.Rational(1/2)*chi*(sp.sympify(3)*rspsiD_dD[0][0]/rcubed_dD[0] - psiD_dD[0][0] ) #

    psig_rhsD[0] = -chip*gg*(2*rcoord/(1-rcoord**2))/chi + gp*(2*rcoord/(1-rcoord**2))/chi - (2*rcoord/(1-rcoord**2))*psigD[0]    +chip*Rprime*psigD[0]/chi + costh*Rprime*psigD[1]/(sinth*Radius**2) -psigD_dD[0][0] + Rprime*psigD_dD[1][1]/Radius**2 + Rprime*psigD_dD[2][2]/(sinth*Radius**2)    # """substitute Evans:"""  + ( -chip*gg/chi +gp/chi -psiD[0] )*Rprime/Radius

    #psi_rhsD[0] = -((chip*gg*Rprime)/(chi*Radius)) + (gp*Rprime)/(chi*Radius) + ( chip*Rprime*psiD[0])/chi - (Rprime*psiD[0])/Radius + Rprime*(costh/sinth)*psiD[1]/Radius**2 - psiD_dD[0][0] + (Rprime*psiD_dD[1][1])/Radius**2 + (Rprime*(sp.sympify(1)/sinth**2)*psiD_dD[2][2])/Radius**2

    #psi_rhsD[0] = chip*psiD[0]/chi -psiD_dD[0][0]  +(sp.sympify(1)/Radius)*(gp/chi -chip*gg/chi -psiD[0]) +psiD[1]*costh/(sinth*Radius**2) +psiD_dD[2][2]/(Radius**2*sinth**2) +psiD_dD[1][1]/Radius**2

    # THE IDEA OF THE EVANS TRICK IS TO GET RID OF THE POWERS OF R IN THE DENOMINATORS BY RESCALING SOME GRIDFUNCTIONS
    # LOOKING AT THE PREVIOUS EXPRESSIONS THE NEEDED GFs TO BE RESCALED ARE psiD_dD[i][j], gg, gp, psiD[i]


    # ORIGINAL EXPRESSIONS
    #gp_rhs = gp_dD[0] -contractedGammahatU[0]*(gp-psiD[0])/sp.sympify(2) + rfm.ghatUU[1][1]*psiD_dD[1][1] + rfm.ghatUU[2][2]*psiD_dD[2][2]  -contractedGammahatU[1]*psiD[1]
    #psi_rhsD[0] = -psiD_dD[0][0] -contractedGammahatU[0]*(gp-psiD[0])/sp.sympify(2) + rfm.ghatUU[1][1]*psiD_dD[1][1] + rfm.ghatUU[2][2]*psiD_dD[2][2]  -contractedGammahatU[1]*psiD[1]

    psig_rhsD[1] = -chip*psigD[1]/(sp.sympify(2)*chi) + gp_dD[1]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psigD_dD[0][1]
    psig_rhsD[2] = -chip*psigD[2]/(sp.sympify(2)*chi) + gp_dD[2]/(sp.sympify(2)*chi*sinth) + sp.Rational(1/2)*psigD_dD[0][2]/sinth
    #psig_rhsD[1] = -chip*gg_dD[1]/(sp.sympify(2)*chi) + gp_dD[1]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psigD_dD[0][1]
    #psig_rhsD[2] = -chip*gg_dD[2]/(sp.sympify(2)*chi*sinth) + gp_dD[2]/(sp.sympify(2)*chi*sinth) + sp.Rational(1/2)*psigD_dD[0][2]/sinth
    #psi_rhsD[1] = ( gp_dD[1]/chi -chip*gg_dD[1]/chi + psiD_dD[0][1] )/sp.sympify(2)
    #psi_rhsD[2] = ( gp_dD[2]/chi -chip*gg_dD[2]/chi + psiD_dD[0][2] )/sp.sympify(2)

    #rspsi_rhsD[0] = rsquared*chip*psiD[0]/chi -rsquared*psiD_dD[0][0] +psiD_dD[1][1] +psiD_dD[2][2]/sinth**2 +radius*(gp/chi -chip*gg/chi -psiD[0]) +psiD[1]*costh/sinth
    #rspsi_rhsD[1] = rsquared*( gp_dD[1]/chi -chip*gg_dD[1]/chi + psiD_dD[0][1] )/sp.sympify(2)
    #rspsi_rhsD[2] = rsquared*( gp_dD[2]/chi -chip*gg_dD[2]/chi + psiD_dD[0][2] )/sp.sympify(2)
    #for i in range(DIM):
        #rspsi_rhsD[i] = rsquared*psi_rhsD[i]

    # THE FOLLOWING RHSs ARE COPIED FROM THE GOOD ONES.. STILL MISSING THE OTHER TERMS
    # CHECK RHSs FOR B FIELDS!!!!!
    bb_rhs = ( -chip*bb/chi + bp/chi + psibD[0] )/sp.sympify(2)

    bp_rhs = ( sp.sympify(1)/(2*Rprime-1) )*( -chip*bb*(2*rcoord/(1-rcoord**2)) + bp*(2*rcoord/(1-rcoord**2)) -chi*(2*rcoord/(1-rcoord**2))*psibD[0]          -bp*chip/(2*chi) +bb*chip**2/(2*chi) -bp*chip*Rprime/chi +bb*chip**2*Rprime/chi -bb*chipp*Rprime +chip*fgauge*Rprime/(2*chi) -fgaugep*Rprime/(2*chi) -chip**2*gg**2*Rprime/(4*chi**2) +chip*gg*gp*Rprime/(2*chi**2) -gp**2*Rprime/(4*chi**2) -chip*psibD[0]/2 +chip*Rprime*psibD[0] +chi*Rprime*costh*psibD[1]/(sinth*Radius**2) -Rprime*psifD[0]/2 +chip*gg*Rprime*psigD[0]/(2*chi) -gp*Rprime*psigD[0]/(2*chi) -Rprime*psigD[0]**2/4 -chip*bb_dD[0] +bp_dD[0] +chi*Rprime*psibD_dD[1][1]/(Radius**2) +chi*Rprime*psibD_dD[2][2]/(sinth*Radius**2)    )     # """substitute Evans:"""  +( -chip*bb +bp -chi*psibD[0]  )*( sp.sympify(1)/(2*Rprime-1) )*Rprime/Radius

    #( sp.sympify(1)/(2*Rprime-1) )* (  -chip*bp/(2*chi) -chipp*bb*Rprime -chip*bp*Rprime/chi + bb*chip**2/(2*chi) + bb*Rprime*chip**2/chi -(chip*psibD[0])/2 + chip*Rprime*psibD[0] + chi*costh*Rprime*psibD[1]/(sinth*Radius**2) -chip*bb_dD[0] + bp_dD[0] + chi*Rprime*psibD_dD[1][1]/(Radius**2) + chi*Rprime*psibD_dD[2][2]/(sinth**2*Radius**2)  )    +( -chip*bb +bp -chi*psibD[0]  )*( sp.sympify(1)/(2*Rprime-1) )*Rprime/Radius             +(sp.sympify(1)/(2*Rprime-1) )*(   chip*fgauge*Rprime/(2*chi) -fgaugep*Rprime/(2*chi) -chip**2*gg**2*Rprime/(4*chi**2) +chip*gg*gp*Rprime/(2*chi**2) -gp**2*Rprime/(4*chi**2)        -Rprime*psifD[0]/2 +chip*gg*Rprime*psigD[0]/(2*chi) -gp*Rprime*psigD[0]/(2*chi) -Rprime*psigD[0]**2/4    )
    #-chip*bb*( 2*rcoord/(1-rcoord**2) ) + bp*( 2*rcoord/(1-rcoord**2) ) -chi*( 2*rcoord/(1-rcoord**2) )*psibD[0]

    psib_rhsD[0] = -chip*bb*(2*rcoord/(1-rcoord**2))/chi + bp*(2*rcoord/(1-rcoord**2))/chi - (2*rcoord/(1-rcoord**2))*psibD[0]         +chip*fgauge*Rprime/(2*chi**2) -fgaugep*Rprime/(2*chi**2) -chip**2*gg**2*Rprime/(4*chi**3) +chip*gg*gp*Rprime/(2*chi**3) -gp**2*Rprime/(4*chi**3)  +chip*Rprime*psibD[0]/chi +Rprime*costh*psibD[1]/(sinth*Radius**2) -Rprime*psifD[0]/(2*chi) +chip*gg*Rprime*psigD[0]/(2*chi**2) -gp*Rprime*psigD[0]/(2*chi**2) -Rprime*psigD[0]**2/(4*chi) -psibD_dD[0][0] +Rprime*psibD_dD[1][1]/Radius**2 +Rprime*psibD_dD[2][2]/(sinth*Radius**2)     # """substitute Evans:"""  + ( -chip*bb/chi +bp/chi -psibD[0] )*Rprime/Radius

    #-chip*bb*( 2*rcoord/(1-rcoord**2) )/chi + bp*( 2*rcoord/(1-rcoord**2) )/chi + chip*Rprime*psibD[0]/chi -(2*rcoord/(1-rcoord**2))*psibD[0] + costh*Rprime*psibD[1]/(sinth*Radius**2) -psibD_dD[0][0] + Rprime*psibD_dD[1][1]/Radius**2 + Rprime*psibD_dD[2][2]/(sinth**2*Radius**2)    + ( -chip*bb/chi +bp/chi -psibD[0] )*Rprime/Radius         +chip*fgauge*Rprime/(2*chi**2) -fgaugep*Rprime/(2*chi**2) -chip**2*gg**2*Rprime/(4*chi**3) +chip*gg*gp*Rprime/(2*chi**3) -gp**2*Rprime/(4*chi**3)    -Rprime*psifD[0]/(2*chi) +chip*gg*Rprime*psigD[0]/(2*chi**2) -gp*Rprime*psigD[0]/(2*chi**2) -Rprime*psigD[0]**2/(4*chi)

    psib_rhsD[1] = -chip*psibD[1]/(sp.sympify(2)*chi) + bp_dD[1]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psibD_dD[0][1]

    psib_rhsD[2] = -chip*psibD[2]/(sp.sympify(2)*chi) + bp_dD[2]/(sp.sympify(2)*chi*sinth) + sp.Rational(1/2)*psibD_dD[0][2]/sinth





    uu_rhs = ( -chip*uu/chi + up/chi + psiuD[0] )/sp.sympify(2)

    up_rhs = ( sp.sympify(1)/(2*Rprime-1) )*( -chip*uu*(2*rcoord/(1-rcoord**2)) + up*(2*rcoord/(1-rcoord**2)) -chi*(2*rcoord/(1-rcoord**2))*psiuD[0]       -chip*up/(2*chi) -Rprime*up/chi -chip*Rprime*up/chi +chip**2*uu/(2*chi) +chip*Rprime*uu/chi +chip**2*Rprime*uu/chi -chipp*Rprime*uu -chip*psiuD[0]/2 -Rprime*psiuD[0] +chip*Rprime*psiuD[0] +chi*Rprime*costh*psiuD[1]/(sinth*Radius**2) +up_dD[0] -chip*uu_dD[0] +chi*Rprime*psiuD_dD[1][1]/(Radius**2) +chi*Rprime*psiuD_dD[2][2]/(sinth*Radius**2)       )   # """substitute Evans:"""  +( -chip*uu +up -chi*psiuD[0]  )*( sp.sympify(1)/(2*Rprime-1) )*Rprime/Radius

    psiu_rhsD[0] = -chip*uu*(2*rcoord/(1-rcoord**2))/chi + up*(2*rcoord/(1-rcoord**2))/chi - (2*rcoord/(1-rcoord**2))*psiuD[0]        -Rprime*up/chi**2 +chip*Rprime*uu/chi**2 -Rprime*psiuD[0]/chi +chip*Rprime*psiuD[0]/chi +Rprime*costh*psiuD[1]/(sinth*Radius**2) -psiuD_dD[0][0] +Rprime*psiuD_dD[1][1]/Radius**2 + Rprime*psiuD_dD[2][2]/(sinth*Radius**2)  # """substitute Evans:"""  + ( -chip*uu/chi +up/chi -psiuD[0] )*Rprime/Radius

    psiu_rhsD[1] = -chip*psiuD[1]/(sp.sympify(2)*chi) + up_dD[1]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psiuD_dD[0][1]

    psiu_rhsD[2] = -chip*psiuD[2]/(sp.sympify(2)*chi) + up_dD[2]/(sp.sympify(2)*chi*sinth) + sp.Rational(1/2)*psiuD_dD[0][2]/sinth





    fgauge_rhs = ( -chip*fgauge/chi + fgaugep/chi + psifD[0] )/sp.sympify(2)

    fgaugep_rhs = ( sp.sympify(1)/(2*Rprime-1) )* ( -chip*fgauge*(2*rcoord/(1-rcoord**2)) + fgaugep*(2*rcoord/(1-rcoord**2)) -chi*(2*rcoord/(1-rcoord**2))*psifD[0]          +chip**2*fgauge/(2*chi) -chip*fgaugep/(2*chi) +chip*fgauge*Rprime/chi +chip**2*fgauge*Rprime/chi -chipp*fgauge*Rprime -fgaugep*Rprime/chi -chip*fgaugep*Rprime/chi -chip**2*gg**2*Rprime/(2*chi**2) +chip*gg*gp*Rprime/(chi**2) -gp**2*Rprime/(2*chi**2) -chip*psifD[0]/2 -Rprime*psifD[0] +chip*Rprime*psifD[0] +chi*Rprime*costh*psifD[1]/(sinth*Radius**2) +chip*gg*Rprime*psigD[0]/chi -gp*Rprime*psigD[0]/chi -Rprime*psigD[0]**2/2 -chip*fgauge_dD[0] +fgaugep_dD[0] +chi*Rprime*psifD_dD[1][1]/(Radius**2) +chi*Rprime*psifD_dD[2][2]/(sinth*Radius**2)   )  # """substitute Evans:"""  +( -chip*fgauge +fgaugep -chi*psifD[0]  )*( sp.sympify(1)/(2*Rprime-1) )*Rprime/Radius

    psif_rhsD[0] = -chip*fgauge*(2*rcoord/(1-rcoord**2))/chi +fgaugep*(2*rcoord/(1-rcoord**2))/chi - (2*rcoord/(1-rcoord**2))*psifD[0]      +chip*fgauge*Rprime/chi**2 -fgaugep*Rprime/chi**2 -chip**2*gg**2*Rprime/(2*chi**3) +chip*gg*gp*Rprime/chi**3 -gp**2*Rprime/(2*chi**3) -Rprime*psifD[0]/chi +chip*Rprime*psifD[0]/chi +Rprime*costh*psifD[1]/(sinth*Radius**2) +chip*gg*Rprime*psigD[0]/chi**2 -gp*Rprime*psigD[0]/chi**2 -Rprime*psigD[0]**2/(2*chi) -psifD_dD[0][0] +Rprime*psifD_dD[1][1]/Radius**2 +Rprime*psifD_dD[2][2]/(sinth*Radius**2)   # """substitute Evans:"""  + ( -chip*fgauge/chi +fgaugep/chi -psifD[0] )*Rprime/Radius

    psif_rhsD[1] = -chip*psifD[1]/(sp.sympify(2)*chi) + fgaugep_dD[1]/(sp.sympify(2)*chi) + sp.Rational(1/2)*psifD_dD[0][1]

    psif_rhsD[2] = -chip*psifD[2]/(sp.sympify(2)*chi) + fgaugep_dD[2]/(sp.sympify(2)*chi*sinth) + sp.Rational(1/2)*psifD_dD[0][2]/sinth


    # Step 7: Generate C code for scalarwave evolution equations,
    #         print output to the screen (standard out, or stdout).
    # fin.FD_outputC("stdout",
    #                [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
    #                 lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)])

    # Add Kreiss-Oliger dissipation to the RHSs: (taken from nrpytutorial/Tutorial-BaikalETK.ipynb)
        #thismodule = "KO_Dissipation"
    diss_strength = sp.Rational(2/100) #FIX DISSIPATION/
    #par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

    gg_dKOD = ixp.declarerank1("gg_dKOD")
    gp_dKOD = ixp.declarerank1("gp_dKOD")
    psigD_dKOD = ixp.declarerank2("psigD_dKOD","nosym")
    bb_dKOD = ixp.declarerank1("bb_dKOD")
    bp_dKOD = ixp.declarerank1("bp_dKOD")
    psibD_dKOD = ixp.declarerank2("psibD_dKOD","nosym")
    uu_dKOD = ixp.declarerank1("uu_dKOD")
    up_dKOD = ixp.declarerank1("up_dKOD")
    psiuD_dKOD = ixp.declarerank2("psiuD_dKOD","nosym")
    fgauge_dKOD = ixp.declarerank1("fgauge_dKOD")
    fgaugep_dKOD = ixp.declarerank1("fgaugep_dKOD")
    psifD_dKOD = ixp.declarerank2("psifD_dKOD","nosym")

    for i in range(DIM): #(DIM):
        gg_rhs += diss_strength*gg_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        gp_rhs += diss_strength*gp_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        bb_rhs += diss_strength*bb_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        bp_rhs += diss_strength*bp_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        uu_rhs += diss_strength*uu_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        up_rhs += diss_strength*up_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        fgauge_rhs += diss_strength*fgauge_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        fgaugep_rhs += diss_strength*fgaugep_dKOD[i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]
        for j in range(DIM):
            psig_rhsD[j] += diss_strength*psigD_dKOD[j][i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]"""
            psib_rhsD[j] += diss_strength*psibD_dKOD[j][i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]"""
            psiu_rhsD[j] += diss_strength*psiuD_dKOD[j][i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]"""
            psif_rhsD[j] += diss_strength*psifD_dKOD[j][i] # *rfm.ReU[i] # ReU[k] = 1/scalefactor_orthog_funcform[k]"""

