import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri                # NRPy+: Functions having to do with numerical grids
#import ScalarWave.ScalarWaveCurvilinear_RHSs_MetricComponents_FT1S_CharacteristicFields  as swrhs
import ScalarWave.ScalarWaveCurvilinear_RHSs_MetricComponents_FT1S_GBUF  as swrhs

thismodule = __name__

def Constraint():

    DIM = 3

# Alex: this should register the quantity - it may have to be modified / completed
    global CDg, CDb, CDu, CDf
    #for i in range(len(gri.glb_gridfcs_list)):
            #if "CDg0" in gri.glb_gridfcs_list[i].name:
                #return
            #if "CDb0" in gri.glb_gridfcs_list[i].name:
                #return
            #if "CDu0" in gri.glb_gridfcs_list[i].name:
                #return
            #if "CDf0" in gri.glb_gridfcs_list[i].name:
                #return

    CDg = ixp.register_gridfunctions_for_single_rank1("AUX", "CDG")
    CDg = ixp.zerorank1()
    CDb = ixp.register_gridfunctions_for_single_rank1("AUX", "CDB")
    CDb = ixp.zerorank1()
    CDu = ixp.register_gridfunctions_for_single_rank1("AUX", "CDU")
    CDu = ixp.zerorank1()
    CDf = ixp.register_gridfunctions_for_single_rank1("AUX", "CDF")
    CDf = ixp.zerorank1()


# Alex: here we give the expression for the constraint

    CDg[0] = -swrhs.gp/swrhs.chi +( sp.sympify(1)/(sp.sympify(2)*swrhs.Rprime-sp.sympify(1)) )*(swrhs.psigD[0] -swrhs.chip*swrhs.gg/swrhs.chi +sp.sympify(2)*swrhs.gg_dD[0] )
    CDg[1] = swrhs.psigD[1]-swrhs.gg_dD[1]
    CDg[2] = swrhs.sinth*swrhs.psigD[2]-swrhs.gg_dD[2]

    CDb[0] = -swrhs.bp/swrhs.chi +( sp.sympify(1)/(sp.sympify(2)*swrhs.Rprime-sp.sympify(1)) )*(swrhs.psibD[0] -swrhs.chip*swrhs.bb/swrhs.chi +sp.sympify(2)*swrhs.bb_dD[0] )
    CDb[1] = swrhs.psibD[1]-swrhs.bb_dD[1]
    CDb[2] = swrhs.sinth*swrhs.psibD[2]-swrhs.bb_dD[2]

    CDu[0] = -swrhs.up/swrhs.chi +( sp.sympify(1)/(sp.sympify(2)*swrhs.Rprime-sp.sympify(1)) )*(swrhs.psiuD[0] -swrhs.chip*swrhs.uu/swrhs.chi +sp.sympify(2)*swrhs.uu_dD[0] )
    CDu[1] = swrhs.psiuD[1]-swrhs.uu_dD[1]
    CDu[2] = swrhs.sinth*swrhs.psiuD[2]-swrhs.uu_dD[2]

    CDf[0] = -swrhs.fgaugep/swrhs.chi +( sp.sympify(1)/(sp.sympify(2)*swrhs.Rprime-sp.sympify(1)) )*(swrhs.psifD[0] -swrhs.chip*swrhs.fgauge/swrhs.chi +sp.sympify(2)*swrhs.fgauge_dD[0] )
    CDf[1] = swrhs.psifD[1]-swrhs.fgauge_dD[1]
    CDf[2] = swrhs.sinth*swrhs.psifD[2]-swrhs.fgauge_dD[2]

