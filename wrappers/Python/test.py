import sys
sys.path.append('../../src')

import InternalFlow
print dir(InternalFlow)
CPS = InternalFlow.CoolPropStateClassSI("Propane")
CPS.update(InternalFlow.iT,300,InternalFlow.iQ,0.5)
print InternalFlow.Lockhart_Martinelli_1949_DPDZ_f(CPS,100,0.1,0.5)

sys.path.append('../../../achp_git/src')
import ACHP
CPS2 = ACHP.CoolPropStateClassSI("Propane")
CPS2.update(ACHP.iT,300,ACHP.iQ,0.5)
print InternalFlow.Lockhart_Martinelli_1949_DPDZ_f(CPS2,100,0.1,0.5)