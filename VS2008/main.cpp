
#include "CoolProp.h"
#include "CPState.h"
#include "TwoPhaseHTDP.h"
#include <iostream>
#include "time.h"


int main()
{
	clock_t t1,t2;

	double Kim2 = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("REFPROP-PROPANE", 100, 0.001, 1500000,1.0,0.7);
	double Kim0 = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("Propane", 100, 0.001, 1500000,1.0,0.7);

	double Kim = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("Propane", 100, 0.001, 1500000,1.0,0.7);
	double Kim1 = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("Propane", 100, 0.001, 1500000,1.0,0.9);

	double p = IProps(iP,iT,40+273.15,iQ,0,get_Fluid_index("R134a"))*1000;
	double Kim4 = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("R134a", 100, 0.0014, p, 1.0, 0.75);
	double Cavallini = Cavallini_2009_AnnularMist_DPDZ_f("R134a", 100, 0.0014, p, 0.05);
	double LM = Lockhart_Martinelli_1949_DPDZ_f("R134a", 100, 0.0014, p, 0.75);
	double Friedl = Friedl_1979_DPDZ_f("R134a", 100, 0.0014, p, 0.75);

	double pp = IProps(iP,iT,40+273.15,iQ,0,get_Fluid_index("R410A"))*1000;
	double Kim5 = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("R410A", 50, 0.001, pp, -1, 0.8);

	double LM2 = Lockhart_Martinelli_1949_DPDZ_f("Propane", 100, 0.0015, 9.5e5, 0.8);
	double Kim6 = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("Propane", 100, 0.0015, 9.5e5, -1, 0.8);

	double Zivi = Zivi_DPDZ_a("R134a",100,1500000,0.9999,0.0001);

	t1 = clock();
	long N = 1000*10;
	double rr = 4;
	
	for (double Q = 0; Q <= 1.0; Q+=1/((double)N-1))
	{	
		//CoolPropStateClass CPS = CoolPropStateClass("Propane");
		//CPS.update(iP,1500,iQ,Q+1e-20*Q);

		//IProps(iT,iP,1500,iQ,Q+1e-20*Q,get_Fluid_index("Propane"));

		rr = Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f("Propane", 100, 0.01, 1500000,-1,Q+1e-20*Q);
	}
	t2 = clock();
	printf("elapsed time %g us/call\n",double(t2-t1)/CLOCKS_PER_SEC/N*1e6);
	return 0;
}