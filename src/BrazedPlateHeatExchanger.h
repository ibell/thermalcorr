#ifndef THERMALCORR_BRAZED_PLATE_HEAT_EXCHANGER_H
#define THERMALCORR_BRAZED_PLATE_HEAT_EXCHANGER_H

#include "../externals/coolprop/CoolProp/CPState.h"

namespace BrazedPlateHeatExchanger
{
	/*
	Based on the single-phase pressure drop and heat transfer correlations
	in VDI Heat Atlas Chapter N6: Pressure Drop and Heat Transfer in Plate Heat 
	Exchangers by Holger Martin DOI: 10.1007/978-3-540-77877-6_66 Springer Verlag
	::
		=============================        
		||   __               __    ||
		||  /  \             /  \   ||
		|| |    |           |    |  ||  ===
		||  \__/             \__/   ||   |
		||                          ||   |
		||             | <-  B   -> ||   |
		||                          ||   |
		||                          ||   |
		||                          ||
		||                          ||
		||             |\           ||
		||             | \          ||   Lp
		||             |  \         ||  
		||             |   \        ||
		||             |phi \       ||
		||             |     \      ||   |
		||                          ||   |
		||   __               __    ||   |
		||  /  \             /  \   ||   |
		|| |    |           |    |  ||  ===
		||  \__/             \__/   ||
		||                          ||
		=============================
		 | -----      Bp  --------- |
	     
		 phi is the inclination angle
	*/

	class BPHEInputs
	{
	public:
		double PlateAmplitude, 
			PlateWavelength, 
			PlateConductivity,
			Bp, 
			Lp, 
			Nplates, 
			PlateLength,
			InclinationAngle,
			Re,
			HTC,
			mdot_channel
			;
		int MoreChannels;
		CoolPropStateClassSI CPS;
	};
	void BPHE_1phase(BPHEInputs BPHE);
};

#endif