#ifndef THERMALCORR_BRAZED_PLATE_HEAT_EXCHANGER_H
#define THERMALCORR_BRAZED_PLATE_HEAT_EXCHANGER_H

#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"

namespace ThermalCorr
{
namespace BrazedPlateHX
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

	class BPHEGeometry
	{
	public:
		BPHEGeometry(){PlateAmplitude = _HUGE; 
					   PlateWavelength = _HUGE; 
					   InclinationAngle = _HUGE; 
					   PlateThickness = _HUGE;
					   Bp = _HUGE; 
					   Lp = _HUGE; };
					   
		double PlateAmplitude; ///< The amplitude of the corrugations of the plates [m]
		double PlateWavelength; ///< The wavelength of the plate corrugations [m]
		double PlateThickness; ///< The thickness of the plates [m]
		double InclinationAngle; ///< The inclination angle of the plates [rad]
		double Bp; ///< The width between the centers of the ports [m]
		double Lp; ///< The length between the centers of the ports [m]
	};
	struct BPHEData
	{
		double Re, ///< Reynolds number [-]
			   HTC, ///< Heat transfer coefficient [W/m^2/K]
			   mdot_per_channel, ///< Mass flow per channel [kg/s]
			   DELTAP; ///< Pressure change [Pa]
		shared_ptr<CoolProp::AbstractState> AS;
	};
	void BPHE_1phase(BPHEGeometry, BPHEData &);
};
};
#endif