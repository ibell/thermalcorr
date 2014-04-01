#ifndef INTERNALFLOW_H
#define INTERNALFLOW_H

#include "CPState.h"

namespace ThermalCorr
{
namespace GeneralInternal
{
/// Shah Condensation 
/// from Shah, M., 1976. A New Correlation for Heat Transfer During Boiling Flow Through Pipes. ASHRAE Transactions 82, 66-86.
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param D Diameter [m]
/// @param x Quality [-]
/// @return HTC Heat transfer coefficient [W/m^2/K]
double Shah_1976_HTC(CoolPropStateClassSI *CPS, double G, double D, double x);

/// Homogeneous equilibrium model using the Adams formulation for the two-phase viscosity
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
double HEM_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x);

/// A. Cavallini, D. Del Col, M. Matkovic, L. Rossetto, "Frictional pressure drop during vapour-liquid flow in minichannels: Modelling
/// and experimental evaluation", International Journal of Heat and Fluid Flow 30 (2009) 131–139
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
/// Erratum: log() in equation 8 for E is actually log_10, or base-10 logarithm, not natural logarithm
double Cavallini_2009_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x);

/// Frictional pressure drop of Lockhart-Martinelli in tubes
/// Lockhart, R.W., Martinelli, R.C., 1949, Proposed Correlation of Data for Isothermal Two-Phase 
/// Two-Component Flow in Pipes. Chemical Engineering Progress. v. 45, 39-48
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
double Lockhart_Martinelli_1949_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x);

/// Friedel, L. (1979) Improved friction pressure drop correlations for horizontal and vertical 
/// two-phase pipe flow. European Two-Phase Flow Group Meeting, Ispra, Italy, paper E2
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
double Friedl_1979_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x);

/// Accelerational pressure drop using the slip ratio definition of Zivi
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param x1 Inlet Quality [-]
/// @param x2 Outlet Quality [-]
/// @return HTC Heat transfer coefficient [Pa/m]
/*!
Description
-----------
From the consideration of two-phase flow analysis, the accelerational presssure drop can be obtained.  It is caused by the change in velocity of the vapor and liquid phases due to phase change, which in boiling creates vapor and accelerates the vapor, or in the case of condensation, reduces the vapor velocity, resulting in a pressure increase.
    
\f[
\left( \frac{\partial p}{\partial z}\right)_A=-G^2\frac{d}{dz}\left[\frac{x^2v_g}{\epsilon}+\frac{(1-x)^2v_f}{1-\epsilon}\right]    
\f]

where \f$\epsilon\f$ is the refrigerant vapor void fraction (typically the symbol \f$\alpha\f$ is used for void fraction, but here we are using that for heat transfer coefficient). Integrating over the length where the quality goes from \f$x_1\f$ to \f$x_2\f$ yields

\f[
\Delta p_A=\int_{0}^{L}\left[\left( \frac{\partial p}{\partial z}\right)_A dz\right]
\f]
    
\f[
\Delta p_A=-G^2L\left[\left(\frac{x_2^2v_g}{\epsilon_2}+\frac{(1-x_2)^2v_f}{1-\epsilon_2}\right) -\left(\frac{x_1^2v_g}{\epsilon_1}+\frac{(1-x_1)^2v_f}{1-\epsilon_1} \right) \right]
\f]
        
where \f$\Delta p_A\f$ is negative if the pressure is dropping in going from quality \f$x_1\f$ to \f$x_2\f$.  If the quality in the term 

\f[
\label{eq:bracketedtermDPa}
\left(\frac{x^2v_g}{\epsilon}+\frac{(1-x)^2v_f}{1-\epsilon} \right)
\f]
    
is 0 or 1, one part is zero and the other is an indeterminate form of 0/0.  One evaluation of L'Hopital's rule can be used to show that if the quality is zero, the term in Equation \ref{eq:bracketedtermDPa} is equal to \f$v_f\f$, or if the quality is 1, this term is equal to \f$v_g\f$.
*/
double Zivi_DPDZ_a(CoolPropStateClassSI *CPS, double G, double x1, double x2);
};
};

#endif