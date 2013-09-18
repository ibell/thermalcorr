#ifndef TWOPHASEHTDP_H
#define TWOPHASEHTDP_H

#define EXPORT_CODE
#define CONVENTION

// Notes: The format for the doxygen string must strictly follow the template in order to be converted to Python. Non-matching elements will be dropped
// Things that are supported: inline math, in both reST and doxygen

/// Shah Condensation 
/// from Shah, M., 1976. A New Correlation for Heat Transfer During Boiling Flow Through Pipes. ASHRAE Transactions 82, 66-86.
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param D Diameter [m]
/// @param p Pressure [Pa]
/// @param x Quality [-]
/// @return HTC Heat transfer coefficient [W/m^2/K]
EXPORT_CODE double CONVENTION Shah_Condensation_HTC(char *Fluid, double G, double D, double p, double x);

/// Sung-Min Kim, Issam Mudawar, "Universal approach to predicting two-phase frictional pressure drop for adiabatic
/// and condensing mini/micro-channel flows", International Journal of Heat and Mass Transfer 55 (2012) 3246–3261
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param beta Aspect ratio (W/H) [-] (<0 for non-rectangular channel)
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]

/*!
Description
-----------
Microchannel pressure drop correlation from Kim and Mudawar \footnote{Sung-Min Kim, Issam Mudawar, "Universal approach to predicting two-phase frictional pressure drop for adiabatic and condensing mini/micro-channel flows", International Journal of Heat and Mass Transfer 55 (2012) 3246–-3261}

\f[Re_{fo} = \frac{GD_h}{\mu_f}\f]
\f[Su_{go} = \frac{\rho\sigma D_h}{\mu_g^2}\f]
\f[Re_f = \frac{G(1-x)D_h}{\mu_f}\f]
\f[Re_g = \frac{GxD_h}{\mu_g}\f]

Friction factor (\f$\beta = w_i/h_i\f$, \f$k\f$ is one of \f$f\f$ or \f$g\f$)

Re | Type | Friction factor 
---|------|----------------
\f$\mathrm{Re}_k < 2000\f$ | Rect. | \f$f_k = \frac{ 24(1-1.3553\beta+1.9467\beta^2-1.7012\beta^3+0.9564\beta^4-0.2537\beta^5)}{\mathrm{Re}_k}\f$
\f$\mathrm{Re}_k < 2000\f$ | Non-Rect. | \f$f_k = \frac{16}{\mathrm{Re}_k}\f$
\f$2000 < \mathrm{Re}_k < 20,000\f$ | Any | \f$f_k = \frac{0.079}{\mathrm{Re}_k^{0.25}}\f$
\f$\mathrm{Re}_k > 20,000\f$ | Any | \f$f_k = \frac{0.046}{\mathrm{Re}_k^{0.2}}\f$

While generally the correlations describe the frictional pressure drop as a positive quantity, it is the author's opinion that use of the pressure gradient is more straightforward to understand.

Thus the single-phase frictional pressure gradients are given by
\f[
\left.\frac{dp}{dz}\right|_{f} = -\frac{2 f_f G^2(1-x)^2}{\rho_f D_h}
\f]

\f[
\left.\frac{dp}{dz}\right|_{g} = -\frac{2 f_g G^2 x^2}{\rho_g D_h}
\f]

And the Lockhart-Martinelli parameter is given by
\f[
X^2 = \frac{\left.\frac{dp}{dz}\right|_{f}}{\left.\frac{dp}{dz}\right|_{g}}
\f]
And the Lockhart-Martinelli $C$ factor is given by Table \ref{tab:LMCfactor}

Lockhart-Martinelli \f$C\f$ factor

\f$\mathrm{Re}_f\f$ | \f$\mathrm{Re}_g\f$ | \f$C\f$ factor
----------------|---------------|-------------
\f$Re_f<2000\f$ | \f$Re_g<2000\f$ | \f$C = 3.5\mathrm{x}10^{-5}\mathrm{Re}_{fo}^{0.44}\mathrm{Su}_{go}^{0.50}(\rho_f/\rho_g)^{0.48}\f$
\f$Re_f<2000\f$ | \f$Re_g>2000\f$ | \f$C = 0.0015\mathrm{Re}_{fo}^{0.59}\mathrm{Su}_{go}^{0.19}(\rho_f/\rho_g)^{0.36}\f$
\f$Re_f>2000\f$ | \f$Re_g<2000\f$ | \f$C = 8.7\mathrm{x}10^{-4}\mathrm{Re}_{fo}^{0.17}\mathrm{Su}_{go}^{0.50}(\rho_f/\rho_g)^{0.14}\f$
\f$Re_f>2000\f$ | \f$Re_g>2000\f$ | \f$C = 0.39\mathrm{Re}_{fo}^{0.03}\mathrm{Su}_{go}^{0.10}(\rho_f/\rho_g)^{0.35}\f$

The two-phase-multiplier is then given by
\f[
\phi_{f} = 1+\frac{C}{X}+\frac{1}{X^2}
\f]
and the frictional pressure gradient by
\f[
\left.\frac{dp}{dz}\right|_F = \phi_f \left.\frac{dp}{dz}\right|_{f}
\f]

*/
EXPORT_CODE double CONVENTION Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f(char *Fluid, double G, double Dh, double p, double beta, double x);

/// Homogeneous equilibrium model using the Adams formulation for the two-phase viscosity
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
EXPORT_CODE double CONVENTION HEM_DPDZ_f(char *Fluid, double G, double Dh, double p, double x);

/// Sung-Min Kim, Issam Mudawar, "Universal approach to predicting two-phase frictional pressure drop for 
/// mini/micro-channel saturated flow boiling", International Journal of Heat and Mass Transfer 58 (2013) 718–734
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param beta Aspect ratio (W/H) [-] (<0 for non-rectangular channel)
/// @param q_fluxH Heat flux based on heated perimeter[W/m^2]
/// @param PH_PF Ratio of heated perimeter to wetted perimeter [-]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
EXPORT_CODE double CONVENTION Kim_Mudawar_2013_Boiling_Microchannel_DPDZ_f(char *Fluid, double G, double Dh, double p, double beta, double q_fluxH, double PH_PF, double x);

/// A. Cavallini, D. Del Col, M. Matkovic, L. Rossetto, "Frictional pressure drop during vapour-liquid flow in minichannels: Modelling
/// and experimental evaluation", International Journal of Heat and Fluid Flow 30 (2009) 131–139
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
/// Erratum: log() in equation 8 for E is actually log_10, or base-10 logarithm, not natural logarithm
EXPORT_CODE double CONVENTION Cavallini_2009_AnnularMist_DPDZ_f(char *Fluid, double G, double Dh, double p, double x);

/// Frictional pressure drop of Lockhart-Martinelli in tubes
/// Lockhart, R.W., Martinelli, R.C., 1949, Proposed Correlation of Data for Isothermal Two-Phase 
/// Two-Component Flow in Pipes. Chemical Engineering Progress. v. 45, 39-48
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
EXPORT_CODE double CONVENTION Lockhart_Martinelli_1949_DPDZ_f(char *Fluid, double G, double Dh, double p, double x);

/// Friedel, L. (1979) Improved friction pressure drop correlations for horizontal and vertical 
/// two-phase pipe flow. European Two-Phase Flow Group Meeting, Ispra, Italy, paper E2
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
EXPORT_CODE double CONVENTION Friedl_1979_DPDZ_f(char *Fluid, double G, double Dh, double p, double x);

/// Stefan S. Bertsch, Eckhard A. Groll, Suresh V. Garimella, 
/// "A composite heat transfer correlation for saturated flow boiling in small channels"
/// International Journal of Heat and Mass Transfer 52 (2009) 2110-2118
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param p Pressure [Pa]
/// @param q_flux Heat flux based on heated perimeter[W/m^2]
/// @param L Length [m]
/// @param x Quality [-]
/// @return HTC Heat transfer coefficient [W/m^2/K]

/*!
Description
===========
The non-dimensional groups are given by:
\f[Re_{fo} = \frac{GD_h}{\mu_f}\f]
\f[Re_{go} = \frac{GD_h}{\mu_g}\f]
\f[\mathrm{Co} = \sqrt{\frac{\sigma}{g(\rho_f-\rho_g)D_h^2}}\f]

The nucleate boiling term is given by
\f[
\tilde h_{nb} = 55\left(\frac{p}{p_c}\right)^{0.12}  \left[-\log_{10}\left(\frac{p}{p_c}\right)\right]^{-0.55} M^{-0.5} (q''_H)^{0.67}
\f]
where \f$q''_H\f$ is for the heated portion of the perimeter.

The convective terms for liquid and vapor are given by
\f[
\tilde h_{conv,f} = \left(3.66+\frac{0.0668\frac{D_h}{L}\mathrm{Re}_{fo}\mathrm{Pr}_f}{1+0.04(\frac{D_h}{L}\mathrm{Re}_{fo}\mathrm{Pr}_f)^{2/3}}\right)\frac{k_f}{D_h}
\f]
\f[
\tilde h_{conv,g} = \left(3.66+\frac{0.0668\frac{D_h}{L}\mathrm{Re}_{go}\mathrm{Pr}_g}{1+0.04(\frac{D_h}{L}\mathrm{Re}_{go}\mathrm{Pr}_g)^{2/3}}\right)\frac{k_g}{D_h}
\f] 
The convective two-phase part is given by
\f[
\tilde h_{conv,2\phi} = \tilde h_{conv,f}(1-x)+\tilde h_{conv,g}x
\f]
Finally the overall two-phase heat transfer coefficient is given by
\f[
\tilde h_{2\phi} = \tilde h_{nb}(1-x)+\tilde h_{conv,2\phi}[1+80(x^2-x^6)\exp(-0.6\mathrm{Co})]
\f]
*/

EXPORT_CODE double CONVENTION Bertsch_2009_Boiling_Microchannel_HTC(char *Fluid, double G, double Dh, double p, double q_flux, double L, double x);

/// Accelerational pressure drop using the slip ratio definition of Zivi
/// @param Fluid The CoolProp-compliant fluid name
/// @param G Mass Flux [kg/m^2/s]
/// @param p Pressure [Pa]
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
EXPORT_CODE double CONVENTION Zivi_DPDZ_a(char *Fluid, double G, double p, double x1, double x2);
#endif