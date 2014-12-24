#ifndef MICROCHANNEL_H
#define MICROCHANNEL_H

#include "CoolProp/CPState.h"

namespace ThermalCorr
{
namespace Microchannel
{

/// Sung-Min Kim, Issam Mudawar, "Universal approach to predicting two-phase frictional pressure drop for adiabatic
/// and condensing mini/micro-channel flows", International Journal of Heat and Mass Transfer 55 (2012) 3246–3261
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param beta Aspect ratio (W/H) [-] (set <0 for non-rectangular channel)
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
/*!
Description
===========
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
double Kim_Mudawar_2012_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double beta, double x);

/// Sung-Min Kim, Issam Mudawar, "Universal approach to predicting two-phase frictional pressure drop for 
/// mini/micro-channel saturated flow boiling", International Journal of Heat and Mass Transfer 58 (2013) 718–734
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
/// @param beta Aspect ratio (W/H) [-] (<0 for non-rectangular channel)
/// @param q_fluxH Heat flux based on heated perimeter[W/m^2]
/// @param PH_PF Ratio of heated perimeter to wetted perimeter [-]
/// @param x Quality [-]
/// @return dpdz Pressure gradient [Pa/m]
double Kim_Mudawar_2013_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double beta, double q_fluxH, double PH_PF, double x);

/// Stefan S. Bertsch, Eckhard A. Groll, Suresh V. Garimella, 
/// "A composite heat transfer correlation for saturated flow boiling in small channels"
/// International Journal of Heat and Mass Transfer 52 (2009) 2110-2118
/// @param CPS CoolPropStateClassSI instance
/// @param G Mass Flux [kg/m^2/s]
/// @param Dh Hydraulic Diameter [m]
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
double Bertsch_2009_HTC(CoolPropStateClassSI *CPS, double G, double Dh, double q_flux, double L, double x);
};
};
#endif