#include "BrazedPlateHeatExchanger.h"
#include "math.h"

void ThermalCorr::BrazedPlateHX::BPHE_1phase(BPHEGeometry geo, BPHEData &BPHE)
{  
    double X = 2*M_PI*geo.PlateAmplitude/geo.PlateWavelength;
    double PHI = 1.0/6.0*(1+sqrt(1+X*X)+4*sqrt(1+X*X/2));
    
    // The projected surface between the ports
    double A0 = geo.Bp*geo.Lp;
    
    // The plane surface of one plate
    double Ap = PHI*A0;
    
    // The volume of one channel
    double Vchannel = geo.Bp*geo.Lp*2*geo.PlateAmplitude;
    
    // Hydraulic diameter
    double dh = 4*geo.PlateAmplitude/PHI;
    
	double rho_g = BPHE.AS->rhomass();
	double k_g = BPHE.AS->conductivity();
	double eta_g = BPHE.AS->viscosity();
    double Pr_g = BPHE.AS->cpmass()*BPHE.AS->viscosity()/BPHE.AS->conductivity();
    
    double eta_g_w = eta_g; // TODO: allow for temperature dependence?
    double w_g = BPHE.mdot_per_channel/BPHE.AS->rhomass()/(2*geo.PlateAmplitude*geo.Bp);
    double Re_g = rho_g*w_g*dh/eta_g;
    
    // Calculate the friction factor zeta
    double phi = geo.InclinationAngle;
    
	double zeta0, zeta1_0, zeta1, a, b, c;
	if (Re_g < 2000)
	{
        zeta0 = 64/Re_g;
        zeta1_0 = 597/Re_g+3.85;
	}
    else
	{
        zeta0 = pow(1.8*log(Re_g)-1.5,-2);
        zeta1_0 = 39/pow(Re_g,0.289);
	}
    
    a = 3.8;
    b = 0.18;
    c = 0.36;
    
    zeta1 = a*zeta1_0;
    // RHS from Equation 18
    double RHS = cos(phi)/sqrt(b*tan(phi)+c*sin(phi)+zeta0/cos(phi))+(1-cos(phi))/sqrt(zeta1);
    double zeta = 1/(RHS*RHS);
    // Hagen number
    double Hg = zeta*Re_g*Re_g/2;
    
    // Constants for Nu correlation
    double c_q = 0.122;
    double q = 0.374;

    // Nusselt number [-]
    double Nu = c_q*pow(Pr_g,1.0/3.0)*pow(eta_g/eta_g_w,1.0/6.0)*pow(2*Hg*sin(2*phi),q);
    
    // Heat transfer coefficient [W/m^2-K]
    BPHE.HTC = Nu*k_g/dh;
    
    // Pressure drop [Pa]
    BPHE.DELTAP = Hg*eta_g*eta_g*geo.Lp/(rho_g*dh*dh*dh);
}