
#include "float.h"
#include "InternalFlow.h"

namespace ThermalCorr{
namespace GeneralInternal{
double HEM_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x)
{
	double f_tp;

	double rho_f = (*CPS).rhoL();
	double mu_f = (*CPS).viscL(); //[kg/m-s]
	double rho_g = (*CPS).rhoV();
	double mu_g = (*CPS).viscV(); //[kg/m-s]
	double v_f = 1/rho_f;
	double v_fg = 1/rho_g-1/rho_f;

	double mu_tp = 1/(x/mu_g+(1-x)/mu_f);
	double Re_tp = G*Dh/mu_tp;
	if (Re_tp < 2000)
	{
		f_tp = 16/Re_tp;
	}
	else if (Re_tp >= 2000 && Re_tp < 20000)
	{
		f_tp = 0.079*pow(Re_tp,-0.25);
	}
	else
	{
		f_tp = 0.046*pow(Re_tp,-0.2);
	}
	return -2*f_tp*v_f*G*G/Dh*(1+x*v_fg/v_f);
}

double Friedl_1979_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x)
{
	double f_fo,f_go;

	double rho_f = (*CPS).rhoL();
	double mu_f = (*CPS).viscL(); //[kg/m-s]
	double rho_g = (*CPS).rhoV();
	double mu_g = (*CPS).viscV(); //[kg/m-s]
	double Re_fo = G*Dh/mu_f;
	double Re_go = G*Dh/mu_g;
	double sigma = (*CPS).keyed_output(iI); //[N/m]
	double rhoH = 1/(x/rho_g+(1-x)/rho_f);

	if (Re_fo < 2000){
		f_fo = 16/Re_fo;
	}
	else{
		if (Re_fo<20000){
			f_fo = 0.079*pow(Re_fo,-0.25);
		}
		else{
			f_fo = 0.046*pow(Re_fo,-0.2);
		}
	}

	if (Re_go < 2000){
		f_go = 16/Re_go;
	}
	else{
		if (Re_go < 20000){
			f_go = 0.079*pow(Re_go,-0.25);
		}
		else{
			f_go = 0.046*pow(Re_go,-0.2);
		} 
	}

	double H = pow(rho_f/rho_g,0.91)*pow(mu_g/mu_f,0.19)*pow(1-mu_g/mu_f,0.7);
	double F = pow(x,0.78)*pow(1-x,0.224);
	double E = pow(1-x,2)+pow(x,2)*(rho_f*f_go)/(rho_g*f_fo);
	double Fr = G*G/(9.813*Dh*rhoH);
	double We = G*G*Dh/(rhoH*sigma);

	double dpdz_f = -2*f_fo*pow(G,2)/rho_f/Dh;
	double two_phase_multiplier = E+3.24*F*H/(pow(Fr,0.045)*pow(We,0.035));

	return dpdz_f*two_phase_multiplier;
}

double Lockhart_Martinelli_1949_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x)
{
	double f_f,f_g,w,dpdz_f,dpdz_g,X,C,phi_g2,phi_f2;

	double rho_f = (*CPS).rhoL();
	double mu_f = (*CPS).viscL(); //[kg/m-s]
	double rho_g = (*CPS).rhoV();
	double mu_g = (*CPS).viscV(); //[kg/m-s] 

	// 1. Find the Reynolds Number for each phase based on the actual flow rate of the individual phase
	double Re_f = G*(1-x)*Dh/mu_f;
	double Re_g = G*x*Dh/mu_g;

    // 2. Friction factor for each phase
	if (fabs(1-x)<10*DBL_EPSILON){ //No liquid
        f_f=0; //Just to be ok until next step, otherwise a divide by zero error
	}
	else if (Re_f<1000){ //#Laminar
        f_f=16/Re_f;
	}
	else if (Re_f>2000){ //Fully-Turbulent
        f_f=0.046/pow(Re_f,0.2);
	}
	else{ //Mixed
        // Weighting factor
        w=(Re_f-1000)/(2000-1000);
        // Linear interpolation between laminar and turbulent
        f_f=w*16.0/Re_f+(1-w)*0.046/pow(Re_f,0.2);
	}

	// 2a. Friction factor for gas
	if (fabs(x)<10*DBL_EPSILON){ //No gas, just liquid
        f_g=0; //Just to be ok until next step
	}
	else if (Re_g<1000){ //Laminar
        f_g=16.0/Re_g;
	}
	else if (Re_g>2000){ //Fully-Turbulent
        f_g=0.046/pow(Re_g,0.2);
	}
	else{ //Mixed
        // Weighting factor
        w = (Re_g-1000)/(2000-1000);
        // Linear interpolation between laminar and turbulent
        f_g = w*16.0/Re_g+(1-w)*0.046/pow(Re_g,0.2);
	}

    // 3. Frictional pressure drop based on actual flow rate of each phase
    dpdz_f = -2*f_f*pow(G*(1-x),2)/rho_f/Dh;
    dpdz_g = -2*f_g*pow(G*x,2)/rho_g/Dh;

	if (x <= 10*DBL_EPSILON){
        // Entirely liquid
        return dpdz_f;
	}
	if (x >= 1-10*DBL_EPSILON){
        // Entirely vapor
        return dpdz_g;
	}

    // 4. Lockhart-Martinelli parameter
    X=sqrt(dpdz_f/dpdz_g);

    // 5. Find the Constant based on the flow Re of each phase
    //    (using 1500 as the transitional Re to ensure continuity)
    
    // Calculate C
	if (Re_f>1500 && Re_g > 1500){
        C=20.0;
	}
	else if (Re_f<1500 && Re_g>1500){
        C=12.0;
	}
	else if (Re_f>1500 && Re_g<1500){
        C=10.0;
	}
	else{
        C=5.0;
	}

    // 6. Two-phase multipliers for each phase
    phi_g2=1+C*X+X*X;
    phi_f2=1+C/X+1/X/X;
        
    // 7. Find gradient
	if (dpdz_g*phi_g2>dpdz_f*phi_f2){
        return dpdz_g*phi_g2;
	}
	else{
        return dpdz_f*phi_f2;
	}
}
double Cavallini_2009_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double x)
{
	double rho_GC,E_new,change;

	double rho_f = (*CPS).rhoL(); // [kg/m^3]
	double mu_f = (*CPS).viscL(); // [kg/m-s]
	double rho_g = (*CPS).rhoV(); // [kg/m^3]
	double mu_g = (*CPS).viscV(); // [kg/m-s]
	double sigma = (*CPS).keyed_output(iI); //[N/m]

	double pr = (*CPS).p() / (*CPS).keyed_output(iPcrit); //[-]

	double Re_fo = G*Dh/mu_f;

	double W = 1.398*pr;
    double f_fo = 0.046*pow(Re_fo,-0.2);
    double H = pow(rho_f/rho_g,1.132)*pow(mu_g/mu_f,0.44)*pow(1-mu_g/mu_f,3.542);
    double g = 9.81;
    double J_G = x*G/sqrt(g*Dh*rho_g*(rho_f-rho_g));

    double Z = (1-x)*(1-x)+x*x*rho_f/rho_g*pow(mu_g/mu_f,0.2);
    double F = pow(x,0.9525)*pow(1-x,0.414);
    double j_G = x*G/rho_g;
	double E = 0.4;
	
	do
	{
		rho_GC = (x+(1-x)*E)/(x/rho_g+(1-x)*E/rho_f);
		E_new = 0.015+0.44*log10(rho_GC/rho_f*pow(mu_f*j_G/sigma,2)*1e4);
		change = E_new-E;
		E = E_new;
	}while (fabs(change)>1e-6);

	if (E>0.95){
		E = 0.95;
	}
	else if (E<0.0){
		E = 0.0;
	}

	double two_phase_multiplier =  Z+3.595*F*H*pow(1-E,W);
	double dpdz_fo = -2*f_fo*G*G/Dh/rho_f;
	return dpdz_fo*two_phase_multiplier;
}

double Shah_1976_HTC(CoolPropStateClassSI *CPS, double G, double D, double x)
{

	double mu_f = (*CPS).viscL(); //[kg/m-s]
	double cp_f = (*CPS).cpL(); //[J/kg-K]
	double k_f = (*CPS).condL(); //[W/m-K]
	double Pr_f = cp_f * mu_f / k_f; //[-]
	double Pstar = (*CPS).p() / (*CPS).keyed_output(iPcrit); //[-]

	// Liquid heat transfer coefficient
	double h_L = 0.023 * pow(G*D/mu_f,0.8) * pow(Pr_f,0.4) * k_f / D; //[W/m^2-K]

	// Liquid heat transfer coefficient multiplied by two-phase multiplier
	double HTC = h_L * (pow(1 - x, 0.8) + (3.8 * pow(x,0.76) * pow(1 - x,0.04)) / pow(Pstar,0.38) ); //[W/m^2-K]
	
	return  HTC;
}

double Zivi_DPDZ_a(CoolPropStateClassSI *CPS, double G, double x1, double x2)
{
	double term1, term2;
	
	double rho_f = (*CPS).rhoL();
	double rho_g = (*CPS).rhoV();
	double v_f = 1/(*CPS).rhoL();
	double v_g = 1/(*CPS).rhoV();
	double S = pow(v_g/v_f,1.0/3.0);
	
	// Void fraction at x1
	double alpha1 = 1/(1+v_f/v_g*S*(1-x1)/x1);
	// Void fraction at x2
	double alpha2 = 1/(1+v_f/v_g*S*(1-x2)/x2);

	if (fabs(x1) < 10*DBL_EPSILON){
		term1 = v_f;
	}
	else if (fabs(x1-1) < 10*DBL_EPSILON){
		term1 = v_g;
	}
	else{
		term1 = x1*x1*v_g/alpha1 + (1-x1)*(1-x1)*v_f/(1-alpha1);
	}

	if (fabs(x2) < 10*DBL_EPSILON){
		term2 = v_f;
	}
	else if (fabs(x2-1) < 10*DBL_EPSILON){
		term2 = v_g;
	}
	else{
		term2 = x2*x2*v_g/alpha2 + (1-x2)*(1-x2)*v_f/(1-alpha2);
	}

	return G*G*(term1-term2);
}
};
};
