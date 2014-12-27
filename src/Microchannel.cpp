
#include "float.h"
#include "Microchannel.h"
#include "InternalFlow.h"

double ThermalCorr::Microchannel::Kim_Mudawar_2013_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double beta, double q_fluxH, double PH_PF, double x)
{
	double f_f,f_g,C,Cnon_boiling;
	bool f_laminar, g_laminar;

	double rho_f = (*CPS).rhoL();
	double mu_f = (*CPS).viscL(); //[kg/m-s]
	double rho_g = (*CPS).rhoV();
	double mu_g = (*CPS).viscV(); //[kg/m-s]
	double sigma = (*CPS).keyed_output(iI); //[N/m]

	double h_fg = ((*CPS).hV()-(*CPS).hL()); // [J/kg]

	double Re_f = G*(1-x)*Dh/mu_f;
	double Re_g = G*x*Dh/mu_g;

	if (Re_f < 2000)
	{
		if (beta>0){
			f_f = 24*(1-1.3553*beta+1.9467*beta*beta-1.7012*pow(beta,3)+0.9564*pow(beta,4)-0.2537*pow(beta,5))/Re_f;
		}
		else{
			f_f = 16/Re_f;
		}
		f_laminar = true;
	}
	else
	{
		f_laminar = false;
		if (Re_f<20000){
			f_f = 0.079*pow(Re_f,-0.25);
		}
		else{
			f_f = 0.046*pow(Re_f,-0.2);
		}
	}

	if (Re_g < 2000)
	{
		if (beta>0){
			f_g = 24*(1-1.3553*beta+1.9467*beta*beta-1.7012*pow(beta,3)+0.9564*pow(beta,4)-0.2537*pow(beta,5))/Re_g;
		}
		else{
			f_g = 16/Re_g;
		}
		g_laminar = true;
	}
	else
	{
		g_laminar = false;
		if (Re_g<20000){
			f_g = 0.079*pow(Re_g,-0.25);
		}
		else{
			f_g = 0.046*pow(Re_g,-0.2);
		}
	}

	double Re_fo = G*Dh/mu_f;
	double Su_go = rho_g*sigma*Dh/pow(mu_g,2);

	double dpdz_f = -2*f_f/rho_f*pow(G*(1-x),2)/Dh;
	double dpdz_g = -2*f_g/rho_g*pow(G*x,2)/Dh;

	double X_squared = dpdz_f/dpdz_g;

	// Find the C coefficient
	if (f_laminar && g_laminar){
		Cnon_boiling = 3.5e-5*pow(Re_fo,0.44)*pow(Su_go,0.50)*pow(rho_f/rho_g,0.48);
	}
	else if (f_laminar && !g_laminar){
		Cnon_boiling = 0.0015*pow(Re_fo,0.59)*pow(Su_go,0.19)*pow(rho_f/rho_g,0.36);
	}
	else if (!f_laminar && g_laminar){
		Cnon_boiling = 8.7e-4*pow(Re_fo,0.17)*pow(Su_go,0.50)*pow(rho_f/rho_g,0.14);
	}
	else if (!f_laminar && !g_laminar){
		Cnon_boiling = 0.39*pow(Re_fo,0.03)*pow(Su_go,0.10)*pow(rho_f/rho_g,0.35);
	}

	double We_fo = G*G*Dh/rho_f/sigma;
	double Bo = q_fluxH/(G*h_fg);

	if (Re_f >= 2000)
	{
		C = Cnon_boiling*(1+60*pow(We_fo,0.32)*pow(Bo*PH_PF,0.78));
	}
	else
	{
		C = Cnon_boiling*(1+530*pow(We_fo,0.52)*pow(Bo*PH_PF,1.09));
	}

	double two_phase_multiplier = 1+C/sqrt(X_squared)+1/X_squared;
	return dpdz_f*two_phase_multiplier;
}
double ThermalCorr::Microchannel::Kim_Mudawar_2012_DPDZ_f(CoolPropStateClassSI *CPS, double G, double Dh, double beta, double x)
{
	double f_f,f_g,C;
	bool f_laminar, g_laminar;

	double rho_f = (*CPS).rhoL();
	double mu_f = (*CPS).viscL(); //[kg/m-s, or Pa-s]
	double rho_g = (*CPS).rhoV();
	double mu_g = (*CPS).viscV(); //[kg/m-s, or Pa-s]
	double sigma = (*CPS).keyed_output(iI); //[N/m]

	double Re_f = G*(1-x)*Dh/mu_f;
	double Re_g = G*x*Dh/mu_g;

	if (Re_f < 2000)
	{
		if (beta>0){
			f_f = 24*(1-1.3553*beta+1.9467*beta*beta-1.7012*pow(beta,3)+0.9564*pow(beta,4)-0.2537*pow(beta,5))/Re_f;
		}
		else{
			f_f = 16/Re_f;
		}
		f_laminar = true;
	}
	else
	{
		f_laminar = false;
		if (Re_f<20000){
			f_f = 0.079*pow(Re_f,-0.25);
		}
		else{
			f_f = 0.046*pow(Re_f,-0.2);
		}
	}

	if (Re_g < 2000)
	{
		if (beta>0){
			f_g = 24*(1-1.3553*beta+1.9467*beta*beta-1.7012*pow(beta,3)+0.9564*pow(beta,4)-0.2537*pow(beta,5))/Re_g;
		}
		else{
			f_g = 16/Re_g;
		}
		g_laminar = true;
	}
	else
	{
		g_laminar = false;
		if (Re_g<20000){
			f_g = 0.079*pow(Re_g,-0.25);
		}
		else{
			f_g = 0.046*pow(Re_g,-0.2);
		} 
	}

	double Re_fo = G*Dh/mu_f;
	double Su_go = rho_g*sigma*Dh/pow(mu_g,2);

	double dpdz_f = -2*f_f/rho_f*pow(G*(1-x),2)/Dh;
	double dpdz_g = -2*f_g/rho_g*pow(G*x,2)/Dh;

	double X_squared = dpdz_f/dpdz_g;

	// Find the C coefficient
	if (f_laminar && g_laminar){
		C = 3.5e-5*pow(Re_fo,0.44)*pow(Su_go,0.50)*pow(rho_f/rho_g,0.48);
	}
	else if (f_laminar && !g_laminar){
		C = 0.0015*pow(Re_fo,0.59)*pow(Su_go,0.19)*pow(rho_f/rho_g,0.36);
	}
	else if (!f_laminar && g_laminar){
		C = 8.7e-4*pow(Re_fo,0.17)*pow(Su_go,0.50)*pow(rho_f/rho_g,0.14);
	}
	else if (!f_laminar && !g_laminar){
		C = 0.39*pow(Re_fo,0.03)*pow(Su_go,0.10)*pow(rho_f/rho_g,0.35);
	}

	double two_phase_multiplier = 1+C/sqrt(X_squared)+1/X_squared;
	return dpdz_f*two_phase_multiplier;
}
double ThermalCorr::Microchannel::Bertsch_2009_HTC(CoolPropStateClassSI *CPS, double G, double Dh, double q_flux, double L, double x)
{
	double rho_f = (*CPS).rhoL();
	double mu_f = (*CPS).viscL(); //[kg/m-s]
	double cp_f = (*CPS).cpL(); //[J/kg-K]
	double k_f = (*CPS).condL(); //[W/m-K]
	double rho_g = (*CPS).rhoV();
	double mu_g = (*CPS).viscV(); //[kg/m-s]
	double cp_g = (*CPS).cpV(); //[J/kg-K]
	double k_g = (*CPS).condV(); //[W/m-K]
	double sigma = (*CPS).keyed_output(iI); //[N/m]
	double Pr_f = cp_f * mu_f / k_f; //[-]
	double Pr_g = cp_g * mu_g / k_g; //[-]
	double pr = (*CPS).p() / (*CPS).keyed_output(iPcrit); //[-]
	double M = (*CPS).keyed_output(iMM); //[kg/kmol]
	double g = 9.81;

	double Re_fo = G*Dh/mu_f;
	double Re_go = G*Dh/mu_g;
	
	// Cooper correlation for the nucleate boiling contribution
    double h_nb = 55*pow(pr,0.12)*pow(-log10(pr),-0.55)*pow(M,-0.5)*pow(q_flux,0.67);
    double h_conv_l = (3.66+(0.0668*Dh/L*Re_fo*Pr_f)/(1+0.04*pow(Dh/L*Re_fo*Pr_f,2.0/3.0)))*k_f/Dh;
    double h_conv_g = (3.66+(0.0668*Dh/L*Re_go*Pr_g)/(1+0.04*pow(Dh/L*Re_go*Pr_g,2.0/3.0)))*k_g/Dh;
    double Co = sqrt(sigma/(g*(rho_f-rho_g)*Dh*Dh));
    
    // Things here depend on quality
	double h_conv_tp = h_conv_l*(1-x)+h_conv_g*x;
    double h_TP = h_nb*(1-x)+h_conv_tp*(1.0+80.0*(pow(x,2)-pow(x,6))*exp(-0.6*Co));

	return h_TP;
}