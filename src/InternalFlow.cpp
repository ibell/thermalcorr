#include "CoolProp.h"
#include "CPState.h"
#include "float.h"

double Kim_Mudawar_2013_Boiling_Microchannel_DPDZ_f(char *Fluid, double G, double Dh, double p, double beta, double q_fluxH, double PH_PF, double x)
{
	double f_f,f_g,C,Cnon_boiling;
	bool f_laminar, g_laminar;

	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;

	CPS.update(iP,p,iQ,x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s]
	double sigma = CPS.keyed_output(iI); //[N/m]

	double h_fg = (CPS.hV()-CPS.hL())*1000; // [J/kg]

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
double Kim_Mudawar_2012_AdiabaticCondensing_Microchannel_DPDZ_f(char *Fluid, double G, double Dh, double p, double beta, double x)
{
	double f_f,f_g,C;
	bool f_laminar, g_laminar;

	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;

	CPS.update(iP, p, iQ, x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s, or Pa-s]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s, or Pa-s]
	double sigma = CPS.keyed_output(iI); //[N/m]

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

EXPORT_CODE double CONVENTION HEM_DPDZ_f(char *Fluid, double G, double Dh, double p, double x)
{
	double f_tp;
	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;

	CPS.update(iP,p,iQ,x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s]
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

double Friedl_1979_DPDZ_f(char *Fluid, double G, double Dh, double p, double x)
{
	double f_fo,f_go;
	CoolPropStateClass CPS = CoolPropStateClass(Fluid);
	///TODO: Update units
	p /= 1000;

	CPS.update(iP, p, iQ, x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s]
	double Re_fo = G*Dh/mu_f;
	double Re_go = G*Dh/mu_g;
	double sigma = CPS.keyed_output(iI); //[N/m]
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

double Lockhart_Martinelli_1949_DPDZ_f(char *Fluid, double G, double Dh, double p, double x)
{
	double f_f,f_g,w,dpdz_f,dpdz_g,X,C,phi_g2,phi_f2;
	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;

	CPS.update(iP, p, iQ, x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s] 

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
double Cavallini_2009_AnnularMist_DPDZ_f(char *Fluid, double G, double Dh, double p, double x)
{
	double rho_GC,E_new,change;

	///TODO: Update units
	p /= 1000;

	CoolPropStateClass CPS = CoolPropStateClass(Fluid);
	CPS.update(iP,p,iQ,x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s]
	double sigma = CPS.keyed_output(iI); //[N/m]

	double pr = p / CPS.keyed_output(iPcrit); //[-]

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

double Shah_Condensation_HTC(char *Fluid, double G, double D, double p, double x)
{
	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;
	
	CPS.update(iP,p,iQ,x);

	double mu_f = CPS.viscL(); //[kg/m-s]
	double cp_f = CPS.cpL()*1000; //[J/kg-K]
	double k_f = CPS.condL()*1000; //[W/m-K]
	double Pr_f = cp_f * mu_f / k_f; //[-]
	double Pstar = p / CPS.keyed_output(iPcrit); //[-]

	// Liquid heat transfer coefficient
	double h_L = 0.023 * pow(G*D/mu_f,0.8) * pow(Pr_f,0.4) * k_f / D; //[W/m^2-K]

	// Liquid heat transfer coefficient multiplied by two-phase multiplier
	double HTC = h_L * (pow(1 - x, 0.8) + (3.8 * pow(x,0.76) * pow(1 - x,0.04)) / pow(Pstar,0.38) ); //[W/m^2-K]
	
	return  HTC;
}

double Bertsch_2009_Boiling_Microchannel_HTC(char *Fluid, double G, double Dh, double p, double q_flux, double L, double x)
{
	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;
	
	CPS.update(iP,p,iQ,x);

	double rho_f = CPS.rhoL();
	double mu_f = CPS.viscL(); //[kg/m-s]
	double cp_f = CPS.cpL()*1000; //[J/kg-K]
	double k_f = CPS.condL()*1000; //[W/m-K]
	double rho_g = CPS.rhoV();
	double mu_g = CPS.viscV(); //[kg/m-s]
	double cp_g = CPS.cpV()*1000; //[J/kg-K]
	double k_g = CPS.condV()*1000; //[W/m-K]
	double sigma = CPS.keyed_output(iI); //[N/m]
	double Pr_f = cp_f * mu_f / k_f; //[-]
	double Pr_g = cp_g * mu_g / k_g; //[-]
	double pr = p / CPS.keyed_output(iPcrit); //[-]
	double M = CPS.keyed_output(iMM); //[kg/kmol]
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

EXPORT_CODE double CONVENTION Zivi_DPDZ_a(char *Fluid, double G, double p, double x1, double x2)
{
	double term1, term2;

	CoolPropStateClass CPS = CoolPropStateClass(Fluid);

	///TODO: Update units
	p /= 1000;
	
	CPS.update(iP, p, iQ, x1);
	
	double rho_f = CPS.rhoL();
	double rho_g = CPS.rhoV();
	double v_f = 1/CPS.rhoL();
	double v_g = 1/CPS.rhoV();
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