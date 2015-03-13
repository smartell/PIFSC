//  ******************************************************************
//  Statistical Catch Age Model (SCAM)
//
//  Created by Steve Martell on 2010-04-09.
//  Copyright (c) 2010. All rights reserved.
//  Comments:
//   - At present there is a bias in the estimate of Rbar b/c
//   - initial-age's should be set up with Rinit
//   
//  TODO:
//       1) build control file for leading parameters
//       2) use of ivectors to index array in surveys
//       3) Uncertainty estimates
//       4) Reference points & forecasts
//  ******************************************************************
DATA_SECTION

	int sim;
	int rseed;
	LOC_CALCS
		sim=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-sim" command line option
		//if it exists the if statement retreives the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
			cout<<"______________________________________________________\n"<<endl;
			cout<<"    **Implementing Simulation--Estimation trial** "<<endl;
			cout<<"______________________________________________________"<<endl;
			cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
			cout<<"______________________________________________________\n"<<endl;
			//if(sim)exit(1);
		}
	END_CALCS
	
	//Read in objects from data file using the init_ prefix
	init_int syr;
	init_int nyr;
	init_int nage;
	vector age(1,nage);
	!! age.fill_seqadd(1,1);
	
	init_number m;
	init_number linf;
	init_number k;
	init_number a;
	init_number b;
	init_number ah;
	init_number gh;
	vector fa(1,nage);		//fecundity-at-age
	vector la(1,nage);		//length-at-age
	vector wa(1,nage);		//weight-at-age
	LOC_CALCS
	  la=linf*(1.-exp(-k*age));
	  wa=a*pow(la,b);
	  fa=elem_prod(plogis(age,ah,gh),wa);
	  cout<<setprecision(2);		//2 decimal places for output
	  cout<<"la\n"<<la<<endl;
	  cout<<"wa\n"<<wa<<endl;
	  cout<<"fa\n"<<fa<<endl;
	  cout<<setprecision(5);
	END_CALCS
	
	//Time series data
	init_vector obs_ct(syr,nyr);
	init_int nit;
	init_ivector iyr(1,nit);
    init_vector it(1,nit);

	//Age comps
	init_int p_sage;
	init_int p_nage;
	init_matrix P(syr,nyr,p_sage,p_nage);
	
	init_int q_sage;
	init_int q_nage;
	init_matrix Q(1,nit,q_sage,q_nage);
	
	//End of data file
	init_int eof;
	int nodes;
	!! nodes=nage-6;
	vector x(1,nodes);
	vector fage(1,nage-1);
	!! x.fill_seqadd(0,1./(1-nodes));  
	!! fage.fill_seqadd(0,1./(1-(nage-1)));
	
	LOC_CALCS
	  if(eof==999){
		cout<<"\n -- END OF DATA SECTION -- \n"<<endl;
	  }else{
		cout<<"\n ** ERROR READING DATA ** \n"<<endl; exit(1);
	  }
	END_CALCS
	
	number tau2_q;			//mle of variance for survey age comps
	number tau2_p;			//mle of variance for fishery age comps
	
	!! ad_comm::change_datafile_name("scam.ctl");
	init_vector lb_theta(1,4);
	init_vector	ub_theta(1,4);
	init_vector in_theta(1,4);
	init_ivector ph_theta(1,4);
	

INITIALIZATION_SECTION
	theta in_theta;
	//
	//log_ro      4.60;
	//log_avgrec  4.60;
	//log_kappa   3.50;
	//log_avg_f  -1.60;
    //

PARAMETER_SECTION
	number log_ro;				//log of unfished recruitment
	number log_kappa;			//log of recruitment compensation
	number log_avgrec;			//log of average recruitment
	number log_avg_f;			//log of average fishing mortality
	
	init_bounded_number_vector theta(1,4,lb_theta,ub_theta,ph_theta);

	init_bounded_dev_vector log_rec_devs(syr-nage-1,nyr,-10.,10,1);
	init_bounded_dev_vector log_ft_devs(syr,nyr,-15.,15.,2);
	init_bounded_dev_vector log_sel_coffs(1,nage-1,-15.,15.,1);
	init_number_vector log_sel_vec(1,2,1);	//survey selectivity parameters
	
	
	objective_function_value f;
    
	number ro;					//unfished age-1 recruits
	number bo;					//unfished spawning stock biomass
	number kappa;
	number so;
	number beta;
	
	vector log_rt(syr-nage-1,nyr);
	vector vax(1,nage);			//survey selectivity coefficients
	vector ct(syr,nyr);			//predicted catch biomass
	vector sbt(syr,nyr+1);		//spawning stock biomass
	
	vector rt(syr+1,nyr); 		//predicted age-1 recruits from S-R curve
	
	vector delta(syr+1,nyr);	//residuals for stock recruitment
	vector epsilon(1,nit);		//residuals for survey abundance index
	
	matrix log_sel(syr,nyr,1,nage);		//selectivity coefficients
	matrix N(syr,nyr+1,1,nage);
	matrix F(syr,nyr,1,nage);
	matrix Z(syr,nyr,1,nage);
	matrix S(syr,nyr,1,nage);
	
	matrix Phat(syr,nyr,p_sage,p_nage);		//predicted fishery age proportions
	matrix Qhat(1,nit,q_sage,q_nage);		//predicted survey age proportions
	
	sdreport_number sd_depletion;
	sdreport_vector sd_noF_sbt(syr,nyr+1);
	
PRELIMINARY_CALCS_SECTION
  //Run the model with input parameters to simulate real data.
  if(sim) simulation_model(rseed);

RUNTIME_SECTION
    maximum_function_evaluations 100,100,500,5000,5000
    convergence_criteria 0.01,0.01,1.e-4,1.e-4

PROCEDURE_SECTION
	initialize_params();
	
	calc_mortality();
	
	calc_numbers_at_age();
	
	calc_fishery_observations();
	
	calc_survey_observations();
	
	calc_stock_recruitment();
	
	calc_objective_function();
	
	if(sd_phase())
	{
		sd_depletion=sbt(nyr+1)/bo;
		//call your function here
		sb_with_noF();
	}

FUNCTION initialize_params
	log_ro     = theta(1);
	log_kappa  = theta(2);
	log_avgrec = theta(3);
	log_avg_f  = theta(4);
	


FUNCTION sb_with_noF
	// Good luck!
	S = mfexp(-m);
	calc_numbers_at_age();
	sd_noF_sbt = N*fa; //elem_div(N*fa,sbt);
	
	S = mfexp(-Z);
	calc_numbers_at_age();
	
FUNCTION calc_mortality
	int i;
	dvar_matrix tmp(syr,nyr-1,1,nage-1);  //infrastructure for time-varying selectivity
	log_sel.initialize();
	tmp.initialize();
	log_sel(syr)(1,nage-1)    = log_sel_coffs;
	log_sel(syr)(nage-1,nage) = log_sel_coffs(nage-1);
		
	for(i=syr;i<nyr;i++)
	{
		log_sel(i+1)(1,nage-1)    = log_sel(i)(1,nage-1)+tmp(i);
		log_sel(i+1)(nage-1,nage) = log_sel(i+1,nage-1);
	}
	
	
	for(i=syr;i<=nyr;i++)
	{
		log_sel(i) -= log(mean(mfexp(log_sel(i))));
		F(i) = mfexp(log_avg_f+log_ft_devs(i)+log_sel(i));
	}
	
	
	Z=m+F;
	S=mfexp(-Z);
	
	
FUNCTION calc_numbers_at_age
	int i,j;
	
	// annual recruitment
	log_rt=log_avgrec+log_rec_devs; 
	for(i=syr;i<=nyr;i++)
	{
		N(i,1) = mfexp(log_rt(i));
	} 
	
	// initial numbers-at-age
	for(j=2;j<=nage;j++)
	{
		N(syr,j)=mfexp(log_rt(syr-j))*exp(-m*(j-1));
	}
	N(syr,nage)/=(1.-exp(-m));
	
	// update numbers-at-age
	for(i=syr;i<=nyr;i++)
	{
		N(i+1)(2,nage)=++elem_prod(N(i)(1,nage-1),S(i)(1,nage-1));
		N(i+1,nage)+=N(i,nage)*S(i,nage);
	}
	
FUNCTION calc_fishery_observations
	int i;
	
	//Catch-at-age
	//C = F/Z*(1-S)*N;
    dvar_matrix C = elem_prod(elem_prod(elem_div(F,Z),1.-S),N);
	
	//Proportions-at-age
	for(i=syr;i<=nyr;i++)
		Phat(i) = C(i)(p_sage,p_nage)/sum(C(i)(p_sage,p_nage));
	
	//predicted total catch
	ct=C*wa;
	

FUNCTION calc_survey_observations
	int i,ii;
	vax = plogis(age,mfexp(log_sel_vec(1)),mfexp(log_sel_vec(2))); 
	
	//survey age comps
	dvar_matrix V(1,nit);
	for(i=1;i<=nit;i++)
	{   
		ii=iyr(i);
		V(i)=elem_prod(N(ii),vax);
		Qhat(i) = V(i)(q_sage,q_nage)/sum(V(i)(q_sage,q_nage));
	}
	
	//survey abudance index
	dvar_vector zt=log(it)-log(V*wa);
	epsilon = zt-mean(zt);
	
FUNCTION calc_stock_recruitment
	/*
	The following code is used to derive unfished
	spawning stock biomass bo and the stock-
	recruitment parameters for the Beverton-Holt
	model.  Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta)
	
	This is associated with slide 
	*/ 
	int i;
	ro=mfexp(log_ro);
	kappa=mfexp(log_kappa)+1.001;
	
	dvector lx(1,nage); lx=1;
	for(i=2;i<=nage;i++) lx(i)=lx(i-1)*exp(-m);
	lx(nage)/=(1.-exp(-m));
	double phib = lx * fa;
	bo = ro*phib;  
	
	sbt=N*fa;
	dvar_vector t1=kappa*ro*sbt(syr,nyr-1);
	dvar_vector t2=(bo+(kappa-1.)*sbt(syr,nyr-1));
	dvar_vector tmp_rt=++elem_div(t1,t2);		//predicted recruits
	
	//residuals in stock-recruitment curve
	rt=exp(log_rt(syr+1,nyr));//trans(N)(1)(syr+1,nyr);
	delta = log(rt)-log(tmp_rt);
	
FUNCTION calc_objective_function
	/*
	There are several components to the objective function
	Likelihoods:
		-1) likelihood of the catch data
		-2) likelihood of the survey abundance index
		-3) likelihood of the survey age comps
		-4) likelihood of the fishery age comps
		-5) likelihood for stock-recruitment relationship
		-6) likelihood for fishery selectivities
	*/
	int i,j;
	dvar_vector lvec(1,6); lvec.initialize();
	
	//1) likelihood of the catch data
	lvec(1)=dnorm(log(obs_ct)-log(ct),0.05);
	
	//2) likelihood of the survey abundance index
	lvec(2)=dnorm(epsilon,0.25);
	
	//3) likelihood of the survey age comps
	lvec(3)=dmvlogistic(Q,Qhat,tau2_q);
	
	//4) likelihood of the survey age-comps
	lvec(4)=dmvlogistic(P,Phat,tau2_p);
	
	//5) likelihood for stock-recruitment relationship
	lvec(5)=dnorm(delta,0.9);
	
	//6) likelihood for fishery selectivity paramters
	dvar_vector df2=first_difference(first_difference(log_sel(syr)));
	lvec(6)=12.5*norm2(df2);
 
	
	/*
	The following are penalties that are invoked in early
	phases to ensure reasonable solutions are obtained,
	and to reduce the sensitivity of initial starting
	conditions.  Penalties include:
		-1) keep average fishing mortality rate near 
			0.2 and in the last phase relax this constraint.
		-2) normal prior on fishing mortality deviates with
			a large standard deviation of 50.
		-3) normal prior for log rec devs with std=50.
	*/
	
	dvar_vector pvec(1,5);
	pvec.initialize();
	
	//Penalties to regularize the solution for fishing mortality rates
	if(last_phase())
	{
		pvec(1)=0.001*square(log_avg_f-log(0.2));
		pvec(2) = dnorm(log_ft_devs,50.);
		pvec(3) = dnorm(log_rec_devs,5.);
	}
	else
	{
		pvec(1)=1000.*square(log_avg_f-log(0.2));
		pvec(2)=100.*norm2(log_ft_devs);
		pvec(3)=100.*norm2(log_rec_devs);
	}
	
	
	//beta prior for steepness
	dvariable h=kappa/(4.+kappa);
	pvec(5) = dbeta((h-0.2)/0.8,5,3); 
		
	f=sum(lvec)+sum(pvec);

FUNCTION void simulation_model(const long& seed)
	/*
	Call this routine to simulate data for simulation testing.
	The random number seed can be used to repeat the same 
	sequence of random number for simulation testing.
	
	Implemented using the "-sim 99" command line option
	
	-This routine will over-write the observations in memory
	with simulated data, where the true parameter values are 
	the initial values.  Change the standard deviations of the 
	random number vectors epsilon (observation error) or 
	recruitment devs wt (process error).
	*/
	int i,j,ii;
	dmatrix C(syr,nyr,1,nage);
	
	//Random number generator
	random_number_generator rng(seed);
	dvector wt(syr-nage-1,nyr);			//recruitment anomalies
	dvector epsilon(1,nit);			//observation errors in survey
	wt.fill_randn(rng); wt *= 0.4;
	epsilon.fill_randn(rng); epsilon *= 0.25;
	
    //Initialize model
	initialize_params();
	
	ro          = exp(log_ro);
	kappa       = exp(log_kappa)+1.001;
	dvector lx  = pow(exp(-m),age-1.);
	lx(nage)   /=(1.-exp(-m));
	double phie = lx*fa;		//this is the sum of products for two vectors
	so          = kappa/phie;
	beta        = (kappa-1.)/(ro*phie);
	
	//Initial numbers-at-age with recruitment devs
	N(syr,1)=exp(log_avgrec+wt(syr));
	for(j=2;j<=nage;j++)
	{
		N(syr,j)=exp(log_avgrec+wt(syr-j))*lx(j);
	}
	
	dvector va=plogis(age,4,1.5);				//fishery selectivity 
	//dvector va=eplogis(age,1/1.5,4,0.1);		//fishery selectivity domed
	dvector vax=plogis(age,4.5,1.8);			//survey selectivity
	for(i=syr;i<=nyr;i++)
	{   
		//Approximate fishing mortality with Popes approximation
		//double btmp = value(N(i))*exp(-m/2.)*elem_prod(va,wa);
		//double ftmp = obs_ct(i)/btmp;
		
		dvector bt = elem_prod(value(N(i)),wa);
		double ft = get_ft(obs_ct(i),m,va,bt);
		dvector zt = m+ft*va;
		
		
		//Update numbers at age
		double et=value(N(i))*fa;
		N(i+1,1)=so*et/(1.+beta*et)*exp(wt(i));
		N(i+1)(2,nage)=++elem_prod(N(i)(1,nage-1),exp(-zt(1,nage-1)));
		N(i+1,nage)+=N(i,nage)*exp(-zt(nage));
		
		//Catch & Catch-at-age
		C(i) = elem_prod(elem_div(ft*va,zt),elem_prod(1.-exp(-zt),value(N(i))));
		
		//Proportions-at-age in commercial samples
		P(i) = rmvlogistic(C(i)(p_sage,p_nage),0.3,i+seed);
		
	}
	ct = C*wa;	//total catch in biomass
	/*
		Calculations for survey time series
	*/
	for(i=1;i<=nit;i++)
	{   
		ii=iyr(i);
		it(i) = (value(N(ii))*elem_prod(wa,vax))*exp(epsilon(i));
		
		//Proportions-at-age in fisheries independent survey
		Q(i) = rmvlogistic(elem_prod(value(N(ii)(q_sage,q_nage)),vax(q_sage,q_nage)),0.2,i+seed);
	}
	
	cout<<"Simulated exploitation rate\n"<<elem_div(ct,(N.sub(syr,nyr)*elem_prod(wa,va)))<<endl;
	
	
REPORT_SECTION
	REPORT(ro);
	double rbar=value(exp(log_avgrec));
	REPORT(rbar);
	REPORT(bo);
	REPORT(kappa);
	
	REPORT(tau2_q);
	REPORT(tau2_p);
	
	ivector yr(syr,nyr);
	ivector yrs(syr,nyr+1);
	yr.fill_seqadd(syr,1); 
	yrs.fill_seqadd(syr,1); 
	REPORT(yr);
	REPORT(yrs);
	REPORT(iyr);
	REPORT(age);
	REPORT(la);
	REPORT(wa);
	REPORT(fa);
	REPORT(log_sel);
	REPORT(vax);
	
	REPORT(obs_ct);
	REPORT(ct);
	dvector ft=value(mfexp(log_avg_f+log_ft_devs));
	REPORT(ft);
	report<<"bt\n"<<N*wa<<endl;
	report<<"sbt\n"<<sbt<<endl;
	REPORT(sd_noF_sbt);
	REPORT(rt);
	REPORT(log_rt(syr,nyr));
	REPORT(delta);
	REPORT(it);
	dvector pit=value(N.sub(syr,nyr)*elem_prod(vax,wa));
	REPORT(pit(iyr));
	REPORT(epsilon);
	REPORT(F);
	REPORT(P); REPORT(Phat);
	REPORT(Q); REPORT(Qhat);


TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
    
	#undef COUT
	#define COUT(object) cout << #object "\n" << object <<endl;

	#include <admodel.h>
	#include <time.h>
	//#include <stats.cxx>
	#include <baranov.cxx>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
	dvariable plogis(const dvariable& x, const double& mu, const dvariable& std)
	{
		return 1./(1.+mfexp((mu-x)/std));
	}

	dvar_vector plogis(const dvector& x, const dvariable& mu, const dvariable& std)
	{
		return 1./(1.+mfexp((mu-x)/std));
	}

	dvector plogis(const dvector& x, const double& mu, const double& std)
	{
		return 1./(1.+mfexp((mu-x)/std));
	}

	dvar_vector plogis(const dvar_vector& x, const dvariable& mu, const dvariable& std)
	{
		return 1./(1.+mfexp((mu-x)/std));
	}




FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;

