DATA_SECTION
	init_adstring datafile;
	!!COUT(datafile);
	!! ad_comm::change_datafile_name(datafile);



	int rseed;
	int mseed;
	LOC_CALCS
		// Command line options.
		rseed = 0;
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			rseed = atoi(ad_comm::argv[on+1]);
		}

		mseed = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt))>-1)
		{
			mseed = atoi(ad_comm::argv[on+1]);
		}

	END_CALCS



	init_int nyrs;
	init_number alpha;
	init_number rho;
	init_number wk;
	init_int agek;


	init_matrix data(1,nyrs,1,4);
	//!! cout<<data<<endl;
	ivector year(1,nyrs);
	vector ct(1,nyrs);
	vector cpue(1,nyrs);
	vector wt(1,nyrs);
	
	LOC_CALCS
		year = ivector(column(data,1));
		ct   = column(data,2);
		cpue = column(data,3);
		wt   = column(data,4);
	END_CALCS
	
	// ADDING MSE FRAMEWORK
	//friend_class OperatingModel;


INITIALIZATION_SECTION
	log_bo   3.65;
	log_reck 2.48;
	log_m   -1.60;
	log_sigma_epsilon  -2.60;
	log_sigma_nu       -2.20;
	log_sigma_delta    -4.60;
	log_tau            -1.60;
	log_fbar           -1.60;
	log_tau            -1.60;
	


PARAMETER_SECTION
	init_number log_bo;
	init_number log_reck;
	init_number log_m;
	init_number log_sigma_epsilon(2);
	init_number log_sigma_nu(4);
	init_number log_sigma_delta(5);
	init_number log_tau(3);
	init_number log_fbar(3);
	
	init_bounded_dev_vector psi(1,nyrs,-5,5,3);
	init_bounded_dev_vector fdev(1,nyrs,-5,5);

	objective_function_value f;
	number bo;
	number ro;
	number no;
	number reck;
	
	number m;
	number sigma_epsilon;
	number sigma_nu;
	number sigma_delta;
	number tau;

	number so;
	number beta;
	number lnq;

	vector ft(1,nyrs);
	vector st(1,nyrs);
	vector bt(1,nyrs);
	vector nt(1,nyrs);
	vector rt(1,nyrs);
	vector qt(1,nyrs);
	vector what(1,nyrs);
	vector chat(1,nyrs);
	vector nu(1,nyrs);
	vector delta(1,nyrs);
	vector epsilon(1,nyrs);

	vector prior_vec(1,8);

	


PRELIMINARY_CALCS_SECTION
	if(rseed != 0)
	{
		cout<<"Simulating fake data with random seed "<<rseed<<endl;
		runSimulationModel();
		//exit(1);
	}


PROCEDURE_SECTION

	initialStates();
	calcFishingMortality();
	populationDynamics();
	observationModels();
	calculatePriors();
	calcObjectiveFunction();


FUNCTION runSimulationModel
	/**
	This simulation model is based on the same routines used 
	in the assessment model.  

	Pseudocode:
		- fill vectors with random normal deviates.
		- call functions to calculate predicted observations
		- add random error to the simulated observations
		- continue with running the assessment and estimate
		  model parameters.
		- write estimated parameters to a file for comparison
		  with the true values (pin file) that were used to 
		  simulate the model data.



	*/
	
	// Random number generator class
	random_number_generator rng(rseed);

	// vectors for random numbers must be data types, not dvar_...
	dvector repsilon(1,nyrs);
	dvector rnu(1,nyrs);
	dvector rdelta(1,nyrs);
	dvector rpsi(1,nyrs);

	// fill with random normal deviates, mean =0, sd = 1
	repsilon.fill_randn(rng);
	rnu.fill_randn(rng);
	rdelta.fill_randn(rng);
	rpsi.fill_randn(rng);

	//COUT(repsilon);

	// call functions to simulate data.
	initialStates();
	calcFishingMortality();
	populationDynamics();
	observationModels();

	// overwrite existing data with simulated values.
	ct   = value(elem_prod( chat,exp(sigma_delta*rdelta) )) ;
	cpue = value(elem_prod( exp(lnq)*bt,exp(sigma_epsilon*repsilon) ));
	wt   = value(what+sigma_nu*rnu);
	psi  = value(tau)*rpsi;
	
	




FUNCTION initialStates
	bo   = mfexp(log_bo);
	reck = mfexp(log_reck) + 1.001;
	m    = mfexp(log_m);
	sigma_nu      = 1.0 / mfexp(log_sigma_nu);
	sigma_epsilon = 1.0 / mfexp(log_sigma_epsilon);
	sigma_delta   = 1.0 / mfexp(log_sigma_delta);
	tau           = 1.0 / mfexp(log_tau);
	//wk   = mfexp(log_wk);

	dvariable s    = exp(-m);
	dvariable wbar = (s*alpha+wk*(1.-s))/(1-rho*s);
	no   = bo/wbar;
	ro   = no*(1.-s);
	so   = reck*ro/bo;
	beta = (reck-1.)/bo;


FUNCTION calcFishingMortality
	ft = mfexp(log_fbar + fdev);

FUNCTION populationDynamics
	// add option for starting off an fished state.
	dvariable s  = exp(-m-ft(1));
	dvariable we = (s*alpha+wk*(1.-s))/(1-rho*s);
	dvariable be = -(we*(wk*so-1.)+s*(alpha+rho*we))/(beta*(s*alpha+s*rho*we-we));
	// Be = -(We*(wk*so-1.)+t2*(alpha+rho*We))/(beta*(t2*alpha+t2*rho*We-We));

	bt(1) = be;
	nt(1) = be/we;
	rt(1,agek) = be*(1-s)/we * exp(psi(1,agek));

	int i;
	for( i = 1; i < nyrs; i++ )
	{

		st(i)   = exp(-m-ft(i));
		bt(i+1) = st(i)*(alpha*nt(i)+rho*bt(i)) + wk * rt(i);
		nt(i+1) = st(i)*nt(i) + rt(i);
		if(i+agek <= nyrs)
		{
			rt(i+agek) = so*bt(i)/(1.+beta*bt(i)) * exp(psi(i));
		}
	}
	//sd_bterm = bt(nyrs);
	// COUT(bt);

FUNCTION observationModels
	// average weight.
	what = elem_div(bt,nt);
	nu = wt - what;

	// cpue
	dvar_vector zt = log(cpue) - log(bt);
	lnq = mean(zt);
	epsilon = zt - lnq;

	// random walk in q for cpue.
	//qt(1) = exp(zt(1));
	//dvar_vector fd_zt = first_difference(zt);
	//dvariable zwbar = mean(fd_zt);
	//epsilon(1,nyrs-1) = fd_zt -zwbar;
	//for(int i = 2; i <= nyrs; i++ )
	//{
	//	qt(i) = qt(i-1) * exp(fd_zt(i-1));
	//}

	// COUT(epsilon(1,5));
	// COUT(log(cpue(1,5))-lnq+log(bt(1,5)));


	// predicted catch
	chat  = elem_prod(ft,bt);
	delta = log(ct) - log(chat); 

	//FUNCTION dvariable mydnorm(dvariable& x, dvariable& mu, dvariable &sd)
	//dvariable nloglike;
	//nloglike = log(sd)+0.5*log(2.*M_PI)+0.5*square(x-mu)/(sd*sd);
	//return(nloglike);
FUNCTION calculatePriors
	prior_vec.initialize();

	prior_vec(1) = dlnorm(bo,1.65,0.2);
	//dvariable h  = reck/(4.+reck);
	//prior_vec(2) = dbeta((h-0.2)/0.8,8.0,3.0);
	prior_vec(2) = dlnorm(reck,log(12),0.2);
	prior_vec(3) = dlnorm(m,log(0.2),0.05);

	prior_vec(4) = dgamma(1.0/square(sigma_epsilon),1.01,1.01);
	prior_vec(5) = dgamma(1.0/square(sigma_nu),5.01,0.31);
	prior_vec(6) = dgamma(1.0/square(sigma_delta),1.01,1.01);
	prior_vec(7) = dgamma(1.0/square(tau),1.01,1.01);



FUNCTION calcObjectiveFunction
	dvar_vector lvec(1,4);
	

	lvec(1) = dnorm(nu,sigma_nu);
	lvec(2) = dnorm(epsilon,sigma_epsilon);
	lvec(3) = dnorm(delta,sigma_delta);
	lvec(4) = dnorm(psi,tau);
	f = sum(lvec) + sum(prior_vec);

	dvariable avgF = mean(ft);
	if( last_phase() )
	{
		f+= 0.001*square(log(avgF/0.2));
	}
	else
	{
		f+= 1000.0*square(log(avgF/0.2));
	}

	//dvariable x = 1;
	//dvariable mu = 0;
	//dvariable sd = 1;
	//COUT(mydnorm(x,mu,sd));

	//double xx = 1;
	//double muu = 0;
	//double sdd = 1;
	//COUT(mydnorm(xx,mu,sd));  //works.
	//COUT(mydnorm(x,muu,sdd));

	

REPORT_SECTION
	REPORT(bo);
	REPORT(reck);
	double h = value(reck/(4+reck));
	REPORT(h);
	REPORT(m);
	REPORT(wk);
	REPORT(ro);
	REPORT(no);
	REPORT(tau);
	REPORT(sigma_epsilon);
	REPORT(sigma_nu);
	REPORT(sigma_delta);
	REPORT(year);
	REPORT(bt);
	REPORT(nt);
	REPORT(rt);
	REPORT(wt);
	REPORT(ft);
	REPORT(fdev);
	REPORT(what);
	REPORT(epsilon);
	REPORT(nu);
	REPORT(delta);
	REPORT(cpue);
	REPORT(psi);

	dvector yt = value(exp(lnq)*bt);
	REPORT(yt);
	REPORT(ct);
	REPORT(chat);
	REPORT(psi);

	if ( rseed != 0 && last_phase() )
	{
		ofstream ofs("SimPars.rep",ios::app);
		adstring tt = "\t";
		ofs<<rseed<<tt<<bo<<tt<<reck<<tt<<m<<tt<<sigma_epsilon<<tt<<sigma_nu<<endl;
	}

	// report results for mse model
	if ( mseed == 0 && last_phase() )
	{
		ofstream ofs("ddmod.res");
		ofs<<bo<<endl;
		ofs<<bt(nyrs)<<endl;

	}

FUNCTION runMSE
	cout<<"Running Management Strategy Evaluation"<<endl;
	mseVariables   mv;
	mv.bo   = value(bo);
	mv.reck = value(reck);
	mv.m    = value(m);
	mv.psi  = value(psi);
	mv.ft   = value(ft);
	mv.lnq  = value(lnq);
	mv.pyrs = 40;

	mseData  md;
	md.nyrs  = nyrs ;
  md.alpha = alpha;
  md.rho   = rho  ;
  md.wk    = wk   ;
  md.agek  = agek ;
  md.year  = year ;
  md.ct    = ct   ;
  md.cpue  = cpue ;
  md.wt    = wt   ;


	OperatingModel om(md,mv,argc,argv);
	//om.runOM();
	om.runOM(mseed);

FINAL_SECTION
	//system("cp DDmod.rep ./saveRuns/DDmod.rep");
	if(mseed !=0)	runMSE();

GLOBALS_SECTION
	//#include "stats.cxx"
	#include "OperatingModel.h"
	#undef COUT
	#define COUT(object) cout<<#object "\n"<<object<<endl;
	#undef REPORT
	#define REPORT(object) report<<#object "\n"<<object<<endl;
	// #include "mylib.cpp"
