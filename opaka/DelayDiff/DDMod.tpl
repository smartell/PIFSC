DATA_SECTION
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
		ct   = column(data,2)/1.0e6;
		cpue = column(data,3);
		wt   = column(data,4);
	END_CALCS

INITIALIZATION_SECTION
	log_bo   3.65;
	log_reck 2.48;
	log_m   -1.60;
	log_sigma_epsilon  -1.60;
	log_sigma_nu       -1.60;
	log_fbar           -1.60;



PARAMETER_SECTION
	init_number log_bo;
	init_number log_reck;
	init_number log_m;
	init_number log_sigma_epsilon;
	init_number log_sigma_nu;
	init_number log_fbar;
	init_bounded_dev_vector fdev(1,nyrs,-5,5);

	objective_function_value f;
	number bo;
	number ro;
	number no;
	number reck;
	number m;
	number sigma_epsilon;
	number sigma_nu;

	number so;
	number beta;
	number lnq;

	vector ft(1,nyrs);
	vector st(1,nyrs);
	vector bt(1,nyrs);
	vector nt(1,nyrs);
	vector rt(1,nyrs);
	vector what(1,nyrs);
	vector chat(1,nyrs);
	vector nu(1,nyrs);
	vector delta(1,nyrs);
	vector epsilon(1,nyrs);

PROCEDURE_SECTION

	initialStates();
	calcFishingMortality();
	populationDynamics();
	observationModels();
	calcObjectiveFunction();



FUNCTION initialStates
	bo   = mfexp(log_bo);
	reck = mfexp(log_reck) + 1.0;
	m    = mfexp(log_m);


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

	bt(1) = bo;
	nt(1) = no;
	rt(1,agek) = ro;

	int i;
	for( i = 1; i < nyrs; i++ )
	{
		st(i)   = exp(-m-ft(i));
		bt(i+1) = st(i)*(alpha*nt(i)+rho*bt(i)) + wk * rt(i);
		nt(i+1) = st(i)*nt(i) + rt(i);
		if(i+agek <= nyrs)
		{
			rt(i+agek) = so*bt(i)/(1.+beta*bt(i));
		}
	}
	// COUT(bt);

FUNCTION observationModels
	// average weight.
	what = elem_div(bt,nt);
	nu = wt - what;

	// cpue
	dvar_vector zt = log(cpue) - log(bt);
	lnq = mean(zt);
	epsilon = zt - lnq;
	// COUT(epsilon(1,5));
	// COUT(log(cpue(1,5))-lnq+log(bt(1,5)));


	// predicted catch
	chat  = elem_prod(ft,bt);
	delta = log(ct) - log(chat); 

FUNCTION calcObjectiveFunction
	dvar_vector lvec(1,3);
	sigma_nu      = mfexp(log_sigma_nu);
	sigma_epsilon = mfexp(log_sigma_epsilon);

	lvec(1) = dnorm(nu,sigma_nu);
	lvec(2) = dnorm(epsilon,sigma_epsilon);
	lvec(3) = dnorm(delta,0.2);
	f = sum(lvec);


REPORT_SECTION
	REPORT(bo);
	REPORT(reck);
	REPORT(m);
	REPORT(ro);
	REPORT(no);
	REPORT(year);
	REPORT(bt);
	REPORT(nt);
	REPORT(rt);
	REPORT(epsilon);
	REPORT(nu);
	REPORT(delta);



GLOBALS_SECTION
	#undef COUT
	#define COUT(object) cout<<#object "\n"<<object<<endl;
	#undef REPORT
	#define REPORT(object) report<<#object "\n"<<object<<endl;

