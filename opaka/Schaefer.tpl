DATA_SECTION
	init_int nyrs;
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
	log_k     2.;
	log_r     -2.5;
	log_q     -5.;
	log_sigma -4.6;


PARAMETER_SECTION
	init_number log_k;
	init_number log_r;
	init_number log_q;
	init_number log_sigma;

	objective_function_value f;
	
	number k;
	number r;
	number q;
	number sigma;
	number fpen;
	vector bt(1,nyrs);
	vector yt(1,nyrs);
	vector epsilon(1,nyrs);

PROCEDURE_SECTION

	// MAIN ROUTINE
	initParameters();
	productionModel();
	observationModel();
	calcObjectiveFunction();

FUNCTION initParameters
	k     = exp(log_k);
	r     = exp(log_r);
	q     = exp(log_q);
	sigma = exp(log_sigma);

FUNCTION productionModel
	int i;
	fpen = 0;
	for( i = 1; i < nyrs; i++ )
	{
		if(i == 1) bt(i) = k;
		bt(i+1) = posfun(bt(i) + r*bt(i)*(1-bt(i)/k) - ct(i),0.01,fpen);
	}

FUNCTION observationModel
	yt = q * bt;
	epsilon = log(cpue) - log(yt);


FUNCTION calcObjectiveFunction
	dvar_vector prior(1,4);
	prior.initialize();
	
	prior(1) = dlnorm(k,log(12),0.10);
	prior(2) = dlnorm(r,log(0.2),0.05);
	prior(3) = -log(q);
	prior(4) = dgamma(1.0/square(sigma),1.01,1.01);

	f = dnorm(epsilon,sigma) + sum(prior) + 1000.*fpen;

REPORT_SECTION
	REPORT(k);
	REPORT(r);
	REPORT(q);
	REPORT(sigma);
	REPORT(bt);
	REPORT(yt);
	REPORT(ct);
	REPORT(cpue);
	REPORT(year);
	REPORT(epsilon);
	dvector ut = value(elem_div(ct,bt));
	REPORT(ut);

GLOBALS_SECTION
	#undef REPORT 
	#define REPORT(object) report<< #object "\n"<<object <<endl;
