DATA_SECTION
	init_int nyrs;
	init_matrix data(1,nyrs,1,4);
	//!! cout<<data<<endl;
	ivector year(1,nyrs);
	vector ct(1,nyrs);
	vector cpue(1,nyrs);
	vector wt(1,nyrs);
	vector mcmcBmsy(1,10000);
	LOC_CALCS
		year = ivector(column(data,1));
		ct   = column(data,2)/1.0e6;
		cpue = column(data,3);
		wt   = column(data,4);
	END_CALCS
	!!CLASS ofstream ofs("MCMC.rep",ios::app);
	!!CLASS ofstream ofs2("MCBT.rep",ios::app);

INITIALIZATION_SECTION
	log_k     3.269017;
	log_r     -1.0;
	log_q     -11.;
	log_sigma -4.;


PARAMETER_SECTION
	init_number log_k;
	init_number log_r;
	init_number log_q;
	init_number log_sigma(2);

	objective_function_value f;
	
	number k;
	number r;
	number q;
	number sigma;
	number fpen;
	vector bt(1,nyrs);
	vector yt(1,nyrs);
	vector epsilon(1,nyrs);

	likeprof_number lp_r;

PROCEDURE_SECTION

	// MAIN ROUTINE
	initParameters();
	productionModel();
	observationModel();
	calcObjectiveFunction();

	//  Calculate the posterior densities for 
	//   FMSY MSY and BMSY.  USING any method.
	if(mceval_phase())
	{
		calcStatusReferencePoints();
	}

FUNCTION calcStatusReferencePoints
	static long nf=0;
	nf++;
	double bmsy = 0.5 * value(k);
	double fmsy = 0.5 * value(r);
	double  msy = fmsy * bmsy;
	double bend = value(bt(nyrs));
	bend = (bend-bmsy)/bmsy;
	if(nf==1)
	{
		ofstream eee("MCBT.rep");
		ofstream fff("MCMC.rep");
		fff<<"Bmsy\t"<<"Fmsy\t"<<"MSY\t"<<"P(Bt<Bmsy)\t"<<endl;
	}
	//ofstream ofs("MCMC.rep",ios::app);
	ofs<<bmsy<<"\t"<<fmsy<<"\t"<<msy<<"\t"<<bend<<endl;

	ofs2<<bt<<endl;

FUNCTION initParameters
	k     = exp(log_k);
	r     = exp(log_r);
	q     = exp(log_q);
	sigma = exp(log_sigma);
	lp_r  = r;

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
	
	prior(1) = dlnorm(k,log(3.269017),0.10);
	prior(2) = dlnorm(r,log(0.2),0.05);
	prior(3) = -log(q);
	prior(4) = dgamma(1.0/square(sigma),1.01,1.01);

	f = dnorm(epsilon,sigma) + sum(prior) + fpen;

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

	//#include "stats.cxx"
	









