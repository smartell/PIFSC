DATA_SECTION
	init_int nobs;
	init_matrix data(1,nobs,1,2);
	vector age(1,nobs);
	vector len(1,nobs);
	!! age = column(data,1);
	!! len = column(data,2);
	//!! cout<<data<<endl;

	// check to read all data.
	init_int eof;
	!! if(eof != 999){cout<<"Error reading data"<<endl; exit(1);}


PARAMETER_SECTION
	init_number linf;
	init_number k;
	init_number to(1);
	init_number log_cv(1);

	!! linf = max(len);
	!! k    = 0.2;
	!! to = -0.5;
	!! log_cv = log(0.1);

	objective_function_value f;

	vector lhat(1,nobs);  //pred len
	vector epsilon(1,nobs);
	vector sd(1,nobs);
	sdreport_number sdla1;

PROCEDURE_SECTION
	lhat = linf*(1.0-exp(-k*(age-to)));
	sd   = exp(log_cv) * lhat;
	sdla1 = linf*(1.0-exp(-k*(1-to)));
	epsilon = len-lhat;
	

	f = dnorm(epsilon,sd);

REPORT_SECTION
	report<<"Linf\n"<<linf<<endl;

	REPORT(linf);


  
GLOBALS_SECTION
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;




















