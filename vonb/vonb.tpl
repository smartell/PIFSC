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

	!! linf = max(len);
	!! k    = 0.2;

	objective_function_value f;

	vector lhat(1,nobs);  //pred len
	vector epsilon(1,nobs);

PROCEDURE_SECTION
	lhat = linf*(1.0-exp(-k*age));
	epsilon = len-lhat;
	f = norm2(epsilon);


REPORT_SECTION



  //system("say you model does not fit, ask Dave Fournier for help");

