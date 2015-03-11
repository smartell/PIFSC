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

PARAMETER_SECTION
	init_number log_bo;
	init_number log_reck;
	init_number log_m;
	init_number log_sigma;
	init_number log_fbar;
	init_bounded_dev_vector fdev(1,nyrs,-5,5);


PROCEDURE_SECTION




REPORT_SECTION




