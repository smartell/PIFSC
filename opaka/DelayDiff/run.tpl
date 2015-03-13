DATA_SECTION

PARAMETER_SECTION
	init_number dum;
	objective_function_value f;

PRELIMINARY_CALCS_SECTION
	runfun();

PROCEDURE_SECTION


FUNCTION runfun
	for(int i=1;i<=20;i++)
	{
		system("./DDmod -sim " + str(i) );
		
	}


