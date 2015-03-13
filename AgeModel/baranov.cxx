/*
   Routines for solving the Baranov catch equation.
*/                                                 

#include<admodel.h>

	double get_ft(const double& ct, const double& m, const dvector& va, const dvector& ba);
	

/** get_ft
  Solving the baranov catch equation using Newtons method
  Catch is based on weight.
*/	
double get_ft(const double& ct, const double& m, const dvector& va, const dvector& ba)
{
	double ft;
	//initial guess for ft
	ft=ct/(va*(ba*exp(-m/2.)));
	//cout<<"ft<<"<<ft;
	for(int i=1;i<=50;i++)
	{
		dvector f = ft*va;
		dvector z = m+f;
		dvector s = exp(-z);
		dvector o = (1.-s);
		
		dvector t1 = elem_div(f,z);
		dvector t2 = elem_prod(t1,o);
		dvector t3 = elem_prod(o,ba);
		//predicted catch
		double pct = t2*ba;
		
		//derivative of catch wrt ft
		double dct = sum(
			elem_div(t3,z) 
			- elem_div(elem_prod(f,t3),square(z))
			+ elem_prod(elem_prod(t1,s),ba));
		
		ft -= (pct-ct)/dct;  //newton step
		if(fabs(pct-ct)<1.e-4) break; //do not use for dvariables
	}
	//cout<<" ft>>"<<ft<<endl;
	return(ft);
}  
