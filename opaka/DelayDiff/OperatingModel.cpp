#include <admodel.h>
#include "DDmod.htp"
#include "OperatingModel.h"



OperatingModel::~OperatingModel(){}

OperatingModel::OperatingModel(mseData _md, mseVariables _mv,int argc,char * argv[])
: model_data(argc,argv), m_mv(_mv), m_md(_md)
{
	pyrs = m_mv.pyrs;
	bo   = m_mv.bo;
	reck = m_mv.reck;
	m    = m_mv.m;
	psi  = m_mv.psi;
	ft   = m_mv.ft;
	lnq  = m_mv.lnq;

	// unpackage the model data struct into variables.
	// some kludge to get things working.
	nyrs  = m_md.nyrs;
	alpha = m_md.alpha;
	rho   = m_md.rho;
	wk    = m_md.wk;
	agek  = m_md.agek;

	year  = m_md.year;
	ct    = m_md.ct; 
	cpue  = m_md.cpue;
	wt    = m_md.wt;




	bt.allocate(1,nyrs+pyrs);
	nt.allocate(1,nyrs+pyrs);
	rt.allocate(1,nyrs+pyrs);
	st.allocate(1,nyrs+pyrs);

	bt.initialize();
	nt.initialize();
	rt.initialize();
	st.initialize();

	refCt.allocate(1,nyrs+pyrs);
	refYt.allocate(1,nyrs+pyrs);
	refWt.allocate(1,nyrs+pyrs);

	refCt.initialize();
	refYt.initialize();
	refWt.initialize();

	cout<<bo<<endl;
	cout<<reck<<endl;
	cout<<m<<endl;
	cout<<"End of constructor\n"<<endl;
}

void OperatingModel::initModel()
{
	double s    = exp(-m);
	double wbar = (s*alpha+wk*(1.-s))/(1-rho*s);
	no   = bo/wbar;
	ro   = no*(1.-s);
	so   = reck*ro/bo;
	beta = (reck-1.)/bo;
}

void OperatingModel::populationModel()
{
	double f;
	double fe = ft(1);
	double s  = exp(-m-fe);
	double we = (s*alpha+wk*(1.-s))/(1-rho*s);
	double be = -(we*(wk*so-1.)+s*(alpha+rho*we))/(beta*(s*alpha+s*rho*we-we));
	

	bt(1) = be;
	nt(1) = be/we;
	rt(1,agek) = be*(1-s)/we * exp(psi(1,agek));

	int i;
	for( i = 1; i < nyrs+pyrs; i++ )
	{
		if(i<=nyrs) // condition model
		{
			f        = ft(i);
			refCt(i) = ct(i);
			refYt(i) = cpue(i);
			refWt(i) = wt(i);
		}
		else        // projection period
		{
			getStockStatus();
			double tac = getTAC();

			// implement fishery 
			if(tac < bt(i))
				f = tac / bt(i);
			else
				f = 0.8;
			cout<<"year "<<i<<" f = "<<f<<endl;
			if (f > 1) exit(1);
			refCt(i) = tac; //add error here.

			// observation model
			refYt(i) = exp(lnq+log(bt(i)));
			// COUT(refYt)

			refWt(i) = bt(i) / nt(i);

			writeDataFile(i);

			runAssessment();
			
		}


		// update population model.
		st(i)   = exp(-m-f);
		bt(i+1) = st(i)*(alpha*nt(i)+rho*bt(i)) + wk * rt(i);
		nt(i+1) = st(i)*nt(i) + rt(i);
		if(i+agek <= nyrs)
		{
			rt(i+agek) = so*bt(i)/(1.+beta*bt(i)) * exp(psi(i));
		}
		if(i+agek > nyrs && i+agek <= nyrs+pyrs)
		{
			rt(i+agek) = so*bt(i)/(1.+beta*bt(i));
		}
	}
	
}

void OperatingModel::runAssessment()
{
	#if defined __APPLE__ || defined __linux
  	// system("./DDmod -ind RunMSE.dat -nox > /dev/null 2>&1");
  	// system("./DDmod -ind RunMSE.dat -nox > /dev/null");
  	system("./DDmod -ind RunMSE.dat -nox -est ");
  #endif

  #if defined _WIN32 || defined _WIN64
	  system("DDmod.exe -ind RunMSE.dat -nox -est ");
  #endif
}

double OperatingModel::getTAC()
{
	double tac = 0.2*estBt;
	return tac;

}

void OperatingModel::getStockStatus()
{
	cifstream ifs("ddmod.res");
	ifs >> estBo;
	ifs >> estBt;
}

void OperatingModel::runOM()
{
	cout<<"RUnning operating Model"<<endl;
	initModel();
	populationModel();

}

void OperatingModel::writeDataFile(int &iyr)
{
	ofstream ofs("MSEdata.dat");
	ofs<<iyr<<endl;
	ofs<<alpha<<endl;
	ofs<<rho<<endl;
	ofs<<wk<<endl;
	ofs<<agek<<endl;
	for(int i = 1; i<= iyr; i++)
	{
		ofs<<min(year)+(i-1)<<"\t";
		ofs<<refCt(i)<<"\t";
		ofs<<refYt(i)<<"\t";
		ofs<<refWt(i)<<"\n";
	}

}