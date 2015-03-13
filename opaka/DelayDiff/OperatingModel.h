#include <admodel.h>
#include "DDmod.htp"
#undef COUT
#define COUT(object) cout<<#object"\n"<<object<<endl;

struct mseVariables
{
	int pyrs;

	double bo;
	double reck;
	double m;
	double lnq;

	dvector psi;
	dvector ft;
};


class OperatingModel: public model_data
{
private:
	
	mseVariables m_mv;
	int pyrs;

	double bo;
	double reck;
	double m;
	double no;
	double ro;
	double so;
	double beta;
	double lnq;

	double estBo;
	double estBt;

	dvector bt;
	dvector nt;
	dvector rt;
	dvector ft;
	dvector st;
	dvector psi;

	dvector refCt;
	dvector refYt;
	dvector refWt;

public:
	OperatingModel(mseVariables _mv,int argc,char * argv[]);
	~OperatingModel();
	
	void runOM();
	void initModel();
	void populationModel();
	void getStockStatus();
	void writeDataFile(int &iyr);
	void runAssessment();
	double getTAC();
};

