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

// A struct is sort of like an R-list, you can pass a whole
// struct just like you would pass a list argument.
struct mseData
{
	int nyrs;
	double alpha;
	double rho;
	double wk;
	int agek;

	ivector year;
	dvector ct;
	dvector cpue;
	dvector wt;
};


class OperatingModel  //: public model_data
{
private:
	
	mseVariables m_mv; // model variable for condition and running om.
	mseData m_md;      // model data for conditioning the om.
	int pyrs;

	int m_mseed;
	int nyrs;
	double alpha;
	double rho;
	double wk;
	int agek;

	ivector year;
	dvector ct;
	dvector cpue;
	dvector wt;

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
	dvector refEpsilon;
	dvector refNu;
	dvector refPsi;
	dvector refDelta;



public:
	OperatingModel(mseData _md, mseVariables _mv,int argc,char * argv[]);
	~OperatingModel();
	
	void runOM();
	void runOM(int mseed);
	void initModel();
	void generateRandomVariables();
	void populationModel();
	void getStockStatus();
	void writeDataFile(int &iyr);
	void runAssessment();
	double getTAC();
};

