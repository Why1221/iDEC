
//using namespace std;

double chi2inv(double x, int n);
double LogGamma(double);
double Gamma(double);
double gamcdf(double x,double a);
double gampdf(double x,double a);
double gaminv(double x,double a);
double kf_lgamma(double z);

static double _kf_gammap(double s, double z);
static double _kf_gammaq(double s, double z);
double kf_gammap(double s, double z);
double kf_gammaq(double s, double z);

