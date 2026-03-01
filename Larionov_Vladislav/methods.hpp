extern double F(double);

#ifdef __NEWTON
extern double F1(double);
#endif

double Round(double X, double Delta);

double BISECT(double Left, double Right, double Eps, int &N);

double ITER(double X0, double Eps, int &N);

double HORDA(double Left, double Right, double Eps, int &N);

double NEWTON(double X, double Eps, int &N);
