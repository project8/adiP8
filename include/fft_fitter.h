

void fit_fft_to_gaussian(double* pars,double* f,double* out,int N);
void fit_fft_to_lorentian(double* pars,double* f,double* pow,int N, int j);
void fit_fft_to_sinc(double* pars,double* f,double* pow,int nP, double f0, int j);
void fit_fft_to_sinc_2nd(double* pars,double* f,double* pow,int N, int j);
void fit_pow_to_cos(double* pars,double* pow_t,double tstep, int nP, double f0, int j);
void add_noise(double* f,double* pow,int N, int j);
