

void fit_fft_to_gaussian(float* pars,double* f,double* out,int N);
void fit_fft_to_lorentian(float* pars,double* f,double* pow,int N, int j);
void fit_fft_to_sinc(float* pars,double* f,double* pow,int nP, double f0, int j);
void fit_fft_to_sinc_2nd(float* pars,double* f,double* pow,int N, int j);
void fit_pow_to_cos(float* pars,double* pow_t,double tstep, int nP, double f0, int j);
void add_noise(double* f,double* pow,int N, int j);
