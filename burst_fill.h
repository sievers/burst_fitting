#include <complex.h>
void fill_model(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq,int nfreq, int n, void *mat);
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in);
double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n);
double calculate_chisq_ampfit(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n);
double calculate_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n);
void calculate_many_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, float *amps, float *RMs, float *phis, int nlike, double *chisq_vec);
double calculate_chisq_qu_rmpow(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n);
