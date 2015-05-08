#include <complex.h>
void fill_model(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq,int nfreq, int n, void *mat);
void fill_model_fast(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in);
void fill_model_double(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in);
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in);
double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n);
double calculate_chisq_ampfit(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n);
double calculate_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n);
void calculate_many_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, float *amps, float *RMs, float *phis, int nlike, double *chisq_vec);
double calculate_chisq_qu_rmpow(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n);
double calculate_chisq_cached(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in);
void fill_model_fast_general(float dt, float fwhm, float tau, float scat, float alpha, float t0, float DM, float DM_ind, float *freq, int nfreq, int n, void *mat_in);
double calculate_chisq_qu_cached(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, void *myscratch);
double calculate_chisq_cached_tau(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in);
double calculate_chisq_cached_dmpow(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in);
double calculate_chisq_qu_rmpow_cached(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, void *myscratch);
double calculate_chisq_linfit_cached(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in);
double calculate_chisq_scatfit_cached(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in);
