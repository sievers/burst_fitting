#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

//gcc-4.8 -O3 -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so
//gcc-4.8 -O3 -fopenmp -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so
void fill_model(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
{



  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;
  float DM0=4150;
  float freq_ref=800.0;
  float sig=fwhm/sqrt(8*log(2))/dt;
  float pi=3.14159265;
  float *dm_lags=(float *)malloc(sizeof(float)*nfreq);
  float *sigvec=(float *)malloc(sizeof(float)*nn);
  for (int i=0;i<nfreq;i++) {
    dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
    dm_lags[i]=dm_lags[i]/dt;
  }
  float lag_min=dm_lags[0];
  for (int i=1;i<nfreq;i++)
    if (dm_lags[i]<lag_min)
      lag_min=dm_lags[i];
  for (int i=0;i<nfreq;i++)
    dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    float fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  float complex *aa2exp=(float complex *)malloc(nn*sizeof(float complex));
  float complex *naa2exp=(float complex *)malloc(nn*sizeof(float complex));
  for (int i=0;i<nn;i++) {
    float complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    float mylag=dm_lags[j]+t0;
    float fac=freq[j]/freq_ref;
    float myscat=scat*(fac*fac*fac*fac);
    
    float mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);

    for (int i=0;i<nn;i++) {


      float aa1=-myscat;
      //float complex aa2=-2*pi*i/n*I;
      //float complex aa=aa1+aa2;
      //float complex aa=-myscat-2*pi*i/n*I;

      //float complex scatvec=(1-cexp(n*aa))/(1-cexp(aa));
      //float complex scatvec=(1-exp(n*aa1)*cexp(n*aa2))/(1-exp(aa1)*cexp(aa2));
      float complex scatvec=(1-exp(n*aa1)*naa2exp[i])/(1-exp(aa1)*aa2exp[i]);
      scatvec=scatvec/mynorm;
      float complex lagvec=cexp(-2*pi*I*i*mylag/n);

      //if ((i==1)&&(j==1)) {
	//printf("aa is %12.5g %12.5g\n",creal(aa),cimag(aa));
	//printf("scatvec is %12.5g %12.5g\n",creal(scatvec),cimag(scatvec));
	//printf("lagvec is %12.5g %12.5g\n",creal(lagvec),cimag(lagvec));
      //}
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}


/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n)
{

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  
  free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
double calculate_chisq_ampfit(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);

  printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));

  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *vec1=(double *)malloc(sizeof(double)*nfreq);
  double *vec2=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    vec1[j]=0;
    vec2[j]=0;
    for (int i=imin;i<nn;i++) {
      float r1=creal(mat[j*nn+i]);
      float i1=cimag(mat[j*nn+i]);
      float r2=creal(dat[j*nn+i]);
      float i2=cimag(dat[j*nn+i]);
      vec1[j]+=r1*r1+i1*i1;
      vec2[j]+=r1*r2+i1*i2;

    }
    
  }
  double denom=0;
  double num=0;
  for (int j=0;j<nfreq;j++) {
    denom+=vec1[j]*weights[j];
    num+=vec2[j]*weights[j];
  }
  float amp=num/denom;
  *amp_fit=amp;

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  
  free(vec1);
  free(vec2);
  free(all_chisq);
  free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
double calculate_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n)
// float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  
  //find that if RM goes up by 1, phi should go down by .143
  //so, to keep phi constant, should be RM*(lambda^2-.143)
  //implies lambda_ref is .3777,and nu_ref is 793 



#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    float lambda=299.79/freq[j];
    float nu_ref=793;
    float lambda_ref =299.79/nu_ref;
    float fac=lambda_ref*lambda_ref+0.010973; //.010973 comes from looking at some likelihoods
    if (j==0)
      printf("fac is %12.6f\n",fac);
    float amp_q=amp*cos(RM*(lambda*lambda-fac)+phi);
    float amp_u=amp*sin(RM*(lambda*lambda-fac)+phi);
    float mycos=cos(RM*(lambda*lambda-fac)+phi);
    float mysin=sin(RM*(lambda*lambda-fac)+phi);
   
    if (j==0) 
      printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));

    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
#if 0
      double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
      double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
      float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
      float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
      double complex deltq=myq-amp*mat[j*nn+i];
      double complex deltu=myu;
#endif
      all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  free(mat);
  return chisq;
}


/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
void calculate_many_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, float *amps, float *RMs, float *phis, int nlike, double *chisq_vec)
// float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  
  
  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  printf("filling model.\n");
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *ref_chisq=(double *)malloc(sizeof(double)*nfreq);

  #pragma omp parallel for
  for (int i=0;i<nfreq;i++) {
    ref_chisq[i]=0;
    for (int j=imin;j<nn;j++) {
      ref_chisq[i]+=creal(dat_q[i*nn+j])*creal(dat_q[i*nn+j])+cimag(dat_q[i*nn+j])*cimag(dat_q[i*nn+j])+creal(dat_u[i*nn+j])*creal(dat_u[i*nn+j])+cimag(dat_u[i*nn+j])*cimag(dat_u[i*nn+j]);
    }
    ref_chisq[i]*=weights[i];
    if (i==1500)
      printf("ref_chisq is %14.4g\n",ref_chisq[i]);
  
  }
  

  //float lambda_ref=299.79/750.0;
  float lambda_ref=299.79/793.0;
  float fac=lambda_ref*lambda_ref+0.010973;
  printf("lambda_ref is %12.4f %12.5f\n",lambda_ref,fac);
  for (int ii=0;ii<nlike;ii++) {
    
    amp=amps[ii];
    RM=RMs[ii];
    phi=phis[ii];
#pragma omp parallel for
    for (int j=0;j<nfreq;j++) {
      float lambda=299.79/freq[j];
      float amp_q=amp*cos(RM*(lambda*lambda-fac)+phi);
      float amp_u=amp*sin(RM*(lambda*lambda-fac)+phi);
      float mycos=cos(RM*(lambda*lambda-fac)+phi);
      float mysin=sin(RM*(lambda*lambda-fac)+phi);
      
      //if (j==0) 
      //printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));
      
      all_chisq[j]=0;
      for (int i=imin;i<nn;i++)  {
#if 0
	double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
	double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
	float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
	float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
	double complex deltq=myq-amp*mat[j*nn+i];
	double complex deltu=myu;
#endif
	all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
      }
      all_chisq[j]*=weights[j];
    }
    //printf("time here is %12.5f\n",omp_get_wtime()-t1);
    double chisq=0;
    for (int i=0;i<nfreq;i++)
      chisq+=all_chisq[i]-ref_chisq[i];
    chisq_vec[ii]=chisq;
  }
  free(all_chisq);
  free(ref_chisq);
  free(mat);
  return;
}

/*--------------------------------------------------------------------------------*/

void calculate_many_chisq_qu_ampphase_fit(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, float *amps, float *RMs, float *phis, int nlike, double *chisq_vec)
// float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  
  
  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  printf("filling model.\n");
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  

  float lambda_ref=299.79/750.0;
  printf("lambda_ref is %12.4f\n",lambda_ref);
  for (int ii=0;ii<nlike;ii++) {
    
    amp=amps[ii];
    RM=RMs[ii];
    phi=phis[ii];
#pragma omp parallel for
    for (int j=0;j<nfreq;j++) {
      float lambda=299.79/freq[j];
      float amp_q=amp*cos(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      float amp_u=amp*sin(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      float mycos=cos(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      float mysin=sin(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      
      //if (j==0) 
      //printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));
      
      all_chisq[j]=0;
      for (int i=imin;i<nn;i++)  {
#if 0
	double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
	double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
	float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
	float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
	double complex deltq=myq-amp*mat[j*nn+i];
	double complex deltu=myu;
#endif
	all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
      }
      all_chisq[j]*=weights[j];
    }
    //printf("time here is %12.5f\n",omp_get_wtime()-t1);
    double chisq=0;
    for (int i=0;i<nfreq;i++)
      chisq+=all_chisq[i];
    chisq_vec[ii]=chisq;
  }
  free(all_chisq);
  free(mat);
  return;
}
