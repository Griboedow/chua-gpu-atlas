/* Independent implementation of Parmar et al., discrete Eqs. 1--4.
   Fixed endpoints; no continuum discretization or added inertia. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline void rhs(const double *s,double *f,double4 p){
    const double rho=p.x, ratio=p.y, alpha=p.z, beta=p.w;
    for(int j=0;j<N-1;j++){
        double left=j?s[j-1]:0, right=j<N-2?s[j+1]:0;
        f[j]=rho*(left-2*s[j]+right+s[N-1+j+1]-s[N-1+j]);
        left=j?s[D+j-1]:0;right=j<N-2?s[D+j+1]:0;
        f[D+j]=rho*(left-2*s[D+j]+right+s[D+N-1+j+1]-s[D+N-1+j]);
    }
    for(int i=0;i<N;i++){
        double left=i?s[i-1]:0,right=i<N-1?s[i]:0;
        double dl=i?s[D+i-1]:0,dr=i<N-1?s[D+i]:0;
        double tc=tanh(s[2*N-1+i]),te=tanh(right-left);
        f[N-1+i]=(-s[N-1+i]+alpha*tc)/ratio;
        f[2*N-1+i]=ratio*(-s[2*N-1+i]+beta*te);
        f[D+N-1+i]=(-s[D+N-1+i]+alpha*(1-tc*tc)*s[D+2*N-1+i])/ratio;
        f[D+2*N-1+i]=ratio*(-s[D+2*N-1+i]+beta*(1-te*te)*(dr-dl));
    }
}
__kernel void advance(__global double *state,__global const double4 *parameters,
 __global double4 *logs,__global double4 *domain,__global int *status,
 int count,int steps,double h,int window){
    int id=get_global_id(0);if(id>=count||status[id])return;
    double s[2*D],k[2*D],tmp[2*D],sum[2*D];
    for(int i=0;i<2*D;i++)s[i]=state[id*2*D+i];
    double4 l=logs[id],bounds=domain[id],p=parameters[id];
    for(int step=0;step<steps;step++){
        rhs(s,k,p);for(int i=0;i<2*D;i++){sum[i]=k[i];tmp[i]=s[i]+h*k[i]/2;}
        rhs(tmp,k,p);for(int i=0;i<2*D;i++){sum[i]+=2*k[i];tmp[i]=s[i]+h*k[i]/2;}
        rhs(tmp,k,p);for(int i=0;i<2*D;i++){sum[i]+=2*k[i];tmp[i]=s[i]+h*k[i];}
        rhs(tmp,k,p);for(int i=0;i<2*D;i++)s[i]+=h*(sum[i]+k[i])/6;
        double norm=0;for(int i=0;i<D;i++)norm+=s[D+i]*s[D+i];norm=sqrt(norm);
        if(!isfinite(norm)||norm<=0){status[id]=2;break;}
        for(int i=0;i<D;i++){s[D+i]/=norm;if(!isfinite(s[i]))status[id]=2;}
        if(status[id])break;
        if(window>=0){
            l[window]+=log(norm);
            for(int i=0;i<N;i++){
                double left=i?s[i-1]:0,right=i<N-1?s[i]:0;
                bounds.x=fmin(bounds.x,1+right-left);
                bounds.y=fmin(bounds.y,1-s[N-1+i]);
                bounds.z=fmin(bounds.z,1+s[2*N-1+i]);
            }
            for(int i=0;i<D;i++)bounds.w=fmax(bounds.w,fabs(s[i]));
        }
    }
    for(int i=0;i<2*D;i++)state[id*2*D+i]=s[i];
    logs[id]=l;domain[id]=bounds;
}
