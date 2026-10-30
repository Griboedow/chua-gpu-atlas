/* Exact NIFTE state and analytic tangent, adaptive DP5(4), physical seconds. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define M (2*N)
inline void rhs(const double s[M],double f[M],const double A[N*N],const double p[3]){
    for(int i=0;i<N;i++){
        f[i]=0.;f[N+i]=0.;
        for(int j=0;j<N;j++){f[i]+=A[i*N+j]*s[j];f[N+i]+=A[i*N+j]*s[N+j];}
    }
    double th=tanh(p[1]*s[1]);f[0]+=p[2]*p[0]*th;
    f[N]+=p[2]*p[0]*p[1]*(1.-th*th)*s[N+1];
}
__kernel void advance(__global double *states,__global const double *matrices,__global const double *par,
    __global double4 *clocks,__global double4 *logs,__global double *stats,
    __global int *status,__global int2 *counts,int lanes,int attempts,double transient,
    double duration,double rtol,double atol,double max_step,int max_attempts){
    int id=get_global_id(0);if(id>=lanes||status[id])return;
    double s[M],A[N*N],p[3],st[M];
    for(int z=0;z<M;z++){s[z]=states[id*M+z];st[z]=stats[id*M+z];}
    for(int z=0;z<N*N;z++)A[z]=matrices[id*N*N+z];
    for(int z=0;z<3;z++)p[z]=par[id*3+z];
    double4 c=clocks[id],l=logs[id];int2 cnt=counts[id];double t=c.s0,h=c.s1;
    for(int j=0;j<attempts&&t<transient+duration;j++){
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        int window=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=window<0?transient:transient+(window+1)*(duration/4.);
        while(target<=t&&window<3){window++;target=transient+(window+1)*(duration/4.);}
        h=fmin(h,fmin(max_step,target-t));
        if(h<1e-13||t+h==t){status[id]=2;break;}
        double k1[M],k2[M],k3[M],k4[M],k5[M],k6[M],k7[M],tmp[M],next[M];
        rhs(s,k1,A,p);
        for(int z=0;z<M;z++)tmp[z]=s[z]+h*k1[z]/5.;rhs(tmp,k2,A,p);
        for(int z=0;z<M;z++)tmp[z]=s[z]+h*(3./40.*k1[z]+9./40.*k2[z]);rhs(tmp,k3,A,p);
        for(int z=0;z<M;z++)tmp[z]=s[z]+h*(44./45.*k1[z]-56./15.*k2[z]+32./9.*k3[z]);rhs(tmp,k4,A,p);
        for(int z=0;z<M;z++)tmp[z]=s[z]+h*(19372./6561.*k1[z]-25360./2187.*k2[z]+64448./6561.*k3[z]-212./729.*k4[z]);rhs(tmp,k5,A,p);
        for(int z=0;z<M;z++)tmp[z]=s[z]+h*(9017./3168.*k1[z]-355./33.*k2[z]+46732./5247.*k3[z]+49./176.*k4[z]-5103./18656.*k5[z]);rhs(tmp,k6,A,p);
        for(int z=0;z<M;z++)next[z]=s[z]+h*(35./384.*k1[z]+500./1113.*k3[z]+125./192.*k4[z]-2187./6784.*k5[z]+11./84.*k6[z]);rhs(next,k7,A,p);
        double err=0.;
        for(int z=0;z<M;z++){
            double error=h*(71./57600.*k1[z]-71./16695.*k3[z]+71./1920.*k4[z]-17253./339200.*k5[z]+22./525.*k6[z]-1./40.*k7[z]);
            double scaled=fabs(error)/(atol+rtol*fmax(fabs(s[z]),fabs(next[z])));
            err=fmax(err,scaled);if(!isfinite(next[z])||!isfinite(scaled))err=INFINITY;
        }
        if(err<=1.){
            for(int z=0;z<M;z++)s[z]=next[z];t+=h;cnt.s0++;
            if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=0.;for(int z=N;z<M;z++)norm+=s[z]*s[z];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);
            for(int z=N;z<M;z++)s[z]/=norm;
            for(int z=0;z<N;z++){
                if(fabs(s[z])>1e6)status[id]=1; /* Numerical bound, NOT physical vessel boundary. */
                if(window>=0){st[2*z]=fmin(st[2*z],s[z]);st[2*z+1]=fmax(st[2*z+1],s[z]);}
            }
            if(status[id])break;
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    c.s0=t;c.s1=h;clocks[id]=c;logs[id]=l;counts[id]=cnt;
    for(int z=0;z<M;z++){states[id*M+z]=s[z];stats[id*M+z]=st[z];}
}
