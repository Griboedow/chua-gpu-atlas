/* Independent deterministic phase/orientation equations, Leyva et al. 2026.
   One work item per independent cluster. Analytic variational equation. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline void rhs(__private const double *s,__private double *out,__global const double4 *pairs,__global const double4 *coef){
    double ct[N],st[N],cc[N],sc[N];
    for(int i=0;i<N;i++){ct[i]=cos(s[i]);st[i]=sin(s[i]);cc[i]=cos(s[N+i]);sc[i]=sin(s[N+i]);}
    for(int i=0;i<AUG;i++)out[i]=0;
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)if(i!=j){
        double4 p=pairs[i*N+j];double w=p.x;
        double co=ct[j]*ct[i]+st[j]*st[i],si=st[j]*ct[i]-ct[j]*st[i];
        double co3=4*co*co*co-3*co,si3=3*si-4*si*si*si;
        double4 a=coef[0],b=coef[1];
        double ft=a.x*co+a.y*si+a.z*co3+a.w*si3;
        double fc=b.x*co+b.y*si+b.z*co3+b.w*si3;
        double dt=-a.x*si+a.y*co-3*a.z*si3+3*a.w*co3;
        double dc=-b.x*si+b.y*co-3*b.z*si3+3*b.w*co3;
        double ca=cc[i]*p.y+sc[i]*p.z,sa=sc[i]*p.y-cc[i]*p.z;
        double cb=cc[j]*p.y+sc[j]*p.z,sb=sc[j]*p.y-cc[j]*p.z;
        double du=s[D+j]-s[D+i],ua=s[D+N+i],ub=s[D+N+j];
        out[i]+=w*ft*ca*cb;out[N+i]+=w*fc*sa*cb;
        out[D+i]+=w*(dt*du*ca*cb-ft*(sa*cb*ua+ca*sb*ub));
        out[D+N+i]+=w*(dc*du*sa*cb+fc*(ca*cb*ua-sa*sb*ub));
    }
    double mean=0,umean=0;for(int i=0;i<N;i++){mean+=out[i]/N;umean+=out[D+i]/N;}
    for(int i=0;i<N;i++){out[i]-=mean;out[D+i]-=umean;}
}
__kernel void advance(__global double *state,__global const double4 *pairs,__global const double4 *coefficients,
 __global double4 *logs,__global double4 *stats,__global int *status,int count,int steps,double h,int window){
    int id=get_global_id(0);if(id>=count||status[id])return;
    double s[AUG],k1[AUG],k2[AUG],k3[AUG],k4[AUG],tmp[AUG];
    for(int i=0;i<AUG;i++)s[i]=state[id*AUG+i];
    double4 l=logs[id],stat=stats[id];__global const double4 *p=pairs+id*N*N;__global const double4 *c=coefficients+id*2;
    for(int step=0;step<steps;step++){
        rhs(s,k1,p,c);for(int i=0;i<AUG;i++)tmp[i]=s[i]+.5*h*k1[i];
        rhs(tmp,k2,p,c);for(int i=0;i<AUG;i++)tmp[i]=s[i]+.5*h*k2[i];
        rhs(tmp,k3,p,c);for(int i=0;i<AUG;i++)tmp[i]=s[i]+h*k3[i];
        rhs(tmp,k4,p,c);for(int i=0;i<AUG;i++)s[i]+=h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
        double norm=0,speed=0,mean=0;for(int i=0;i<N;i++)mean+=s[D+i]/N;
        for(int i=0;i<N;i++)s[D+i]-=mean;
        for(int i=0;i<D;i++){norm+=s[D+i]*s[D+i];speed+=k1[i]*k1[i];}
        norm=sqrt(norm);if(!isfinite(norm)||!(norm>0)){status[id]=2;break;}
        for(int i=0;i<D;i++){s[D+i]/=norm;s[i]=remainder(s[i],6.2831853071795864769);}
        if(window>=0){l[window]+=log(norm);stat[window]=fmax(stat[window],sqrt(speed));}
    }
    for(int i=0;i<AUG;i++)state[id*AUG+i]=s[i];logs[id]=l;stats[id]=stat;
}
