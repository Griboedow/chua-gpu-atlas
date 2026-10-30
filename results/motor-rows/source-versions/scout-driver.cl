/* Exact constrained Alouges N-row field and analytic Jv. N is compiled per
   system size. No change of feedback law and no clipping of accepted states. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#ifndef N
#define N 8
#endif
#define D (3*N-1)
#define K (2*D+1)
#define ST (4*D+4)

inline void field(const double *x, double *f, const double *p) {
    double ms=0.,mds=0.,last=0.,dlast=0.,sumc=0.;
    for(int j=0;j<N;j++){ms+=x[N+j];mds+=x[D+N+j];sumc+=x[j];}
    ms/=N;mds/=N;
    for(int j=0;j<N-1;j++){last-=x[2*N+j];dlast-=x[D+2*N+j];}
    for(int j=0;j<N;j++) {
        double X=j==N-1?last:x[2*N+j];
        double dX=j==N-1?dlast:x[D+2*N+j];
        double v=-p[0]*X+p[1]*(x[N+j]-ms);
        double dv=-p[0]*dX+p[1]*(x[D+N+j]-mds);
        f[j]=(-x[j]-v*x[N+j]+1.)/p[3];
        f[N+j]=(-x[N+j]+v*x[j]+p[2])/p[3];
        f[D+j]=(-x[D+j]-dv*x[N+j]-v*x[D+N+j])/p[3];
        f[D+N+j]=(-x[D+N+j]+dv*x[j]+v*x[D+j])/p[3];
        if(j<N-1){f[2*N+j]=v/p[3];f[D+2*N+j]=dv/p[3];}
    }
    f[2*D]=(-2.*N-p[0]*(N-1)+p[1]*(1.-1./N)*sumc)/p[3];
}

__kernel void advance(__global double *states,__global const double *parameters,
    __global double4 *clocks,__global double4 *logs,__global double *statistics,
    __global int *status,__global int2 *counts,int lanes,int work,
    double transient,double duration,double rtol,double atol,double step_factor,int max_steps) {
    int id=get_global_id(0);if(id>=lanes||status[id])return;
    double x[K],p[5],st[ST],tmp[K],next[K],k1[K],k2[K],k3[K],k4[K],k5[K],k6[K],k7[K];
    for(int j=0;j<K;j++)x[j]=states[id*K+j];
    for(int j=0;j<5;j++)p[j]=parameters[id*5+j];
    for(int j=0;j<ST;j++)st[j]=statistics[id*ST+j];
    double4 clock=clocks[id],l=logs[id];int2 cnt=counts[id];double t=clock.s0,h=clock.s1;
    for(int it=0;it<work&&t<transient+duration;it++) {
        if(cnt.s0+cnt.s1>=max_steps){status[id]=3;break;}
        int win=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=win<0?transient:transient+(win+1)*duration/4.;
        while(target<=t&&win<3){win++;target=transient+(win+1)*duration/4.;}
        double ms=0.,last=0.;for(int j=0;j<N;j++)ms+=x[N+j];ms/=N;
        for(int j=0;j<N-1;j++)last-=x[2*N+j];
        double row=p[0]+2.*p[1];
        for(int j=0;j<N;j++) {
            double X=j==N-1?last:x[2*N+j];
            double v=-p[0]*X+p[1]*(x[N+j]-ms);
            double coupling=2.*p[1]+p[0]*(j==N-1?N-1:1);
            row=fmax(row,1.+fabs(v)+coupling*fmax(fabs(x[j]),fabs(x[N+j])));
        }
        h=fmin(h,fmin(step_factor*p[3]*fmin(.05,.5/fmax(1.,row)),target-t));
        if(h<1e-15||t+h==t){status[id]=2;break;}
        field(x,k1,p);
        for(int j=0;j<K;j++)tmp[j]=x[j]+h*(k1[j]/5.);field(tmp,k2,p);
        for(int j=0;j<K;j++)tmp[j]=x[j]+h*(3./40.*k1[j]+9./40.*k2[j]);field(tmp,k3,p);
        for(int j=0;j<K;j++)tmp[j]=x[j]+h*(44./45.*k1[j]-56./15.*k2[j]+32./9.*k3[j]);field(tmp,k4,p);
        for(int j=0;j<K;j++)tmp[j]=x[j]+h*(19372./6561.*k1[j]-25360./2187.*k2[j]+64448./6561.*k3[j]-212./729.*k4[j]);field(tmp,k5,p);
        for(int j=0;j<K;j++)tmp[j]=x[j]+h*(9017./3168.*k1[j]-355./33.*k2[j]+46732./5247.*k3[j]+49./176.*k4[j]-5103./18656.*k5[j]);field(tmp,k6,p);
        for(int j=0;j<K;j++)next[j]=x[j]+h*(35./384.*k1[j]+500./1113.*k3[j]+125./192.*k4[j]-2187./6784.*k5[j]+11./84.*k6[j]);
        field(next,k7,p);double err=0.;int finite=1;
        for(int j=0;j<2*D;j++) {
            double error=h*(71./57600.*k1[j]-71./16695.*k3[j]+71./1920.*k4[j]-17253./339200.*k5[j]+22./525.*k6[j]-k7[j]/40.);
            err=fmax(err,fabs(error)/(atol+rtol*fmax(fabs(x[j]),fabs(next[j]))));
            finite=finite&&isfinite(next[j]);
        }
        if(!finite||!isfinite(err))err=INFINITY;
        if(err<=1.) {
            double divint=next[2*D];
            for(int j=0;j<K;j++)x[j]=next[j];x[2*D]=0.;t+=h;cnt.s0++;
            if(target-t<1e-13*fmax(1.,target))t=target;
            double norm=0.;for(int j=0;j<D;j++)norm+=x[D+j]*x[D+j];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(win>=0)l[win]+=log(norm);
            for(int j=0;j<D;j++)x[D+j]/=norm;
            double disk=0.,mx=0.;last=0.;
            for(int j=0;j<N;j++)disk=fmax(disk,x[j]*x[j]+x[N+j]*x[N+j]);
            for(int j=0;j<N-1;j++){last-=x[2*N+j];mx=fmax(mx,fabs(x[2*N+j]));}mx=fmax(mx,fabs(last));
            st[4*D]=fmax(st[4*D],disk);st[4*D+1]=fmax(st[4*D+1],mx);
            if(disk>(1.+p[2]*p[2])*(1.+1e-6)||mx>p[4]*(1.+1e-6)+1e-8){status[id]=1;break;}
            if(win>=0) {
                st[4*D+2]+=divint;
                for(int j=0;j<D;j++){st[j]=fmin(st[j],x[j]);st[D+j]=fmax(st[D+j],x[j]);
                    if(win==3){st[2*D+j]=fmin(st[2*D+j],x[j]);st[3*D+j]=fmax(st[3*D+j],x[j]);}}
            }
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    for(int j=0;j<K;j++)states[id*K+j]=x[j];
    for(int j=0;j<ST;j++)statistics[id*ST+j]=st[j];
    clock.s0=t;clock.s1=h;clocks[id]=clock;logs[id]=l;counts[id]=cnt;
}
