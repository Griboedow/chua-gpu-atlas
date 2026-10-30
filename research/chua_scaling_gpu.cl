/* Exact coordinate/time scaling of Chua, not a new circuit model.
Parameters alpha,kappa,m0,m1; state x,y,Z=z/alpha, tangent in same coordinates.
DP tableau inherited from adaptive_kernels.cl; targets supplied in original time. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double real;
typedef double4 real4;
typedef double8 real8;
inline real8 rhs(real8 s, real4 p) {
    real x=s.s0, y=s.s1, Z=s.s2;
    real h=p.s3*x+(p.s2-p.s3)*clamp(x,(real)-1,(real)1);
    real slope=fabs(x)<1?p.s2:p.s3;
    return (real8)(y-x-h, Z+(x-y)/p.s0, -p.s1*y,
                   s.s4-(1+slope)*s.s3, s.s5+(s.s3-s.s4)/p.s0, -p.s1*s.s4, 0, 0);
}
/* Dormand-Prince 5(4), adaptive state + tangent integration in FP64.
   Compiled after kernels.cl. Each lane owns an independent integration clock.
   Short bounded submissions keep the display GPU responsive. */
__kernel void adaptive_advance(__global real8 *state, __global const real4 *params,
        __global real4 *clocks, __global real4 *logs, __global real8 *stats,
        __global int *status, __global int2 *counts, int n, int attempts,
        real target, int window, real rtol, real atol, real escape, real max_step, int max_attempts) {
    int id=get_global_id(0);
    if(id>=n || status[id]) return;
    real8 s=state[id], st=stats[id]; real4 p=params[id], c=clocks[id], l=logs[id];
    int2 cnt=counts[id]; real t=c.s0, h=c.s1;
    target *= p.s0;
    real hmax=max_step;
    for(int j=0;j<attempts && t<target;j++) {
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        h=fmin(h,fmin(hmax,target-t));
        if(h<(real)1e-16) {status[id]=2;break;}
        real8 k1=rhs(s,p);
        real8 k2=rhs(s+h*((real)(1./5.)*k1),p);
        real8 k3=rhs(s+h*((real)(3./40.)*k1+(real)(9./40.)*k2),p);
        real8 k4=rhs(s+h*((real)(44./45.)*k1-(real)(56./15.)*k2+(real)(32./9.)*k3),p);
        real8 k5=rhs(s+h*((real)(19372./6561.)*k1-(real)(25360./2187.)*k2+(real)(64448./6561.)*k3-(real)(212./729.)*k4),p);
        real8 k6=rhs(s+h*((real)(9017./3168.)*k1-(real)(355./33.)*k2+(real)(46732./5247.)*k3+(real)(49./176.)*k4-(real)(5103./18656.)*k5),p);
        real8 next=s+h*((real)(35./384.)*k1+(real)(500./1113.)*k3+(real)(125./192.)*k4-(real)(2187./6784.)*k5+(real)(11./84.)*k6);
        real8 k7=rhs(next,p);
        real8 error=h*((real)(71./57600.)*k1-(real)(71./16695.)*k3+(real)(71./1920.)*k4-(real)(17253./339200.)*k5+(real)(22./525.)*k6-(real)(1./40.)*k7);
        real8 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        real err=fmax(fmax(scaled.s0,scaled.s1),fmax(fmax(scaled.s2,scaled.s3),fmax(scaled.s4,scaled.s5)));
        if(!isfinite(err) || !isfinite(next.s0) || !isfinite(next.s1) || !isfinite(next.s2))err=INFINITY;
        if(err<=1) {
            /* Autonomous field: sub-ULP increments in a large absolute clock
               are valid. Retain them with Kahan compensation in c.s2. */
            s=next;
            real increment=h-c.s2, next_t=t+increment;
            c.s2=(next_t-t)-increment;t=next_t;cnt.s0++;
            if(target-t < (real)1e-12*fmax((real)1,target)){t=target;c.s2=0;}
            real amp=fmax(fabs(s.s0),fmax(fabs(s.s1),fabs(s.s2)));
            if(amp>escape){status[id]=1;break;}
            real norm=sqrt(s.s3*s.s3+s.s4*s.s4+s.s5*s.s5);
            if(!(norm>0)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);
            s.s3/=norm;s.s4/=norm;s.s5/=norm;
            if(window>=0){
                st.s0=fmin(st.s0,s.s0);st.s1=fmax(st.s1,s.s0);
                st.s2=fmin(st.s2,s.s1);st.s3=fmax(st.s3,s.s1);
                st.s4=fmin(st.s4,s.s2);st.s5=fmax(st.s5,s.s2);
                if(window==3){st.s6=fmin(st.s6,s.s0);st.s7=fmax(st.s7,s.s0);}
            }
        }else cnt.s1++;
        h*=err==0 ? (real)5 : clamp((real).9*pow(err,(real)-.2),(real).2,(real)5);
    }
    c.s0=t;c.s1=h;state[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
