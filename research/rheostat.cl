/* Full rheostat oscillator, Li et al2025 DOI10.3390/polym17050617.
   Exact normalization, parameters c,r,g. No asymptotic thermal elimination. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double8 rhs(double8 s,double4 p){
 double x=s.s0,v=s.s1,y=s.s2,c=p.s0,r=p.s1,g=p.s2;
 return (double8)(v,-x-c*v-fabs(v)*v+y,-r*(y+g*x),
   s.s4,-s.s3-(c+2*fabs(v))*s.s4+s.s5,-r*(s.s5+g*s.s3),0,0);
}
__kernel void advance(__global double8 *state, __global const double4 *params,
        __global double4 *clocks, __global double4 *logs, __global double8 *stats,
        __global int *status, __global int2 *counts,
        int n, int attempts, double transient, double duration, double rtol, double atol,
        double step_factor, int max_attempts) {
    int id=get_global_id(0);
    if(id>=n || status[id]) return;
    double8 s=state[id], st=stats[id]; double4 p=params[id]; double4 c=clocks[id], l=logs[id];
    int2 cnt=counts[id]; double t=c.s0, h=c.s1;
    double hmax=step_factor*fmin((double).05,(double).3/sqrt(1+p.s2+p.s1*p.s1));
    for(int j=0;j<attempts && t<transient+duration;j++) {
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        int window=t<transient ? -1 : min(3,(int)floor((t-transient)/(duration/4)));
        double target=window<0 ? transient : transient+(window+1)*(duration/4);
        while(target<=t && window<3){window++;target=transient+(window+1)*(duration/4);}
        h=fmin(h,fmin(hmax,target-t));
        if(h<(double)1e-13 || t+h==t) {status[id]=2;break;}
        double8 k1=rhs(s,p);
        double8 k2=rhs(s+h*((double)(1./5.)*k1),p);
        double8 k3=rhs(s+h*((double)(3./40.)*k1+(double)(9./40.)*k2),p);
        double8 k4=rhs(s+h*((double)(44./45.)*k1-(double)(56./15.)*k2+(double)(32./9.)*k3),p);
        double8 k5=rhs(s+h*((double)(19372./6561.)*k1-(double)(25360./2187.)*k2+(double)(64448./6561.)*k3-(double)(212./729.)*k4),p);
        double8 k6=rhs(s+h*((double)(9017./3168.)*k1-(double)(355./33.)*k2+(double)(46732./5247.)*k3+(double)(49./176.)*k4-(double)(5103./18656.)*k5),p);
        double8 next=s+h*((double)(35./384.)*k1+(double)(500./1113.)*k3+(double)(125./192.)*k4-(double)(2187./6784.)*k5+(double)(11./84.)*k6);
        double8 k7=rhs(next,p);
        double8 error=h*((double)(71./57600.)*k1-(double)(71./16695.)*k3+(double)(71./1920.)*k4-(double)(17253./339200.)*k5+(double)(22./525.)*k6-(double)(1./40.)*k7);
        double8 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        double err=fmax(fmax(scaled.s0,scaled.s1),fmax(fmax(scaled.s2,scaled.s3),fmax(scaled.s4,scaled.s5)));
        if(!isfinite(err) || !isfinite(next.s0) || !isfinite(next.s1) || !isfinite(next.s2))err=INFINITY;
        if(err<=1) {
            s=next;t+=h;cnt.s0++;
            if(target-t < (double)1e-12*fmax((double)1,target))t=target;
            double amp=fmax(fabs(s.s0),fmax(fabs(s.s1),fabs(s.s2)));
            if(amp>1e5){status[id]=1;break;}
            double norm=sqrt(s.s3*s.s3+s.s4*s.s4+s.s5*s.s5);
            if(!(norm>0)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);
            s.s3/=norm;s.s4/=norm;s.s5/=norm;
            if(window>=0){
                st.s0=fmin(st.s0,s.s0);st.s1=fmax(st.s1,s.s0);
                st.s2=fmin(st.s2,s.s2);st.s3=fmax(st.s3,s.s2);
                st.s4=fmax(st.s4,fabs(s.s1));st.s5=fmax(st.s5,amp);
                if(window==3){st.s6=fmin(st.s6,s.s0);st.s7=fmax(st.s7,s.s0);}
            }
        }else cnt.s1++;
        h*=err==0 ? (double)5 : clamp((double).9*pow(err,(double)-.2),(double).2,(double)5);
    }
    c.s0=t;c.s1=h;state[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
