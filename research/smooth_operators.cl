/* Explicit logistic optical-profile extension of Guo et al.'s ideal relays.
   Balanced state (xm, vm/sqrt(k), ym, xp, vp/sqrt(k), yp), then one tangent. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double sigmoid(double x){if(x>=0.)return 1./(1.+exp(-x));double a=exp(x);return a/(1.+a);}
inline double16 rhs(double16 s,double8 p){
    double root=sqrt(p.s0),sm=sigmoid((s.s3-p.s3)/p.s4),sp=sigmoid((p.s2-s.s0)/p.s4);
    return (double16)(root*s.s1,root*(s.s2-s.s0)-p.s1*s.s1,sm-s.s2,
        root*s.s4,root*(s.s5-s.s3)-p.s1*s.s4,sp-s.s5,
        root*s.s7,root*(s.s8-s.s6)-p.s1*s.s7,sm*(1.-sm)/p.s4*s.s9-s.s8,
        root*s.sa,root*(s.sb-s.s9)-p.s1*s.sa,-sp*(1.-sp)/p.s4*s.s6-s.sb,0.,0.,0.,0.);
}
__kernel void advance(__global double16 *states,__global const double8 *parameters,
    __global double4 *clocks,__global double4 *logs,__global double8 *stats,
    __global int *status,__global int2 *counts,int n,int attempts,double transient,
    double duration,double rtol,double atol,double step_factor,int max_attempts){
    int id=get_global_id(0);if(id>=n||status[id])return;
    double16 s=states[id];double8 p=parameters[id],st=stats[id];double4 c=clocks[id],l=logs[id];
    int2 cnt=counts[id];double t=c.s0,h=c.s1,root=sqrt(p.s0);
    double base=step_factor*fmin(.05,.3/root);
    for(int j=0;j<attempts&&t<transient+duration;j++){
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        int window=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=window<0?transient:transient+(window+1)*(duration/4.);
        while(target<=t&&window<3){window++;target=transient+(window+1)*(duration/4.);}
        /* Extra resolution near optical transitions, also when velocity is small.
           This cap is a numerical safeguard, not a rigorous no-missed-root theorem. */
        double dm=fabs(s.s0-p.s2)+p.s4,dp=fabs(s.s3-p.s3)+p.s4;
        double optical=step_factor*.4*fmin(dm/(root*(fabs(s.s1)+sqrt(dm))),dp/(root*(fabs(s.s4)+sqrt(dp))));
        h=fmin(h,fmin(fmin(base,optical),target-t));
        if(h<1e-13||t+h==t){status[id]=2;break;}
        double16 k1=rhs(s,p);
        double16 k2=rhs(s+h*(1./5.*k1),p);
        double16 k3=rhs(s+h*(3./40.*k1+9./40.*k2),p);
        double16 k4=rhs(s+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double16 k5=rhs(s+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double16 k6=rhs(s+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double16 next=s+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double16 k7=rhs(next,p);
        double16 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double16 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        double err=0.;for(int z=0;z<12;z++)err=fmax(err,scaled[z]);
        if(!all(isfinite(next))||!all(isfinite(scaled)))err=INFINITY;
        if(err<=1.){
            s=next;t+=h;cnt.s0++;
            if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=0.;for(int z=6;z<12;z++)norm+=s[z]*s[z];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);
            for(int z=6;z<12;z++)s[z]/=norm;
            if(s.s2 < -1e-7||s.s2>1.+1e-7||s.s5 < -1e-7||s.s5>1.+1e-7){status[id]=1;break;}
            if(window>=0){
                st.s0=fmin(st.s0,s.s0);st.s1=fmax(st.s1,s.s0);st.s2=fmin(st.s2,s.s3);st.s3=fmax(st.s3,s.s3);
                st.s4=fmin(st.s4,fmin(s.s2,s.s5));st.s5=fmax(st.s5,fmax(s.s2,s.s5));
                if(window==3){st.s6=fmin(st.s6,s.s0);st.s7=fmax(st.s7,s.s0);}
            }
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    c.s0=t;c.s1=h;states[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
