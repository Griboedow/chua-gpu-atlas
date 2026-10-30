/* Exact quadratic-drag pair, state+tangent, time tau=sigma/scale. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double8 rhs(double8 s,double4 p){
    double m=p.s0,k=p.s1;
    return p.s2*(double8)(s.s1,-(1.+k)*s.s0-k*s.s2+m*(1.-fabs(s.s1))*s.s1,
        s.s3,-k*s.s0-(1.+k)*s.s2+m*(1.-fabs(s.s3))*s.s3,
        s.s5,-(1.+k)*s.s4-k*s.s6+m*(1.-2.*fabs(s.s1))*s.s5,
        s.s7,-k*s.s4-(1.+k)*s.s6+m*(1.-2.*fabs(s.s3))*s.s7);
}
__kernel void advance(__global double8 *states,__global const double4 *parameters,
    __global double4 *clocks,__global double4 *logs,__global double8 *stats,
    __global double4 *modal,__global int *status,__global int2 *counts,
    int n,int attempts,double transient,double duration,double rtol,double atol,
    double max_step_sigma,int max_attempts){
    int id=get_global_id(0);if(id>=n||status[id])return;
    double8 s=states[id],st=stats[id];double4 p=parameters[id],c=clocks[id],l=logs[id],md=modal[id];
    int2 cnt=counts[id];double t=c.s0,h=c.s1;
    for(int j=0;j<attempts&&t<transient+duration;j++){
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        int window=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=window<0?transient:transient+(window+1)*(duration/4.);
        while(target<=t&&window<3){window++;target=transient+(window+1)*(duration/4.);}
        h=fmin(h,fmin(max_step_sigma/p.s2,target-t));
        if(h<1e-14||t+h==t){status[id]=2;break;}
        double8 k1=rhs(s,p),k2=rhs(s+h*(1./5.*k1),p);
        double8 k3=rhs(s+h*(3./40.*k1+9./40.*k2),p);
        double8 k4=rhs(s+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double8 k5=rhs(s+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double8 k6=rhs(s+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double8 next=s+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double8 k7=rhs(next,p);
        double8 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double8 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        double err=0.;for(int z=0;z<8;z++)err=fmax(err,scaled[z]);
        if(!all(isfinite(next))||!all(isfinite(scaled)))err=INFINITY;
        if(err<=1.){
            s=next;t+=h;cnt.s0++;if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=length(s.hi);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);s.hi/=norm;
            for(int z=0;z<4;z++)if(fabs(s[z])>1e6){status[id]=1;break;}
            if(status[id])break;
            if(window>=0){
                for(int z=0;z<4;z++){st[2*z]=fmin(st[2*z],s[z]);st[2*z+1]=fmax(st[2*z+1],s[z]);}
                double sx=(s.s0+s.s2)/2.,sv=(s.s1+s.s3)/2.,ax=(s.s0-s.s2)/2.,av=(s.s1-s.s3)/2.;
                double es=sv*sv+(1.+2.*p.s1)*sx*sx,ea=av*av+ax*ax;
                double fraction=fmin(es,ea)/(es+ea+1e-300);
                md.s0+=h*es;md.s1+=h*ea;md.s2+=h*fraction;md.s3=fmax(md.s3,fraction);
            }
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    c.s0=t;c.s1=h;states[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;modal[id]=md;counts[id]=cnt;
}
