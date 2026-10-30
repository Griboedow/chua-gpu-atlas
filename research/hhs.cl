/* Exact autonomous single-cell S2: 5 states + 5 tangent components, no PDE terms. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double16 rhs(double16 s,double4 p){
    double m=s.s0,e=s.s1,d=s.s2,y=s.s3,x=s.s4,mu=p.s0,td=p.s1,tl=p.s2;
    double16 f=(double16)(0.);
    f.s0=1.-4.*m+m*m*e;f.s1=3.*m-m*m*e-d*e;f.s2=(.1-d+mu*d*x)/td;
    f.s3=(e-3.8626272360560576-y)/tl;f.s4=y-x;
    f.s8=(-4.+2.*m*e)*s.s8+m*m*s.s9;
    f.s9=(3.-2.*m*e)*s.s8+(-m*m-d)*s.s9-e*s.sa;
    f.sa=(-1.+mu*x)/td*s.sa+mu*d/td*s.sc;
    f.sb=(s.s9-s.sb)/tl;f.sc=s.sb-s.sc;return f;
}
__kernel void advance(__global double16 *states,__global const double4 *parameters,
    __global double4 *clocks,__global double4 *logs,__global double16 *stats,
    __global int *status,__global int2 *counts,
    int n,int attempts,double transient,double duration,double rtol,double atol,
    double max_step,int max_attempts){
    int id=get_global_id(0);if(id>=n||status[id])return;
    double16 s=states[id],st=stats[id];double4 p=parameters[id],c=clocks[id],l=logs[id];
    int2 cnt=counts[id];double t=c.s0,h=c.s1;
    for(int j=0;j<attempts&&t<transient+duration;j++){
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        int window=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=window<0?transient:transient+(window+1)*(duration/4.);
        while(target<=t&&window<3){window++;target=transient+(window+1)*(duration/4.);}
        h=fmin(h,fmin(max_step,target-t));
        if(h<1e-14||t+h==t){status[id]=2;break;}
        double16 k1=rhs(s,p),k2=rhs(s+h*(1./5.*k1),p);
        double16 k3=rhs(s+h*(3./40.*k1+9./40.*k2),p);
        double16 k4=rhs(s+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double16 k5=rhs(s+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double16 k6=rhs(s+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double16 next=s+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double16 k7=rhs(next,p);
        double16 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double16 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        double err=0.;for(int z=0;z<5;z++){err=fmax(err,scaled[z]);err=fmax(err,scaled[8+z]);}
        if(!all(isfinite(next))||!all(isfinite(scaled)))err=INFINITY;
        if(err<=1.){
            s=next;t+=h;cnt.s0++;if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=0.;for(int z=0;z<5;z++)norm+=s[8+z]*s[8+z];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);for(int z=0;z<5;z++)s[8+z]/=norm;
            for(int z=0;z<5;z++)if(fabs(s[z])>1e6)status[id]=1;
            for(int z=0;z<3;z++)if(s[z]<-1e-8)status[id]=5;
            if(status[id])break;
            if(window>=0)for(int z=0;z<5;z++){st[2*z]=fmin(st[2*z],s[z]);st[2*z+1]=fmax(st[2*z+1],s[z]);}
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    c.s0=t;c.s1=h;states[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
