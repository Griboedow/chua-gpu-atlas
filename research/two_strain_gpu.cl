/* Gavish arXiv:2608.13995v1 Eq.2. State coordinates S,ln(I1/mu),
   ln((I2+I12)/mu),I12/(I2+I12),R1; tangent in original compartments
   divided by fixed scales (1,mu,mu,mu,1). No log-coordinate LLE. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double16 rhs(double16 s,double8 p){
    double mu=p.s0,b1=p.s1,b2=p.s2,g1=p.s3,g2=p.s4,g12=p.s5,sigma=p.s6,eta=p.s7;
    double a=1.+(eta-1.)*s.s3,ui=exp(s.s1),uj=exp(s.s2),i1=mu*ui,total=mu*uj,force=b2*total*a;
    return (double16)(mu*(1.-s.s0)-b1*i1*s.s0-force*s.s0,
        b1*s.s0-mu-g1,b2*a*(s.s0+sigma*s.s4)-mu-g2-(g12-g2)*s.s3,
        b2*a*(sigma*s.s4-s.s3*(s.s0+sigma*s.s4))+(g2-g12)*s.s3*(1.-s.s3),
        g1*i1-sigma*force*s.s4-mu*s.s4,
        (-mu-b1*i1-force)*s.s5-b1*s.s0*mu*s.s6-b2*s.s0*mu*(s.s7+eta*s.s8),
        b1*ui*s.s5+(b1*s.s0-mu-g1)*s.s6,
        b2*uj*a*s.s5+(b2*s.s0-mu-g2)*s.s7+b2*eta*s.s0*s.s8,
        sigma*b2*s.s4*(s.s7+eta*s.s8)-(mu+g12)*s.s8+sigma*b2*uj*a*s.s9,
        g1*mu*s.s6-sigma*b2*s.s4*mu*(s.s7+eta*s.s8)-(sigma*force+mu)*s.s9,
        0.,0.,0.,0.,0.,0.);
}
__kernel void advance(__global double16 *states,__global const double8 *parameters,
    __global double4 *clocks,__global double4 *logs,__global double8 *stats,
    __global int *status,__global int2 *counts,int n,int attempts,double transient,
    double duration,double rtol,double atol,double max_step,int max_attempts){
    int id=get_global_id(0);if(id>=n||status[id])return;
    double16 s=states[id];double8 p=parameters[id],st=stats[id];double4 c=clocks[id],l=logs[id];
    int2 cnt=counts[id];double t=c.s0,h=c.s1;
    for(int j=0;j<attempts&&t<transient+duration;j++){
        if(cnt.s0+cnt.s1>=max_attempts){status[id]=3;break;}
        int window=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=window<0?transient:transient+(window+1)*(duration/4.);
        while(target<=t&&window<3){window++;target=transient+(window+1)*(duration/4.);}
        h=fmin(h,fmin(max_step,target-t));
        if(h<1e-13||t+h==t){status[id]=2;break;}
        double16 k1=rhs(s,p),k2=rhs(s+h*(1./5.*k1),p);
        double16 k3=rhs(s+h*(3./40.*k1+9./40.*k2),p);
        double16 k4=rhs(s+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double16 k5=rhs(s+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double16 k6=rhs(s+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double16 next=s+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double16 k7=rhs(next,p);
        double16 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double16 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        double err=0.;for(int z=0;z<10;z++)err=fmax(err,scaled[z]);
        if(!all(isfinite(next))||!all(isfinite(scaled)))err=INFINITY;
        if(err<=1.){
            s=next;t+=h;cnt.s0++;
            if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=0.;for(int z=5;z<10;z++)norm+=s[z]*s[z];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);
            for(int z=5;z<10;z++)s[z]/=norm;
            double pop=s.s0+p.s0*(exp(s.s1)+exp(s.s2))+s.s4;
            if(s.s0< -1e-9||s.s4< -1e-9||s.s3< -1e-9||s.s3>1.+1e-9||pop>1.+1e-9){status[id]=1;break;}
            if(window>=0){
                st.s0=fmin(st.s0,s.s1);st.s1=fmax(st.s1,s.s1);st.s2=fmin(st.s2,s.s2);st.s3=fmax(st.s3,s.s2);
                st.s4=fmin(st.s4,s.s0);st.s5=fmax(st.s5,s.s0);st.s6=fmin(st.s6,1.-pop);st.s7=fmax(st.s7,pop);
            }
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    c.s0=t;c.s1=h;states[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
