/* Exact Oriola collective-motor core. Balanced coordinates (Q=h/r-Y,Y,Z),
   one tangent in those coordinates; s=min(1,rho)*tau is computational time. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double8 rhs(double8 x, double4 p) {
    double rho=p.s0,r=p.s1,rate=p.s2,H=x.s0+x.s1,v=1.-x.s2;
    return (double8)(-rho*H, r*H*v-x.s1, r*H*x.s1-x.s2,
        -rho*(x.s3+x.s4), r*v*(x.s3+x.s4)-x.s4-r*H*x.s5,
        r*x.s1*x.s3+r*(x.s0+2.*x.s1)*x.s4-x.s5,
        -rho-2.+r*v, 0.)/rate;
}
__kernel void advance(__global double8 *state,__global const double4 *parameters,
    __global double4 *clock,__global double4 *logs,__global double16 *stats,
    __global int *status,__global int2 *counts,int n,int work,
    double transient,double duration,double rtol,double atol,double step_factor,
    int max_steps) {
    int id=get_global_id(0); if(id>=n||status[id])return;
    double8 x=state[id];double4 p=parameters[id],c=clock[id],l=logs[id];
    double16 st=stats[id];int2 cnt=counts[id];double t=c.s0,h=c.s1;
    for(int it=0;it<work&&t<transient+duration;it++) {
        if(cnt.s0+cnt.s1>=max_steps){status[id]=3;break;}
        int win=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=win<0?transient:transient+(win+1)*(duration/4.);
        while(target<=t&&win<3){win++;target=transient+(win+1)*(duration/4.);}
        double H=x.s0+x.s1,v=1.-x.s2;
        double row=fmax(2.*p.s0,fmax(2.*p.s1*fabs(v)+1.+p.s1*fabs(H),
            p.s1*(fabs(x.s1)+fabs(x.s0+2.*x.s1))+1.));
        double cap=step_factor*p.s2*fmin(.05,.5/fmax(1.,row));
        h=fmin(h,fmin(cap,target-t));
        if(h<1e-15||t+h==t){status[id]=2;break;}
        double8 k1=rhs(x,p);
        double8 k2=rhs(x+h*(1./5.*k1),p);
        double8 k3=rhs(x+h*(3./40.*k1+9./40.*k2),p);
        double8 k4=rhs(x+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double8 k5=rhs(x+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double8 k6=rhs(x+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double8 next=x+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double8 k7=rhs(next,p);
        double8 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double err=0.;
        for(int j=0;j<6;j++)err=fmax(err,fabs(error[j])/(atol+rtol*fmax(fabs(x[j]),fabs(next[j]))));
        if(!all(isfinite(next))||!isfinite(err))err=INFINITY;
        if(err<=1.) {
            double divergence=h*(35./384.*k1.s6+500./1113.*k3.s6+125./192.*k4.s6-2187./6784.*k5.s6+11./84.*k6.s6);
            x=next;x.s6=0.;x.s7=0.;t+=h;cnt.s0++;
            if(target-t<1e-13*fmax(1.,target))t=target;
            double norm=sqrt(x.s3*x.s3+x.s4*x.s4+x.s5*x.s5);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(win>=0)l[win]+=log(norm);
            x.s3/=norm;x.s4/=norm;x.s5/=norm;
            double disk=x.s1*x.s1+(1.-x.s2)*(1.-x.s2);
            st.s6=fmax(st.s6,disk);st.s7=fmax(st.s7,fabs(x.s0));
            if(disk>1.+1e-6||fabs(x.s0)>1.+1e-6){status[id]=1;break;}
            if(win>=0){
                st.s0=fmin(st.s0,x.s0);st.s1=fmax(st.s1,x.s0);
                st.s2=fmin(st.s2,x.s1);st.s3=fmax(st.s3,x.s1);
                st.s4=fmin(st.s4,x.s2);st.s5=fmax(st.s5,x.s2);
                st.s8+=divergence;
                if(win==3){st.s9=fmin(st.s9,x.s0+x.s1);st.sa=fmax(st.sa,x.s0+x.s1);
                    st.sb=fmin(st.sb,x.s1);st.sc=fmax(st.sc,x.s1);st.sd=fmin(st.sd,x.s2);st.se=fmax(st.se,x.s2);}
            }
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    c.s0=t;c.s1=h;state[id]=x;clock[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
