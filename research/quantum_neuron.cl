/* Independent Eq. 10 implementation, Potter et al. (2024). */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double8 rhs(double8 s, double8 p, __global const double8 *table) {
    double q=(p.s6-0.7071067811865475244*s.s3+32.0)/TABLE_DX;
    if(q<0 || q>=TABLE_N-1) return (double8)(NAN);
    int i=(int)floor(q); double u=q-i,u2=u*u,u3=u2*u;
    double8 a=table[i],b=table[i+1];
    double4 f=(2*u3-3*u2+1)*a.lo+(u3-2*u2+u)*TABLE_DX*a.hi
              +(-2*u3+3*u2)*b.lo+(u3-u2)*TABLE_DX*b.hi;
    double4 df=((6*u2-6*u)*a.lo+(3*u2-4*u+1)*TABLE_DX*a.hi
              +(-6*u2+6*u)*b.lo+(3*u2-2*u)*TABLE_DX*b.hi)/TABLE_DX*(-0.7071067811865475244);
    double x=s.s0,y=s.s1,z=s.s2,v=s.s3;
    double g=.5*((1+z)*f.s0+2*x*f.s1+(1-z)*f.s2);
    double gv=.5*((1+z)*df.s0+2*x*df.s1+(1-z)*df.s2);
    double dv=p.s0-(1+p.s1*g)*v;
    double duv=-p.s1*v*f.s1*s.s4-p.s1*v*.5*(f.s0-f.s2)*s.s6+(-1-p.s1*g-p.s1*v*gv)*s.s7;
    return (double8)(p.s3*y+2*dv*z-p.s2*x,-p.s3*x-p.s2*y,
      -2*dv*x-p.s4*p.s2*(z-p.s5),dv,
      p.s3*s.s5+2*duv*z+2*dv*s.s6-p.s2*s.s4,-p.s3*s.s4-p.s2*s.s5,
      -2*duv*x-2*dv*s.s4-p.s4*p.s2*s.s6,duv);
}
__kernel void advance(__global double8 *state, __global const double8 *params,
        __global double4 *clocks, __global double4 *logs, __global double8 *stats,
        __global int *status, __global int2 *counts, __global const double8 *table, int n, int attempts,
        double target, int window, double rtol, double atol) {
    int id=get_global_id(0);
    if(id>=n || status[id]) return;
    double8 s=state[id], st=stats[id]; double8 p=params[id]; double4 c=clocks[id], l=logs[id];
    int2 cnt=counts[id]; double t=c.s0, h=c.s1;
    double hmax=fmin((double).1, (double).5/p.s3);
    for(int j=0;j<attempts && t<target;j++) {
        h=fmin(h,fmin(hmax,target-t));
        if(h<(double)1e-13 || t+h==t) {status[id]=2;break;}
        double8 k1=rhs(s,p,table);
        double8 k2=rhs(s+h*((double)(1./5.)*k1),p,table);
        double8 k3=rhs(s+h*((double)(3./40.)*k1+(double)(9./40.)*k2),p,table);
        double8 k4=rhs(s+h*((double)(44./45.)*k1-(double)(56./15.)*k2+(double)(32./9.)*k3),p,table);
        double8 k5=rhs(s+h*((double)(19372./6561.)*k1-(double)(25360./2187.)*k2+(double)(64448./6561.)*k3-(double)(212./729.)*k4),p,table);
        double8 k6=rhs(s+h*((double)(9017./3168.)*k1-(double)(355./33.)*k2+(double)(46732./5247.)*k3+(double)(49./176.)*k4-(double)(5103./18656.)*k5),p,table);
        double8 next=s+h*((double)(35./384.)*k1+(double)(500./1113.)*k3+(double)(125./192.)*k4-(double)(2187./6784.)*k5+(double)(11./84.)*k6);
        double8 k7=rhs(next,p,table);
        double8 error=h*((double)(71./57600.)*k1-(double)(71./16695.)*k3+(double)(71./1920.)*k4-(double)(17253./339200.)*k5+(double)(22./525.)*k6-(double)(1./40.)*k7);
        double8 scaled=fabs(error)/(atol+rtol*fmax(fabs(s),fabs(next)));
        double err=fmax(fmax(scaled.s0,scaled.s1),fmax(fmax(scaled.s2,scaled.s3),fmax(fmax(scaled.s4,scaled.s5),fmax(scaled.s6,scaled.s7))));
        if(!isfinite(err) || !isfinite(next.s0) || !isfinite(next.s1) || !isfinite(next.s2))err=INFINITY;
        if(err<=1) {
            s=next;t+=h;cnt.s0++;
            if(target-t < (double)1e-12*fmax((double)1,target))t=target;
            double amp=dot(s.lo.xyz,s.lo.xyz);
            if(amp>1.00001 || s.s3 < -1e-6 || s.s3 > p.s0+1e-6){status[id]=3;break;}
            double norm=length(s.hi);
            if(!(norm>0)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);
            s.hi/=norm;
            if(window>=0){
                st.s0=fmin(st.s0,s.s3);st.s1=fmax(st.s1,s.s3);
                st.s2=fmax(st.s2,amp);st.s3=fmax(st.s3,length(k7.lo));st.s4=fmax(st.s4,amp);
                if(window==3){st.s6=fmin(st.s6,s.s3);st.s7=fmax(st.s7,s.s3);}
            }
        }else cnt.s1++;
        h*=err==0 ? (double)5 : clamp((double).9*pow(err,(double)-.2),(double).2,(double)5);
    }
    c.s0=t;c.s1=h;state[id]=s;clocks[id]=c;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
