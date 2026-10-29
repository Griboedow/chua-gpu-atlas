#pragma OPENCL EXTENSION cl_khr_fp64 : enable

inline double3 segment(double3 y, double4 p, int on, double t, int homogeneous) {
    double k=p.x,b=p.y,I=p.z,d=p.w,r=1.+I*on,q=I*on/r;
#ifdef NORMALIZED_PHOTO
    q=(double)on;
#endif
    if(homogeneous){q=0.;d=0.;}
    double omega=sqrt(k-b*b*.25),A=-k*(y.z-q)/(r*r-b*r+k),eq=d-q;
    double c=y.x-eq-A,s=(y.y+r*A+b*c*.5)/omega;
    double er=exp(-r*t),eb=exp(-b*t*.5),co=cos(omega*t),si=sin(omega*t);
    return (double3)(eq+A*er+eb*(c*co+s*si),
       -r*A*er+eb*((-b*c*.5+omega*s)*co+(-b*s*.5-omega*c)*si),q+(y.z-q)*er);
}

inline int crossing(double3 *state,double3 *tangent,double4 p,int on,int subdivisions,
                    double *elapsed,double *minimum_speed,double *maximum_u) {
    double direction=on?1.:-1.;
    if(direction*state->y<=1e-10)return 2;
    double step=6.2831853071795864769/(sqrt(p.x-p.y*p.y*.25)*subdivisions);
    double lo=0.,hi=0.;int found=0;
    for(int i=1;i<=8192;i++){
        hi=i*step;double3 y=segment(*state,p,on,hi,0);
        *maximum_u=fmax(*maximum_u,fabs(y.x));
        if(direction*y.x<0.){found=1;break;}
        lo=hi;
    }
    if(!found)return 1;
    for(int j=0;j<48;j++){
        double mid=(lo+hi)*.5;double u=segment(*state,p,on,mid,0).x;
        if(direction*u>0.)lo=mid;else hi=mid;
    }
    double t=(lo+hi)*.5;double3 y=segment(*state,p,on,t,0);
    if(fabs(y.y)<1e-10)return 2;
    double3 v=segment(*tangent,p,on,t,1);
    double3 f=(double3)(y.y,p.x*(p.w-y.x-y.z)-p.y*y.y,p.z*on*(1.-y.z)-y.z);
#ifdef NORMALIZED_PHOTO
    f.z=(1.+p.z)*on-(1.+p.z*on)*y.z;
#endif
    v-=f*(v.x/y.y);v.x=0.;y.x=0.;
    if(!isfinite(y.y)||!isfinite(y.z)||!isfinite(v.y)||!isfinite(v.z))return 3;
    *state=y;*tangent=v;*elapsed+=t;*minimum_speed=fmin(*minimum_speed,fabs(y.y));return 0;
}

__kernel void advance(__global double *states,__global const double4 *parameters,
                      __global double *logs,__global double *durations,__global int *status,
                      int lanes,int count,int window,int subdivisions){
    int i=get_global_id(0);if(i>=lanes||status[i])return;
    __global double *s=states+6*i;double4 p=parameters[i];
    double3 state=(double3)(0.,s[0],s[1]),tangent=(double3)(0.,s[2],s[3]);
    double minimum_speed=s[4]>0.?s[4]:INFINITY,maximum_u=s[5],total=0.,stretch=0.;
    for(int j=0;j<count;j++){
        double elapsed=0.;int code=crossing(&state,&tangent,p,1,subdivisions,&elapsed,&minimum_speed,&maximum_u);
        if(!code)code=crossing(&state,&tangent,p,0,subdivisions,&elapsed,&minimum_speed,&maximum_u);
        if(code){status[i]=code;break;}
        double norm=hypot(tangent.y,tangent.z);
        if(!(norm>0.)||!isfinite(norm)){status[i]=4;break;}
        tangent/=norm;total+=elapsed;stretch+=log(norm);
    }
    if(window>=0){logs[4*i+window]+=stretch;durations[4*i+window]+=total;}
    s[0]=state.y;s[1]=state.z;s[2]=tangent.y;s[3]=tangent.z;s[4]=minimum_speed;s[5]=maximum_u;
}
