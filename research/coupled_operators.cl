#pragma OPENCL EXTENSION cl_khr_fp64 : enable

inline double3 segment3(double3 y,double k,double b,double q,double t){
    double omega=sqrt(k-b*b*.25),a=k*(y.z-q)/(1.-b+k),c=y.x-q-a;
    double s=(y.y+a+b*c*.5)/omega,eh=exp(-t),eb=exp(-b*t*.5),co=cos(omega*t),si=sin(omega*t);
    return (double3)(q+a*eh+eb*(c*co+s*si),
        -a*eh+eb*((-b*c*.5+omega*s)*co+(-b*s*.5-omega*c)*si),q+(y.z-q)*eh);
}

inline int event(double3 *ym,double3 *yp,double3 *tm,double3 *tp,
    int *mode_m,int *mode_p,double4 p,int subdivisions,double *elapsed,
    double *min_speed,double *min_time,double *max_x,double *min_heat,double *max_heat,int *which){
    double dt=fmin(.1,6.283185307179586/(sqrt(p.x-p.y*p.y*.25)*subdivisions));
    double sm=1.-2.*(*mode_p),sp=2.*(*mode_m)-1.;
    double lo=0.,hi=0.,root_m=INFINITY,root_p=INFINITY;int found=0;
    for(int i=1;i<=16384;i++){
        hi=i*dt;
        double3 m=segment3(*ym,p.x,p.y,*mode_m,hi),v=segment3(*yp,p.x,p.y,*mode_p,hi);
        *max_x=fmax(*max_x,fmax(fabs(m.x),fabs(v.x)));
        int cm=sm*(m.x-p.z)<0.,cp=sp*(v.x-p.w)<0.;
        if(cm||cp){
            if(cm){double l=lo,h=hi;for(int j=0;j<42;j++){double mid=(l+h)*.5;
                if(sm*(segment3(*ym,p.x,p.y,*mode_m,mid).x-p.z)>0.)l=mid;else h=mid;}root_m=(l+h)*.5;}
            if(cp){double l=lo,h=hi;for(int j=0;j<42;j++){double mid=(l+h)*.5;
                if(sp*(segment3(*yp,p.x,p.y,*mode_p,mid).x-p.w)>0.)l=mid;else h=mid;}root_p=(l+h)*.5;}
            found=1;break;
        }
        lo=hi;
    }
    if(!found)return 1;
    if(fabs(root_m-root_p)<1e-10)return 5;
    double t=fmin(root_m,root_p);int j=root_p<root_m;
    double3 m=segment3(*ym,p.x,p.y,*mode_m,t),v=segment3(*yp,p.x,p.y,*mode_p,t);
    double3 a=segment3(*tm,p.x,p.y,0.,t),c=segment3(*tp,p.x,p.y,0.,t);
    double speed=j?v.y:m.y;if(fabs(speed)<1e-9||t<1e-11)return 2;
    if(j){v.x=p.w;int next=1-*mode_m;a.z+=(next-*mode_m)*c.x/speed;*mode_m=next;}
    else {m.x=p.z;int next=1-*mode_p;c.z+=(next-*mode_p)*a.x/speed;*mode_p=next;}
    if(!all(isfinite(m))||!all(isfinite(v))||!all(isfinite(a))||!all(isfinite(c)))return 3;
    *ym=m;*yp=v;*tm=a;*tp=c;*elapsed=t;*which=j;
    *min_speed=fmin(*min_speed,fabs(speed));*min_time=fmin(*min_time,t);
    *min_heat=fmin(*min_heat,fmin(m.z,v.z));*max_heat=fmax(*max_heat,fmax(m.z,v.z));return 0;
}

__kernel void advance(__global double *states,__global const double4 *parameters,
    __global double *logs,__global double *durations,__global int *status,
    int lanes,int count,int window,int subdivisions){
    int i=get_global_id(0);if(i>=lanes||status[i])return;
    __global double *s=states+20*i;double4 p=parameters[i];
    double3 m=vload3(0,s),v=vload3(1,s),a=vload3(2,s),c=vload3(3,s);
    int mm=(int)s[12],mp=(int)s[13];double min_speed=s[14],min_time=s[15],max_x=s[16],min_heat=s[17],max_heat=s[18],count_m=s[19];
    double total=0.,stretch=0.;
    for(int j=0;j<count;j++){
        double elapsed=0.;int which=0;
        int code=event(&m,&v,&a,&c,&mm,&mp,p,subdivisions,&elapsed,&min_speed,&min_time,&max_x,&min_heat,&max_heat,&which);
        if(code){status[i]=code;break;}
        double norm=sqrt(a.x*a.x+a.y*a.y/p.x+a.z*a.z+c.x*c.x+c.y*c.y/p.x+c.z*c.z);
        if(!(norm>0.)||!isfinite(norm)){status[i]=4;break;}
        a/=norm;c/=norm;stretch+=log(norm);total+=elapsed;count_m+=which==0;
    }
    if(window>=0){logs[4*i+window]+=stretch;durations[4*i+window]+=total;}
    vstore3(m,0,s);vstore3(v,1,s);vstore3(a,2,s);vstore3(c,3,s);
    s[12]=mm;s[13]=mp;s[14]=min_speed;s[15]=min_time;s[16]=max_x;s[17]=min_heat;s[18]=max_heat;s[19]=count_m;
}
