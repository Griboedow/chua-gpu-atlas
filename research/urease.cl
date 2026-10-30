/* Exact four-species network in linear balance coordinates; seconds and 1e-4 M.
   L-stable SDIRK2 base with step doubling and Richardson extrapolation. The
   accepted extrapolate is NOT unconditionally A-stable; adaptive error control
   and independent stiff-solver checks are required. The tangent
   is the derivative of the implicit discrete map, not a finite displacement. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define GAMMA 0.2928932188134524756
inline int physical(const double *y){
    return y[0]>=0. && y[3]>0. && y[1]+y[3]>=0. && y[2]-y[3]>=0.;
}
inline void field(const double *y,const double *p,double *f,double *j){
    double s=y[0],d=y[1],z=y[2],h=y[3],km=p[2]/1e-4,e1=p[3]/1e-4,e2=p[4]/1e-4;
    double den=h+h*h/e1+e2,act=h/den,vm=p[8]/1e-4;
    double v=vm*s/(km+s)*act,vs=vm*km/((km+s)*(km+s))*act;
    double vh=vm*s/(km+s)*(e2-h*h/e1)/(den*den);
    double a=p[0]*1e-4,r=p[1],ks=p[5],kh=p[6],k=p[7],se=p[9]/1e-4,he=p[10]/1e-4;
    f[0]=-v+ks*(se-s);f[1]=2*v-k*d+(kh-k)*h-kh*he;
    f[2]=kh*he-k*z+(k-kh)*h;f[3]=-a*(d+h)*h+r*(z-h)+kh*(he-h);
    for(int i=0;i<16;i++)j[i]=0.;
    j[0]=-vs-ks;j[3]=-vh;j[4]=2*vs;j[5]=-k;j[7]=2*vh+kh-k;
    j[10]=-k;j[11]=k-kh;j[13]=-a*h;j[14]=r;j[15]=-a*(d+2*h)-r-kh;
}
inline int linear(const double *jac,double gh,const double *b,double *out){
    double a[20];
    for(int i=0;i<4;i++){for(int j=0;j<4;j++)a[5*i+j]=(i==j?1.:0.)-gh*jac[4*i+j];a[5*i+4]=b[i];}
    for(int k=0;k<4;k++){
        int pivot=k;for(int i=k+1;i<4;i++)if(fabs(a[5*i+k])>fabs(a[5*pivot+k]))pivot=i;
        if(!isfinite(a[5*pivot+k])||fabs(a[5*pivot+k])<1e-30)return 0;
        if(pivot!=k)for(int j=k;j<5;j++){double v=a[5*k+j];a[5*k+j]=a[5*pivot+j];a[5*pivot+j]=v;}
        for(int i=k+1;i<4;i++){double q=a[5*i+k]/a[5*k+k];for(int j=k+1;j<5;j++)a[5*i+j]-=q*a[5*k+j];}
    }
    for(int i=3;i>=0;i--){double v=a[5*i+4];for(int j=i+1;j<4;j++)v-=a[5*i+j]*out[j];out[i]=v/a[5*i+i];if(!isfinite(out[i]))return 0;}
    return 1;
}
inline int stage(const double *base,const double *guess,const double *p,double gh,double rtol,double atol,double *y,double *f,double *jac){
    for(int i=0;i<4;i++)y[i]=guess[i];
    for(int iter=0;iter<12;iter++){
        field(y,p,f,jac);double residual[4],step[4],err=0.;
        for(int i=0;i<4;i++){residual[i]=base[i]+gh*f[i]-y[i];err=fmax(err,fabs(residual[i])/(atol+rtol*fmax(fabs(base[i]),fabs(y[i]))));}
        if(err<.01 && physical(y))return 1;
        if(!linear(jac,gh,residual,step))return 0;
        double factor=1.,trial[4];int okay=0;
        for(int back=0;back<20;back++){
            for(int i=0;i<4;i++)trial[i]=y[i]+factor*step[i];
            if(physical(trial)){okay=1;break;}factor*=.5;
        }
        if(!okay)return 0;
        double scaled=0.;for(int i=0;i<4;i++){y[i]=trial[i];scaled=fmax(scaled,fabs(factor*step[i])/(atol+rtol*fabs(y[i])));}
        if(scaled<.005){field(y,p,f,jac);return physical(y);}
    }
    return 0;
}
inline int sdirk(const double *input,const double *p,double h,double rtol,double atol,double *output){
    double a[4],b[4],fa[4],fb[4],ja[16],jb[16],base[4],va[4],vb[4];
    if(!stage(input,input,p,GAMMA*h,rtol,atol,a,fa,ja))return 0;
    if(!linear(ja,GAMMA*h,input+4,va))return 0;
    for(int i=0;i<4;i++)base[i]=input[i]+(1.-GAMMA)*h*fa[i];
    if(!stage(base,a,p,GAMMA*h,rtol,atol,b,fb,jb))return 0;
    double tang[4];
    for(int i=0;i<4;i++){double v=input[4+i];for(int j=0;j<4;j++)v+=(1.-GAMMA)*h*ja[4*i+j]*va[j];tang[i]=v;}
    if(!linear(jb,GAMMA*h,tang,vb))return 0;
    for(int i=0;i<4;i++){output[i]=b[i];output[i+4]=vb[i];}
    return 1;
}
__kernel void advance(__global double *states,__global const double *params,__global double *clocks,
    __global double *logs,__global double *ranges,__global int *status,__global int *counts,
    int n,int budget,double transient,double duration,double rtol,double atol,double maxstep,int maxattempts){
    int id=get_global_id(0);if(id>=n||status[id])return;
    double s[8],p[11],l[4],range[8];for(int i=0;i<8;i++){s[i]=states[8*id+i];range[i]=ranges[8*id+i];}
    for(int i=0;i<11;i++)p[i]=params[11*id+i];for(int i=0;i<4;i++)l[i]=logs[4*id+i];
    double t=clocks[2*id],h=clocks[2*id+1];int accepted=counts[2*id],rejected=counts[2*id+1];
    for(int attempt=0;attempt<budget && t<transient+duration;attempt++){
        if(accepted+rejected>=maxattempts){status[id]=3;break;}
        int w=t<transient?-1:min(3,(int)floor((t-transient)/(duration/4.)));
        double target=w<0?transient:transient+(w+1)*duration/4.;
        while(target<=t&&w<3){w++;target=transient+(w+1)*duration/4.;}
        h=fmin(h,fmin(maxstep,target-t));if(h<1e-12||t+h==t){status[id]=2;break;}
        double full[8],mid[8],doubled[8],out[8],err=0.;
        int okay=sdirk(s,p,h,rtol,atol,full)&&sdirk(s,p,h*.5,rtol,atol,mid)&&sdirk(mid,p,h*.5,rtol,atol,doubled);
        if(okay){
            for(int i=0;i<8;i++){
                out[i]=(4.*doubled[i]-full[i])/3.;
                err=fmax(err,fabs(doubled[i]-full[i])/(3.*(atol+rtol*fmax(fabs(out[i]),fabs(s[i])))));
                if(!isfinite(out[i]))okay=0;
            }
            if(!physical(out))okay=0;
        }
        if(!okay)err=INFINITY;
        if(err<=1.){
            double norm=0.;for(int i=4;i<8;i++)norm+=out[i]*out[i];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(w>=0)l[w]+=log(norm);
            for(int i=0;i<8;i++)s[i]=i<4?out[i]:out[i]/norm;
            accepted++;t+=h;if(target-t<1e-12*fmax(1.,target))t=target;
            if(w>=0){double species[4]={s[0],s[3],s[1]+s[3],s[2]-s[3]};for(int i=0;i<4;i++){range[2*i]=fmin(range[2*i],species[i]);range[2*i+1]=fmax(range[2*i+1],species[i]);}}
        }else rejected++;
        h*=err==0.?3.:clamp(.85*pow(err,-1./3.),.15,3.);
    }
    for(int i=0;i<8;i++){states[8*id+i]=s[i];ranges[8*id+i]=range[i];}for(int i=0;i<4;i++)logs[4*id+i]=l[i];
    clocks[2*id]=t;clocks[2*id+1]=h;counts[2*id]=accepted;counts[2*id+1]=rejected;
}
