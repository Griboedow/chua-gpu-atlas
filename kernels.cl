/* Each work item integrates one parameter tuple and one initial condition. */
#ifdef USE_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double real;
typedef double4 real4;
typedef double8 real8;
#else
typedef float real;
typedef float4 real4;
typedef float8 real8;
#endif

inline real8 rhs(real8 s, real4 p) {
    real x = s.s0, y = s.s1, z = s.s2;
    real h = p.s3*x + ((real)0.5)*(p.s2-p.s3)*(fabs(x+1)-fabs(x-1));
    real slope = fabs(x) < 1 ? p.s2 : p.s3;
    return (real8)(p.s0*(y-x-h), x-y+z, -p.s1*y,
                   p.s0*(s.s4-(1+slope)*s.s3),
                   s.s3-s.s4+s.s5, -p.s1*s.s4, 0, 0);
}

inline real8 rk4(real8 s, real4 p, real dt) {
    real8 a = rhs(s,p);
    real8 b = rhs(s+dt*((real)0.5)*a,p);
    real8 c = rhs(s+dt*((real)0.5)*b,p);
    real8 d = rhs(s+dt*c,p);
    return s + (dt/6)*(a+2*b+2*c+d);
}

__kernel void advance(__global real8 *state, __global const real4 *params,
                      __global real4 *log_sums, __global real8 *stats,
                      __global int *status, int n, int steps, real dt,
                      int window, int renorm, real escape) {
    int id = get_global_id(0);
    if (id >= n || status[id]) return;
    real8 s = state[id], st = stats[id];
    real4 p = params[id], logs = log_sums[id];
    real acc = 0;
    for (int k=0; k<steps; ++k) {
        s = rk4(s,p,dt);
        real amp = fmax(fabs(s.s0),fmax(fabs(s.s1),fabs(s.s2)));
        if (!isfinite(s.s0) || !isfinite(s.s1) || !isfinite(s.s2) || amp>escape) {
            status[id] = 1;
            break;
        }
        if ((k+1)%renorm==0) {
            real norm = sqrt(s.s3*s.s3+s.s4*s.s4+s.s5*s.s5);
            if (!(norm>0) || !isfinite(norm)) {status[id]=2; break;}
            if (window>=0) acc += log(norm);
            s.s3 /= norm; s.s4 /= norm; s.s5 /= norm;
        }
        if (window>=0) {
            st.s0 = fmin(st.s0,s.s0); st.s1 = fmax(st.s1,s.s0);
            st.s2 += s.s0; st.s3 += s.s0*s.s0;
            st.s4 = fmax(st.s4,amp);
            if (window==3) {st.s6=fmin(st.s6,s.s0); st.s7=fmax(st.s7,s.s0);}
        }
    }
    if (window>=0) logs[window] += acc;
    state[id]=s; log_sums[id]=logs; stats[id]=st;
}

/* Small diagnostic batches: retain every integration point for CPU comparison. */
__kernel void trace(__global real8 *state, __global const real4 *params,
                    __global real4 *output, int n, int steps, real dt) {
    int id=get_global_id(0);
    if (id>=n) return;
    real8 s=state[id]; real4 p=params[id];
    for(int k=0;k<steps;++k) {
        s=rk4(s,p,dt);
        output[k*n+id]=(real4)(s.s0,s.s1,s.s2,0);
    }
    state[id]=s;
}
