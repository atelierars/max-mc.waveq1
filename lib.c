#include"ext.h"            // standard Max include, always required (except in Jitter)
#include"ext_obex.h"        // required for "new" style objects
#include"z_dsp.h"            // required for MSP objects
#include"commonsyms.h"
#include<Accelerate/Accelerate.h>
#include<simd/simd.h>
typedef struct {
    t_pxobject const super;
	uintptr_t const inputs[2];
	uintptr_t const length;
	uintptr_t const elapse;
	double * const memory;
	double const damping;
	simd_float2 const boundary;
} t_waveeq1d;
C74_HIDDEN static t_class const * class = NULL;
C74_HIDDEN void*new(t_symbol const * const symbol, long const argc, t_atom const * const argv) {
	t_waveeq1d * const this = (t_waveeq1d * const)object_alloc((t_class*const)class);
    if ( this ) {
		// attr
		attr_args_process(this, argc, (t_atom*const)argv);
		*(double const**const)(&this->memory) = (double*const)sysmem_newptr(0);
		*(double*const)(&this->damping) = 0;
		*(simd_float2*const)(&this->boundary) = simd_make_float2(0, 0);
		
		// init
        z_dsp_setup((t_pxobject*const)&this->super, 2);
		((t_pxobject*const)&this->super)->z_misc |= Z_MC_INLETS;//|Z_NO_INPLACE;
        outlet_new(this, "multichannelsignal");
    }
    return this;
}
C74_HIDDEN void del(t_waveeq1d const * const this) {
    dsp_free((t_pxobject*const)&this->super);
	sysmem_freeptr(this->memory);
}
C74_HIDDEN bool inputchanged(t_waveeq1d const * const this, long const index, long const count) {
	switch (index) {
		case 0:
		case 1:
			return*(uintptr_t*const)(this->inputs + index) = count;
		default:
			assert(false);
			return false;
	}
}
C74_HIDDEN long multichanneloutputs(t_waveeq1d const * const this) {
	return this->length;
}
C74_HIDDEN void perfroutine64(t_waveeq1d * const this, t_object const * const dsp64, double const * const * const ins, long const numins, double * const * const outs, long const numouts, long const framecount, long const flags, void const * const userparam) {
	register uintptr_t const N = this->length;
	register uintptr_t const X = this->inputs[0];
	register uintptr_t const Z = this->inputs[1];
	C74_ASSERT(numouts == N);
	C74_ASSERT(numins == X + Z);
	double * const w = this->memory;
	double * const u[3] = {
		this->memory + 1 * N,
		this->memory + 2 * N,
		this->memory + 3 * N,
	};
	double * const a[3] = {
		this->memory + 4 * N,
		this->memory + 5 * N,
		this->memory + 6 * N,
	};
	double * const m[3] = {
		this->memory + 7 * N,
		this->memory + 8 * N,
		this->memory + 9 * N,
	};
	double const r = pow(10, 0.05 * this->damping / (uintptr_t const)userparam);
	register simd_float2 const b = 1 - this->boundary;
	register uintptr_t tau = this->elapse % 3 + 3;
	double const * const * const x = ins;
	double const * const * const z = ins + X;
	for ( register uintptr_t t = 0, T = framecount ; t < T ; ++ t, ++ tau ) {
		register uintptr_t const i = ( tau + 2 ) % 3;
		register uintptr_t const j = ( tau + 3 ) % 3;
		register uintptr_t const k = ( tau + 4 ) % 3;
		
		// a[k] <- in1
		// w    <- in2
		for ( register uintptr_t n = 0 ; n < N ; ++ n ) {
			w[n] = z[n % Z][t];
			a[k][n] = isnormal(w[n]) ? x[n % X][t] / w[n] : 0;
		}
		
		// m[0] <- weight_{k+1} + weight_{k}
		// m[1] <- weight_{k+1} - weight_{k}
		vDSP_vaddsubD(w + 0, 1,
					  w + 1, 1,
					  m[0], 1,
					  m[1], 1,
					  N - 1);
		
		// q <- ( weight_{k+1} - weight_{k} ) / ( weight_{k+1} + weight_{k} )
		vDSP_vdivD(m[0], 1,
				   m[1], 1,
				   m[2], 1, 
				   N - 1);
		
		for ( register uintptr_t n = 0 ; n < N ; ++ n )
			m[2][n] = isnormal(m[2][n]) ? m[2][n] : 0;
		
		// m[0] <- 1 + m[2]
		// m[1] <- 1 - m[2]
		vDSP_vaddsubD(m[2], 1,
					  (double const[]){1}, 0,
					  m[0], 1,
					  m[1], 1,
					  N - 1);
		
		// p <- ∂u/∂x
		vDSP_vsubD(u[j] + 0, 1,
				   u[j] + 1, 1,
				   m[2], 1,
				   N - 1); // ∂u/∂x
		
		// u[f] <- m_{k+1} * ( ∂u/∂x )_{k+1} - m_{k} ( ∂u/∂x )_{k}
		vDSP_vmmsbD(m[2] + 1, 1, m[0] + 1, 1,
					m[2] + 0, 1, m[1] + 0, 1,
					u[k] + 1, 1,
					N - 2);
		
		// boundaries
//		u[k][  0] = u[j][  1] + b[0] * u[i][  0] - 2 * u[j][  0];
//		u[k][N-1] = b[1] * u[i][N-1] + u[j][N-2] - 2 * u[j][N-1];
		u[k][  0] = ( u[j][  1] - u[j][  0] ) * m[0][0  ] - ( u[j][  0] - b[0] * u[i][  0] );
		u[k][N-1] = ( b[1] * u[i][N-1] - u[j][N-1] ) - ( u[j][N-1] - u[j][N-2] ) * m[1][N-2];
		
		// u[k] <- u[k] + 2 * u[j]
		vDSP_vsmaD(u[j], 1,
				   (double const[1]){2},
				   u[k], 1,
				   u[k], 1,
				   N);
		
		// u[k] <- u[k] - u[i]
		vDSP_vsubD(u[i], 1,
				   u[k], 1,
				   u[k], 1,
				   N);
		
		// u[f] <- u[f] + a[f]
		vDSP_vaddD(a[k], 1, u[k], 1, u[k], 1, N);
		
		// u[f] <- u[f] - a[f]
		vDSP_vsubD(a[i], 1, u[k], 1, u[k], 1, N);
		
		// u <- u * r
		vDSP_vsmulD(u[j], 1, &r, u[j], 1, N);
		vDSP_vsmulD(u[k], 1, &r, u[k], 1, N);
		
		// a <- a * r
		vDSP_vsmulD(a[j], 1, &r, a[j], 1, N);
		vDSP_vsmulD(a[k], 1, &r, a[k], 1, N);
		
		for ( register uintptr_t n = 0 ; n < N ; ++ n ) {
			outs[n][t] = w[n] * u[k][n];
		}
		
	}
	*(uintptr_t*const)&this->elapse = tau % 3;
}
C74_HIDDEN void dsp64(t_waveeq1d const * const this, t_object const * const dsp64, short const * const count, double const samplerate, long const maxvectorsize, long const flags) {
	*(double const**const)&this->memory = (double*const)sysmem_resizeptrclear(this->memory, 10 * this->length * sizeof(double const));
	*(uintptr_t*const)&this->elapse = 0;
	dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)perfroutine64, 0, (uintptr_t const)samplerate);
}
C74_HIDDEN void damping(t_waveeq1d const * const this, double const value) {
	*(double*const)&this->damping = value;
}
C74_HIDDEN void boundary(t_waveeq1d const * const this, double const l, double const r) {
	*(simd_float2*const)&this->boundary = simd_make_float2(l, r);
}
C74_EXPORT void ext_main(void * const _) {
    if ( !class ) {
        t_class * const obj = class_new("mc.waveq1~", (method const)new, (method const)del, sizeof(t_waveeq1d), NULL, A_GIMME, 0);
        class_addmethod(obj, (method const)inputchanged, "inputchanged", A_CANT, 0);
        class_addmethod(obj, (method const)multichanneloutputs, "multichanneloutputs", A_CANT, 0);
        class_addmethod(obj, (method const)dsp64, "dsp64", A_CANT, 0);
		class_addmethod(obj, (method const)damping, "damping", A_FLOAT, 0);
		class_addmethod(obj, (method const)boundary, "boundary", A_FLOAT, A_FLOAT, 0);
		class_addattr(obj, attr_offset_new("chans", gensym("long"), 0, (method const)0L, (method const)0L, offsetof(t_waveeq1d, length)));
        class_dspinit(obj);
        class_register(CLASS_BOX, obj);
        class = obj;
    }
}
