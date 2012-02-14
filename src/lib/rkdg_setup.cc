#include "miniwarpx.h"
#include "rkdg_algo.h"

/**
   Sets up 'rd' with RKDG algorithm specific data. This essentially
   means allocating memory for the 'rkdg_step' routine and setting the
   weights and absicca for the Gaussian quadrature.

*/
void
rkdg_setup(Run_Data &rd)
{

    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;
    int mwave = rd.mwave;
    int sp_order = rd.sp_order;

    rd.ncoeffs = sp_order; // RKDG needs as many storage level as spatial order
    rd.step = rkdg_step; // function to advance solution by dt
    
    // allocate workspace for RKDG algorithm: the assignment
    // operator rellocates memory for the needed arrays.
    
    RKDG_Workspace *rws = new RKDG_Workspace();

    // compute weights and abscissa for Gaussian quadrature
    rws->w = FArray<double>(Range(1,sp_order),0.0);
    rws->x = FArray<double>(Range(1,sp_order),0.0);
    gauleg(sp_order,-1.0,1.0,rws->x,rws->w);

    // compute Legendre polynomials and their derivatives of order
    // 0..sp_order-1 at abscissa locations x[1..sp_order]
    rws->pmx = FArray<double>(Range(1,sp_order, 1,sp_order), 0.0);
    rws->ppmx = FArray<double>(Range(1,sp_order, 1,sp_order), 0.0);

    for(int m=1; m<=sp_order; m++)
        for(int i=1; i<=sp_order; i++)
        {
            rws->pmx(m,i) = legendre_poly(m-1,rws->x(i));
            rws->ppmx(m,i) = legendre_poly_deriv(m-1,rws->x(i));
        }

    // set normalization constants for basis function
    rws->Cconst = FArray<double>(Range(1,sp_order), 0.0);
    for(int m=1; m<=sp_order; m++)
        rws->Cconst(m) = 1./(2.*(m-1)+1);

    rws->s    = FArray<double>(Range(1-mbc,mx+mbc, 1,mwave),0.0);
    rws->wave = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,mwave),0.0);
    rws->amdq = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->apdq = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->f    = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->fl   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->fr   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->ql   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->qr   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->df   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    rws->fedge= FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);

    rws->q1  = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,sp_order),0.0);
    rws->rhs = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,sp_order),0.0);

    rws->wave1 = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,mwave),0.0);
    rws->wave2 = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,mwave),0.0);
    
    rd.work = (void*) rws;

}
