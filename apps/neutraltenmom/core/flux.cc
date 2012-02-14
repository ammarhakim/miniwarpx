#include "miniwarpx.h"
#include "tenmoment.h"


// flux function for the ten moment equations
void
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double rhoe,mxe,mye,mze,Pxxe,Pxye,Pxze,Pyye,Pyze,Pzze;
    double ue,ve,we,pxxe,pxye,pxze,pyye,pyze,pzze;

    double rhoi,mxi,myi,mzi,Pxxi,Pxyi,Pxzi,Pyyi,Pyzi,Pzzi;
    double ui,vi,wi,pxxi,pxyi,pxzi,pyyi,pyzi,pzzi;

    double Ex,Ey,Ez,Bx,By,Bz;

    TenMoment_Vars *tmv = (TenMoment_Vars*) rd.mvar;
    // speed of light
    double c0  = tmv -> c0;
    double c02; 

    c02 = c0*c0;



    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        rhoi = q(i,1);
        mxi  = q(i,2);
        myi  = q(i,3);
        mzi  = q(i,4);
        Pxxi = q(i,5);
        Pxyi = q(i,6);
        Pxzi = q(i,7);
        Pyyi = q(i,8);
        Pyzi = q(i,9);
        Pzzi = q(i,10);

        rhoe = q(i,11);
        mxe  = q(i,12);
        mye  = q(i,13);
        mze  = q(i,14);
        Pxxe = q(i,15);
        Pxye = q(i,16);
        Pxze = q(i,17);
        Pyye = q(i,18);
        Pyze = q(i,19);
        Pzze = q(i,20);

        Ex  = q(i,21);
        Ey  = q(i,22);
        Ez  = q(i,23);
        Bx  = q(i,24);
        By  = q(i,25);
        Bz  = q(i,26);

        ui   = mxi/rhoi;
        vi   = myi/rhoi;
        wi   = mzi/rhoi;
        pxxi = Pxxi - rhoi*ui*ui;
        pxyi = Pxyi - rhoi*ui*vi;
        pxzi = Pxzi - rhoi*ui*wi;
        pyyi = Pyyi - rhoi*vi*vi;
        pyzi = Pyzi - rhoi*vi*wi;
        pzzi = Pzzi - rhoi*wi*wi;

        ue   = mxe/rhoe;
        ve   = mye/rhoe;
        we   = mze/rhoe;
        pxxe = Pxxe - rhoe*ue*ue;
        pxye = Pxye - rhoe*ue*ve;
        pxze = Pxze - rhoe*ue*we;
        pyye = Pyye - rhoe*ve*ve;
        pyze = Pyze - rhoe*ve*we;
        pzze = Pzze - rhoe*we*we;



        // compute flux
        fx(i,1) = rhoi*ui;
        fx(i,2) = rhoi*ui*ui + pxxi;
        fx(i,3) = rhoi*ui*vi + pxyi;
        fx(i,4) = rhoi*ui*wi + pxzi;
        fx(i,5) = rhoi*ui*ui*ui + 3*ui*pxxi;
        fx(i,6) = rhoi*ui*vi*ui + 2*ui*pxyi +   vi*pxxi;
        fx(i,7) = rhoi*ui*wi*ui + 2*ui*pxzi +   wi*pxxi;
        fx(i,8) = rhoi*ui*vi*vi +   ui*pyyi + 2*vi*pxyi;
        fx(i,9) = rhoi*ui*wi*vi +   ui*pyzi +   vi*pxzi + wi*pxyi;
        fx(i,10)= rhoi*ui*wi*wi +   ui*pzzi + 2*wi*pxzi;

        fx(i,11) = rhoe*ue;
        fx(i,12) = rhoe*ue*ue + pxxe;
        fx(i,13) = rhoe*ue*ve + pxye;
        fx(i,14) = rhoe*ue*we + pxze;
        fx(i,15) = rhoe*ue*ue*ue + 3*ue*pxxe;
        fx(i,16) = rhoe*ue*ve*ue + 2*ue*pxye +   ve*pxxe;
        fx(i,17) = rhoe*ue*we*ue + 2*ue*pxze +   we*pxxe;
        fx(i,18) = rhoe*ue*ve*ve +   ue*pyye + 2*ve*pxye;
        fx(i,19) = rhoe*ue*we*ve +   ue*pyze +   ve*pxze + we*pxye;
        fx(i,20) = rhoe*ue*we*we +   ue*pzze + 2*we*pxze;

        fx(i,21) = 0; 
        fx(i,22) = +c02*Bz;
        fx(i,23) = -c02*By;

        fx(i,24) = 0; 
        fx(i,25) = -Ez;
        fx(i,26) = +Ey;

    }
}
