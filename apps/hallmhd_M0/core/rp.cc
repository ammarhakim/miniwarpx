#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"
#include "hallmhd.h"

void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    double dx = rd.dx;

    FArray<double> delta(Range(1,8));
    FArray<double> alpha(Range(1,8));
    FArray<double> L(Range(1,8,1,8));
    FArray<double> R(Range(1,8,1,8));

    // speed of light
    Hallmhd_Vars *hv = (Hallmhd_Vars*) rd.mvar;
    double c0 = hv->c0;
    double ur = hv->ur;
    double rli = hv->rli;
    double gas_gamma = hv->gas_gamma;
    double ne,ee,gradpex;
    double ex,bx,by,bz,b;

    int start[4]; // for use in indicating starting indices of sliced arrays
    // arrays for use in slicing
    FArray<double> ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s;

    start[0] = 1-mbc; // first index of all sliced arrays is 1-mbc
    start[1] = 1; // second start index is 1
    start[2] = 1; // third start index is 1

    //
    // call RP solver for Euler equations for ion fluid
    //

    // slice arrays to send to Euler RP solver
    ql_s = ql.slice( Range(1-mbc,mx+mbc, 1,5), start );
    qr_s = qr.slice( Range(1-mbc,mx+mbc, 1,5), start );
    df_s = df.slice( Range(1-mbc,mx+mbc, 1,5), start );

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 1,5, 1,3), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 1,5), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 1,5), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 1,3), start );

    rd.meqn = 5;
    rd.mwave = 3;

    rp_euler(rd,ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s);

    // compute waves
    rd.mwave = 8;
    rd.meqn = 8;

    // Compute waves for electron fluid and Maxwell's eqns
    for(int i=2-mbc; i<=mx+mbc-1; i++)
    {
        ne = (qr(i-1,6)+ql(i,6))/2.;        
        ee = (qr(i-1,7)+ql(i,7))/2.;
        gradpex = (gas_gamma-1.)*(ql(i,7)-qr(i-1,7))/dx;

        ex = (qr(i-1,8)+ql(i,8))/2.;
        bx = (qr(i-1,11)+ql(i,11))/2.;
        by = (qr(i-1,12)+ql(i,12))/2.;
        bz = (qr(i-1,13)+ql(i,13))/2.;
        b  = sqrt(bx*bx + by*by + bz*bz);

        // specify eigenvalues
        s(i,4) = 0.;
        s(i,5) = 0.;
        s(i,6) = ex/b;
        s(i,7) = gas_gamma*ex/b + gas_gamma*gradpex*rli/(ne*b); 
        s(i,8) = -c0/ur;
        s(i,9) = -c0/ur;
        s(i,10) = c0/ur;
        s(i,11) = c0/ur;

        // specify left eigenvectors
        L(1,1)=0;

        L(1,2)=0;

        L(1,3)=0;

        L(1,4)=0;

        L(1,5)=0;

        L(1,6)=1;

        L(1,7)=0;

        L(1,8)=0;

        L(2,1)=0;

        L(2,2)=0;

        L(2,3)=1;

        L(2,4)=0;

        L(2,5)=0;

        L(2,6)=0;

        L(2,7)=0;

        L(2,8)=0;

        if (ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) L(3,1)=0;
        else
            L(3,1)=(ee*gas_gamma*gradpex*
                    rli)/(ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli));

        L(3,2)=0;

        if (ex*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) L(3,3)=0;
        else
            L(3,3)=(ee*gas_gamma*gradpex*
                    rli)/(ex*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli));

        if (b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                           pow(ex,2)*pow(ur,2))==0.) L(3,4)=0;
        else
            L(3,4)=(bz*ee*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*
                    pow(ur,2))/(b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                               pow(ex,2)*pow(ur,2)));

        if (b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                           pow(ex,2)*pow(ur,2))==0.) L(3,5)=0;
        else
            L(3,5)=-((by*ee*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*
                      pow(ur,2))/(b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                                 pow(ex,2)*pow(ur,2))));

        if (pow(b,2)*ex*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) L(3,6)=0;
        else
            L(3,6)=-((bx*ee*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli))/(pow(b,2)*ex*
                      ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)));

        if (pow(b,2)*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                  pow(ex,2)*pow(ur,2))==0.) L(3,7)=0;
        else
            L(3,7)=(by*ee*ex*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*pow(ur,2))/(pow(b,2)*
                      ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                   pow(ex,2)*pow(ur,2)));

        if (pow(b,2)*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                  pow(ex,2)*pow(ur,2))==0.) L(3,8)=0;
        else
            L(3,8)=(bz*ee*ex*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*pow(ur,2))/(pow(b,2)*
                      ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                   pow(ex,2)*pow(ur,2)));

        if (ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) L(4,1)=0;
        else
            L(4,1)=-((ee*gas_gamma*gradpex*rli)/(ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)));

        L(4,2)=1;

        if (ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli==0.) L(4,3)=0;
        else
            L(4,3)=(ee*(-1 + gas_gamma)*ne)/(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli);

        if (b*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                    pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))==0.) L(4,4)=0;
        else
            L(4,4)=(bz*ee*(-1 + gas_gamma)*gas_gamma*ne*pow(ex*ne + gradpex*rli,2)*pow(ur,2))/(b*(ex*(-1 + gas_gamma)*ne + 
                       gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                           pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2)));

        if (b*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                                        pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))==0.) L(4,5)=0;
        else
            L(4,5)=-((by*ee*(-1 + gas_gamma)*gas_gamma*ne*pow(ex*ne + gradpex*rli,2)*pow(ur,2))/(b*(ex*(-1 + gas_gamma)*ne + 
                       gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                           pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))));

        if (pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) L(4,6)=0;
        else
            L(4,6)=-((bx*ee*(-1 + gas_gamma)*(ex*ne + gradpex*rli))/(pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)));

        if (pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                      pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))==0.) L(4,7)=0;
        else
            L(4,7)=(by*ee*(-1 + gas_gamma)*pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,3)*
                pow(ur,2))/(pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                      pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2)));

        if (pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                      pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))==0.) L(4,8)=0;
        else
            L(4,8)=(bz*ee*(-1 + gas_gamma)*pow(gas_gamma,2)* pow(ex*ne + gradpex*rli,3)*
                pow(ur,2))/(pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                      pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2)));

        L(5,1)=0;

        L(5,2)=0;

        L(5,3)=0;

        L(5,4)=-ur/(2.*c0);

        L(5,5)=0;

        L(5,6)=0;

        L(5,7)=0;

        L(5,8)=0.5;

        L(6,1)=0;

        L(6,2)=0;

        L(6,3)=0;

        L(6,4)=0;

        L(6,5)=ur/(2.*c0);

        L(6,6)=0;


        L(6,7)=0.5;

        L(6,8)=0;

        L(7,1)=0;

        L(7,2)=0;

        L(7,3)=0;

        L(7,4)=ur/(2.*c0);

        L(7,5)=0;

        L(7,6)=0;

        L(7,7)=0;

        L(7,8)=0.5;

        L(8,1)=0;

        L(8,2)=0;

        L(8,3)=0;

        L(8,4)=0;

        L(8,5)=-ur/(2.*c0);

        L(8,6)=0;

        L(8,7)=0.5;

        L(8,8)=0;

        // specify right eigenvectors
        if (pow(b,2)*ex==0.) R(1,1)=0;
        else R(1,1)=(bx*(ex*ne + gradpex*rli))/(pow(b,2)*ex);

        if (ex==0.) R(1,2)=0;
        else R(1,2)=-(ne/ex);

        if (ee*gas_gamma*gradpex*rli==0.) R(1,3)=0;
        else R(1,3)=(ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli))/(ee*gas_gamma*gradpex*rli);

        R(1,4)=0;

        if (pow(b,2)*(b*c0 + ex*ur)==0.) R(1,5)=0;
        else R(1,5)=(bz*(ex*ne + gradpex*rli)*ur)/(pow(b,2)*(b*c0 + ex*ur));

        if (pow(b,2)*(b*c0 + ex*ur)==0.) R(1,6)=0;
        else R(1,6)=(by*(ex*ne + gradpex*rli)*ur)/(pow(b,2)*(b*c0 + ex*ur));

        if (pow(b,2)*(b*c0 - ex*ur)==0.) R(1,7)=0;
        else R(1,7)=-((bz*(ex*ne + gradpex*rli)*ur)/(pow(b,2)*(b*c0 - ex*ur)));

        if (pow(b,2)*(b*c0 - ex*ur)==0.) R(1,8)=0;
        else R(1,8)=-((by*(ex*ne + gradpex*rli)*ur)/(pow(b,2)*(b*c0 - ex*ur)));

        if (pow(b,2)*ex*ne==0.) R(2,1)=0;
        else R(2,1)=(bx*ee*(ex*ne + gradpex*rli))/(pow(b,2)*ex*ne);

        if (ex==0.) R(2,2)=0;
        else R(2,2)=-(ee/ex);

        R(2,3)=1;

        R(2,4)=1;

        if (pow(b,2)*ne*(b*c0 + ex*ur)*(b*c0*ne + gas_gamma*(ex*ne + gradpex*rli)*ur)==0.) R(2,5)=0;
        else
            R(2,5)=(bz*ee*gas_gamma*(ex*ne + gradpex*rli)*
                    ur*(b*c0*ne + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*ne*(b*c0 + ex*ur)*(b*c0*ne + 
                                                          gas_gamma*(ex*ne + gradpex*rli)*ur));

        if (pow(b,2)*ne*(b*c0 + ex*ur)*(b*c0*ne + gas_gamma*(ex*ne + gradpex*rli)*ur)==0.) R(2,6)=0;
        else
            R(2,6)=(by*ee*gas_gamma*(ex*ne + gradpex*rli)*
                    ur*(b*c0*ne + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*ne*(b*c0 + ex*ur)*(b*c0*ne + 
                                                          gas_gamma*(ex*ne + gradpex*rli)*ur));

        if (pow(b,2)*ne*(b*c0 - ex*ur)*(b*c0*ne - gas_gamma*(ex*ne + gradpex*rli)*ur)==0.) R(2,7)=0;
        else
            R(2,7)=(bz*ee*gas_gamma*(ex*ne + gradpex*rli)*ur*(-(b*c0*ne) + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*
                                     ne*(b*c0 - ex*ur)*(b*c0*ne - gas_gamma*(ex*ne + gradpex*rli)*ur));

        if (pow(b,2)*ne*(b*c0 - ex*ur)*(b*c0*ne - gas_gamma*(ex*ne + gradpex*rli)*ur)==0.) R(2,8)=0;
        else
            R(2,8)=(by*ee*gas_gamma*(ex*ne + gradpex*rli)*ur*(-(b*c0*ne) + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*
                                     ne*(b*c0 - ex*ur)*(b*c0*ne - gas_gamma*(ex*ne + gradpex*rli)*ur));

        R(3,1)=0;

        R(3,2)=1;

        R(3,3)=0;

        R(3,4)=0;

        R(3,5)=0;

        R(3,6)=0;

        R(3,7)=0;

        R(3,8)=0;

        R(4,1)=0;

        R(4,2)=0;

        R(4,3)=0;

        R(4,4)=0;

        R(4,5)=-(c0/ur);

        R(4,6)=0;

        R(4,7)=c0/ur;

        R(4,8)=0;

        R(5,1)=0;

        R(5,2)=0;

        R(5,3)=0;

        R(5,4)=0;

        R(5,5)=0;

        R(5,6)=c0/ur;

        R(5,7)=0;

        R(5,8)=-(c0/ur);

        R(6,1)=1;

        R(6,2)=0;

        R(6,3)=0;

        R(6,4)=0;

        R(6,5)=0;

        R(6,6)=0;

        R(6,7)=0;

        R(6,8)=0;

        R(7,1)=0;

        R(7,2)=0;

        R(7,3)=0;

        R(7,4)=0;

        R(7,5)=0;

        R(7,6)=1;

        R(7,7)=0;

        R(7,8)=1;

        R(8,1)=0;

        R(8,2)=0;

        R(8,3)=0;

        R(8,4)=0;

        R(8,5)=1;

        R(8,6)=0;

        R(8,7)=1;

        R(8,8)=0;


        for (int n1=1; n1<=rd.meqn; n1++)
           delta(n1) = df(i,n1+5);

        // compute coefficients of the 13 eigenvectors
        for (int n1=1; n1<=rd.mwave; n1++)
        {
            alpha(n1) = 0.;
            for (int n2=1; n2<=rd.meqn; n2++)
                alpha(n1) = alpha(n1) + L(n1,n2)*delta(n2);
        }

        // compute waves

        for (int n1=1; n1<=rd.mwave; n1++)
           for (int n2=1; n2<=rd.meqn; n2++)
           {
               wave(i,n2+5,n1+3) = alpha(n1)*R(n2,n1);
           }

    }

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 6,13, 4,11), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 6,13), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 6,13), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 4,11), start );
    // compute fluctuations
    eval_fluctuations(rd,wave_s,s_s,amdq_s,apdq_s);
    rd.meqn  =13;
    rd.mwave =11;
}
