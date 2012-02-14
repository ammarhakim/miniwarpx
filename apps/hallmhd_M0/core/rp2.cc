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
    int mwave = rd.mwave;
    int mbc = rd.mbc;
    double dx = rd.dx;

    FArray<double> delta(Range(1,13));
    FArray<double> alpha(Range(1,13));
    FArray<double> L(Range(1,13,1,13));
    FArray<double> R(Range(1,13,1,13));

    // speed of light
    Hallmhd_Vars *hv = (Hallmhd_Vars*) rd.mvar;
    double c0 = hv->c0;
    double ur = hv->ur;
    double rli = hv->rli;
    double gas_gamma = hv->gas_gamma;
    double ne,ee,gradpex,alp,wa;
    double ni,ui,vi,wi,ei,ex,bx,by,bz,b,rhsqrtl,rhsqrtr,rhsq2,a;


    // compute waves
    for(int i=2-mbc; i<=mx+mbc-1; i++)
    {
        rhsqrtl = sqrt(qr(i-1,1));
        rhsqrtr = sqrt(ql(i,1));
        rhsq2 = rhsqrtl + rhsqrtr;
        ni = (qr(i-1,1)+ql(i,1))/2.;
        ui = (qr(i-1,2)/rhsqrtl + ql(i,2)/rhsqrtr) / rhsq2;
        vi = (qr(i-1,3)/rhsqrtl + ql(i,3)/rhsqrtr) / rhsq2;
        wi = (qr(i-1,4)/rhsqrtl + ql(i,4)/rhsqrtr) / rhsq2;
        ei = (qr(i-1,5)+ql(i,5))/2.;

        ne = (qr(i-1,6)+ql(i,6))/2.;        
        ee = (qr(i-1,7)+ql(i,7))/2.;
        gradpex = (gas_gamma-1.)*(ql(i,7)-qr(i-1,7))/dx;

        ex = (qr(i-1,8)+ql(i,8))/2.;
        bx = (qr(i-1,11)+ql(i,11))/2.;
        by = (qr(i-1,12)+ql(i,12))/2.;
        bz = (qr(i-1,13)+ql(i,13))/2.;
        b  = sqrt(bx*bx + by*by + bz*bz);

        a  = sqrt(gas_gamma*(gas_gamma-1)*(2*ei-ni*ui*ui)/(2*ni));

        // specify eigenvalues
        s(i,1) = 0.;
        s(i,2) = 0.;
        s(i,3) = ex/b;
        s(i,4) = gas_gamma*ex/b + gas_gamma*gradpex*rli/(ne*b); 
        s(i,5) = ui; 
        s(i,6) = ui; 
        s(i,7) = ui;
        s(i,8) = ui-a;
        s(i,9) = ui+a;
        s(i,10) = -c0/ur;
        s(i,11) = -c0/ur;
        s(i,12) = c0/ur;
        s(i,13) = c0/ur;

        // specify left eigenvectors
        L(1,1)=0;
        L(1,2)=0;
        L(1,3)=0;
        L(1,4)=0;
        L(1,5)=0;
        L(1,6)=0;
        L(1,7)=0;
        L(1,8)=0;
        L(1,9)=0;
        L(1,10)=0;
        L(1,11)=1;
        L(1,12)=0;
        L(1,13)=0;

        L(2,1)=0;
        L(2,2)=0;
        L(2,3)=0;
        L(2,4)=0;
        L(2,5)=0;
        L(2,6)=0;
        L(2,7)=0;
        L(2,8)=1;
        L(2,9)=0;
        L(2,10)=0;
        L(2,11)=0;
        L(2,12)=0;
        L(2,13)=0;

        if  (ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.)  L(3,1)=0;
        else L(3,1)=(ee*gas_gamma*gradpex*
                     rli)/(ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli));

        L(3,2)=0;
        L(3,3)=0;
        L(3,4)=0;
        L(3,5)=0;
        L(3,6)=0;
        L(3,7)=0;

        if (ex*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) L(3,8)=0;
        else L(3,8)=(ee*gas_gamma*gradpex*
                rli)/(ex*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli));

        if (b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - pow(ex,2)*pow(ur,2))==0.) 
            L(3,9)=0;
        else L(3,9)=(bz*ee*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*
                     pow(ur,2))/(b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                                        pow(ex,2)*pow(ur,2)));

        if (b*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - pow(ex,2)*pow(ur,2))==0.)
            L(3,10)=0;
        else L(3,10)=-((by*ee*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*
                        pow(ur,2))/(b*ne*(ex*(-1 + gas_gamma)*ne + 
                                          gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                  pow(ex,2)*pow(ur,2))));

        if (pow(b,2)*ex*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.)
            L(3,11)=0;
        else  L(3,11)=-((bx*ee*gas_gamma*gradpex*
                         rli*(ex*ne + gradpex*rli))/(pow(b,2)*ex*
                                                     ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)));
        
        if(pow(b,2)*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                         pow(ex,2)*pow(ur,2))==0.)
            L(3,12)=0;
        else L(3,12)=(by*ee*ex*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*
                      pow(ur,2))/(pow(b,2)*
                                  ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                                       pow(ex,2)*pow(ur,2)));
           
        if (pow(b,2)*ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                          pow(ex,2)*pow(ur,2))==0.)
            L(3,13)=0;
        else L(3,13)=(bz*ee*ex*gas_gamma*gradpex*rli*(ex*ne + gradpex*rli)*
                      pow(ur,2))/(pow(b,2)*
                                  ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2) - 
                                                                                       pow(ex,2)*pow(ur,2)));
           
        if (ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.) 
            L(4,1)=0;
        else L(4,1)=-((ee*gas_gamma*gradpex*
                       rli)/(ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)));
        L(4,2)=1;
        L(4,3)=0;
        L(4,4)=0;
        L(4,5)=0;
        L(4,6)=0;
        L(4,7)=0;

        if (ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli==0.) 
            L(4,8)=0;
        else L(4,8)=(ee*(-1 + gas_gamma)*ne)/(ex*(-1 + gas_gamma)*ne + 
                                              gas_gamma*gradpex*rli);

        if (b*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                     pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))==0.)
            L(4,9)=0;
        else L(4,9)=(bz*ee*(-1 + gas_gamma)*gas_gamma*ne*pow(ex*ne + gradpex*rli,2)*
                     pow(ur,2))/(b*(ex*(-1 + gas_gamma)*ne + 
                                    gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                                            pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*
                                                            pow(ur,2)));

        if (b*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                     pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*pow(ur,2))==0.)
            L(4,10)=0;
        else L(4,10)=-((by*ee*(-1 + gas_gamma)*gas_gamma*ne*pow(ex*ne + gradpex*rli,2)*
                     pow(ur,2))/(b*(ex*(-1 + gas_gamma)*ne + 
                                    gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                                            pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*
                                                            pow(ur,2))));

        if (pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)==0.)
            L(4,11)=0;
        else L(4,11)=-((bx*ee*(-1 + gas_gamma)*(ex*ne + gradpex*rli))/(pow(b,
                        2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)));

        if (pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                                            pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*
                                                            pow(ur,2))==0.)
            L(4,12)=0;
        else L(4,12)=(by*ee*(-1 + gas_gamma)*pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,3)*
                      pow(ur,2))/(pow(b,2)*(ex*(-1 + gas_gamma)*ne + 
                      gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                              pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*
                                              pow(ur,2)));

        if (pow(b,2)*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                                            pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*
                                                            pow(ur,2))==0.)
            L(4,13)=0;
        else L(4,13)=(bz*ee*(-1 + gas_gamma)*pow(gas_gamma,2)*
                      pow(ex*ne + gradpex*rli,3)*pow(ur,2))/(pow(b,2)*(ex*(-1 + gas_gamma)*ne + 
                      gas_gamma*gradpex*rli)*(pow(b,2)*pow(c0,2)*pow(ne,2) - 
                                              pow(gas_gamma,2)*pow(ex*ne + gradpex*rli,2)*
                                              pow(ur,2)));

        L(5,1)=0;
        L(5,2)=0;
        L(5,3)=(2*ei*gas_gamma*pow(ui,2) - (1 + gas_gamma)*ni*pow(ui,4))/(4*ei*
                                                                      gas_gamma - 2*gas_gamma*ni*pow(ui,2));
        L(5,4)=(ni*pow(ui,3))/(2*ei*gas_gamma - gas_gamma*ni*pow(ui,2));
        L(5,5)=0;
        L(5,6)=0;
        L(5,7)=(ni*pow(ui,2))/(gas_gamma*(-2*ei + ni*pow(ui,2)));
        L(5,8)=0;
        L(5,9)=0;
        L(5,10)=0;
        L(5,11)=0;
        L(5,12)=0;
        L(5,13)=0;
        
        L(6,1)=0;
        L(6,2)=0;
        L(6,3)=(ni*pow(ui,2)*
                wi)/(gas_gamma*(-2*ei + ni*pow(ui,2)));
        L(6,4)=(2*ni*ui*wi)/(2*ei*gas_gamma - 
                             gas_gamma*ni*pow(ui,2));
        L(6,5)=0;
        L(6,6)=1;
        L(6,7)=(-2*ni*wi)/(2*ei*gas_gamma - 
                           gas_gamma*ni*pow(ui,2));
        L(6,8)=0;
        L(6,9)=0;
        L(6,10)=0;
        L(6,11)=0;
        L(6,12)=0;
        L(6,13)=0;
        
        L(7,1)=0;
        L(7,2)=0;
        L(7,3)=(ni*pow(ui,2)*
                vi)/(gas_gamma*(-2*ei + ni*pow(ui,2)));
        L(7,4)=(2*ni*ui*vi)/(2*ei*gas_gamma - 
                             gas_gamma*ni*pow(ui,2));
        L(7,5)=1;
        L(7,6)=0;
        L(7,7)=(-2*ni*vi)/(2*ei*gas_gamma - 
                           gas_gamma*ni*pow(ui,2));
        L(7,8)=0;
        L(7,9)=0;
        L(7,10)=0;
        L(7,11)=0;
        L(7,12)=0;
        L(7,13)=0;
        
        L(8,1)=0;
        L(8,2)=0;
        L(8,3)=(ui*(-((-1 + gas_gamma)*ni*
                      pow(ui,
                            2)*((1 + gas_gamma)*ni*ui - 
                                2*sqrt(2)*
                                sqrt(-((-1 + gas_gamma)*gas_gamma*
                                       ni*(-2*ei + ni*pow(ui,2)))))) - 
                    2*ei*gas_gamma*(ni*(ui - gas_gamma*ui) + 
                                sqrt(2)*
                                sqrt(-((-1 + gas_gamma)*gas_gamma*
                                       ni*(-2*ei + ni*pow(ui,2)))))))/(4.*(-1 + gas_gamma)*
                                                                         gas_gamma*ni*(-2*ei + ni*pow(ui,2)));
        L(8,4)=(-4*pow(ei,2)*pow(gas_gamma,2) + 
                2*ei*gas_gamma*(-3 + 4*gas_gamma)*ni*pow(ui,2) - 
                ni*pow(ui,
                         3)*(3*(-1 + gas_gamma)*gas_gamma*ni*ui + 
                             sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                            ni*(-2*ei + ni*pow(ui,2))))))/(2.*sqrt(2)*
                                                                             gas_gamma*(2*ei - ni*pow(ui,2))*
                                                                             sqrt(-((-1 + gas_gamma)*gas_gamma*
                                                                                    ni*(-2*ei + ni*pow(ui,2)))));
        L(8,5)=0;
        L(8,6)=0;
        L(8,7)=(2*ei*gas_gamma - 
                ui*((-1 + gas_gamma)*ni*ui + 
                    sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                   ni*(-2*ei + ni*pow(ui,2))))))/(4*ei*gas_gamma - 
                                                                    2*gas_gamma*ni*pow(ui,2));
        L(8,8)=0;
        L(8,9)=0;
        L(8,10)=0;
        L(8,11)=0;
        L(8,12)=0;
        L(8,13)=0;
        
        L(9,1)=0;
        L(9,2)=0;
        L(9,3)=(ui*(-4*ei*gas_gamma + ui*(2*gas_gamma*ni*ui + sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                     ni*(-2*ei + ni*pow(ui,2))))))*(2*ei*gas_gamma + ui*(ni*(ui - gas_gamma*ui) 
                     +sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*ni*(-2*ei + ni*pow(ui,2)))))))/(4.*sqrt(2)*
                     gas_gamma*(2*ei - ni*pow(ui,2))*sqrt(-((-1 + gas_gamma)*gas_gamma*ni*(-2*ei + ni*pow(ui,2)))));
        L(9,4)=(4*pow(ei,2)*pow(gas_gamma,2) - 2*ei*gas_gamma*(-3 + 4*gas_gamma)*ni*pow(ui,2) + 
                ni*pow(ui,3)*(3*(-1 + gas_gamma)*gas_gamma*ni*ui - sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                ni*(-2*ei + ni*pow(ui,2))))))/(2.*sqrt(2)*gas_gamma*(2*ei - ni*pow(ui,2))*
                sqrt(-((-1 + gas_gamma)*gas_gamma*ni*(-2*ei + ni*pow(ui,2)))));
        L(9,5)=0;
        L(9,6)=0;
        L(9,7)=(2*ei*gas_gamma + 
                ui*(ni*(ui - gas_gamma*ui) + 
                    sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                   ni*(-2*ei + ni*pow(ui,2))))))/(4*ei*gas_gamma - 
                                                                    2*gas_gamma*ni*pow(ui,2));
        L(9,8)=0;
        L(9,9)=0;
        L(9,10)=0;
        L(9,11)=0;
        L(9,12)=0;
        L(9,13)=0;

        L(10,1)=0;
        L(10,2)=0;
        L(10,3)=0;
        L(10,4)=0;
        L(10,5)=0;
        L(10,6)=0;
        L(10,7)=0;
        L(10,8)=0;
        L(10,9)=-ur/(2.*c0);
        L(10,10)=0;
        L(10,11)=0;
        L(10,12)=0;
        L(10,13)=0.5;

        L(11,1)=0;
        L(11,2)=0;
        L(11,3)=0;
        L(11,4)=0;
        L(11,5)=0;
        L(11,6)=0;
        L(11,7)=0;
        L(11,8)=0;
        L(11,9)=0;
        L(11,10)=ur/(2.*c0);
        L(11,11)=0;
        L(11,12)=0.5;
        L(11,13)=0;

        L(12,1)=0;
        L(12,2)=0;
        L(12,3)=0;
        L(12,4)=0;
        L(12,5)=0;
        L(12,6)=0;
        L(12,7)=0;
        L(12,8)=0;
        L(12,9)=ur/(2.*c0);
        L(12,10)=0;
        L(12,11)=0;
        L(12,12)=0;
        L(12,13)=0.5;

        L(13,1)=0;
        L(13,2)=0;
        L(13,3)=0;
        L(13,4)=0;
        L(13,5)=0;
        L(13,6)=0;
        L(13,7)=0;
        L(13,8)=0;
        L(13,9)=0;
        L(13,10)=-ur/(2.*c0);
        L(13,11)=0;
        L(13,12)=0.5;
        L(13,13)=0;

        // specify right eigenvectors
        if (pow(b,2)*ex==0.) R(1,1)=0;
        else R(1,1)=(bx*(ex*ne + gradpex*rli))/(pow(b,2)*ex);

        if (ex==0.) R(1,2)=0;
        else R(1,2)=-(ne/ex);

        if (ee*gas_gamma*gradpex*rli==0.) R(1,3)=0;
        else R(1,3)=(ne*(ex*(-1 + gas_gamma)*ne + gas_gamma*gradpex*rli))/(ee*gas_gamma*
                                                              gradpex*rli);
        R(1,4)=0;
        R(1,5)=0;
        R(1,6)=0;
        R(1,7)=0;
        R(1,8)=0;
        R(1,9)=0;
        if (pow(b,2)*(b*c0 + ex*ur)==0.) R(1,10)=0;
        else R(1,10)=(bz*(ex*ne + gradpex*rli)*
                      ur)/(pow(b,2)*(b*c0 + ex*ur));
        if (pow(b,2)*(b*c0 + ex*ur)==0.) R(1,11)=0.;
        else R(1,11)=(by*(ex*ne + gradpex*rli)*
                      ur)/(pow(b,2)*(b*c0 + ex*ur));
        if (pow(b,2)*(b*c0 - ex*ur)==0.) R(1,12)=0;
        else R(1,12)=-((bz*(ex*ne + gradpex*rli)*
                        ur)/(pow(b,2)*(b*c0 - ex*ur)));
        if (pow(b,2)*(b*c0 - ex*ur)==0.) R(1,13)=0;
        else R(1,13)=-((by*(ex*ne + gradpex*rli)*
                        ur)/(pow(b,2)*(b*c0 - ex*ur)));

        if (pow(b,2)*ex*ne==0.) R(2,1)=0;
        else R(2,1)=(bx*ee*(ex*ne + gradpex*rli))/(pow(b,2)*ex*ne);
        if (ex==0.) R(2,2)=0;
        else R(2,2)=-(ee/ex);
        R(2,3)=1;
        R(2,4)=1;
        R(2,5)=0;
        R(2,6)=0;
        R(2,7)=0;
        R(2,8)=0;
        R(2,9)=0;
        if (pow(b,2)*ne*(b*c0 + ex*ur)*(b*c0*ne + gas_gamma*(ex*ne + gradpex*rli)*ur)==0.)
            R(2,10)=0;
        else
            R(2,10)=(bz*ee*gas_gamma*(ex*ne + gradpex*rli)*
                 ur*(b*c0*ne + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*
                                                           ne*(b*c0 + ex*ur)*(b*c0*ne + 
                                                                              gas_gamma*(ex*ne + gradpex*rli)*ur));
        if (pow(b,2)*ne*(b*c0 + ex*ur)*(b*c0*ne + gas_gamma*(ex*ne + gradpex*rli)*ur)==0.)
            R(2,11)=0;
        else
            R(2,11)=(by*ee*gas_gamma*(ex*ne + gradpex*rli)*
                 ur*(b*c0*ne + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*
                                                           ne*(b*c0 + ex*ur)*(b*c0*ne + 
                                                                             gas_gamma*(ex*ne + gradpex*rli)*ur));
        if (pow(b,2)*ne*(b*c0 - ex*ur)*(b*c0*ne - gas_gamma*(ex*ne + gradpex*rli)*ur)==0.)
            R(2,12)=0;
        else
            R(2,12)=(bz*ee*gas_gamma*(ex*ne + gradpex*rli)*
                 ur*(-(b*c0*ne) + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*
                                                              ne*(b*c0 - ex*ur)*(b*c0*ne - 
                                                                                gas_gamma*(ex*ne + gradpex*rli)*ur));
        if (pow(b,2)*ne*(b*c0 - ex*ur)*(b*c0*ne - gas_gamma*(ex*ne + gradpex*rli)*ur)==0.)
            R(2,13)=0;
        else
            R(2,13)=(by*ee*gas_gamma*(ex*ne + gradpex*rli)*
                 ur*(-(b*c0*ne) + ex*ne*ur + gradpex*rli*ur))/(pow(b,2)*
                                                              ne*(b*c0 - ex*ur)*(b*c0*ne - 
                                                                                gas_gamma*(ex*ne + gradpex*rli)*ur));

        R(3,1)=0;
        R(3,2)=0;
        R(3,3)=0;
        R(3,4)=0;
        R(3,5)=2/pow(ui,2);
        R(3,6)=0;
        R(3,7)=0;
        R(3,8)=(-2*ni)/(-2*ei*gas_gamma + 
                        ui*((-1 + gas_gamma)*ni*ui + 
                            sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                           ni*(-2*ei + ni*pow(ui,2))))));
        R(3,9)=(-2*ni)/(-2*ei*gas_gamma + 
                        ui*((-1 + gas_gamma)*ni*ui - 
                            sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                           ni*(-2*ei + ni*pow(ui,2))))));
        R(3,10)=0;
        R(3,11)=0;
        R(3,12)=0;
        R(3,13)=0;

        R(4,1)=0;
        R(4,2)=0;
        R(4,3)=0;
        R(4,4)=0;
        R(4,5)=2/ui;
        R(4,6)=0;
        R(4,7)=0;
        R(4,8)=(-2*ni*ui + 
                sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*ni*(-2*ei + ni*pow(ui,2)))))/(-2*
                                                                                  ei*gas_gamma + 
                                                                                  ui*((-1 + gas_gamma)*ni*ui + 
                                                                                      sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                                                                                     ni*(-2*ei + ni*pow(ui,2))))));
        R(4,9)=(2*ni*ui + 
                sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*ni*(-2*ei + ni*pow(ui,2)))))/(2*
                                                                                  ei*gas_gamma + 
                                                                                  ui*(ni*(ui - gas_gamma*ui) + 
                                                                                      sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                                                                                     ni*(-2*ei + ni*pow(ui,2))))));
        R(4,10)=0;
        R(4,11)=0;
        R(4,12)=0;
        R(4,13)=0;

        R(5,1)=0;
        R(5,2)=0;
        R(5,3)=0;
        R(5,4)=0;
        R(5,5)=0;
        R(5,6)=0;
        R(5,7)=1;
        R(5,8)=(-2*ni*vi)/(-2*ei*gas_gamma + 
                           ui*((-1 + gas_gamma)*ni*ui + 
                               sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                              ni*(-2*ei + ni*pow(ui,2))))));
        R(5,9)=(-2*ni*vi)/(-2*ei*gas_gamma + 
                           ui*((-1 + gas_gamma)*ni*ui - 
                               sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                              ni*(-2*ei + ni*pow(ui,2))))));
        R(5,10)=0;
        R(5,11)=0;
        R(5,12)=0;
        R(5,13)=0;

        R(6,1)=0;
        R(6,2)=0;
        R(6,3)=0;
        R(6,4)=0;
        R(6,5)=0;
        R(6,6)=1;
        R(6,7)=0;
        R(6,8)=(-2*ni*wi)/(-2*ei*gas_gamma + 
                           ui*((-1 + gas_gamma)*ni*ui + 
                               sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                              ni*(-2*ei + ni*pow(ui,2))))));
        R(6,9)=(-2*ni*wi)/(-2*ei*gas_gamma + 
                           ui*((-1 + gas_gamma)*ni*ui - 
                               sqrt(2)*sqrt(-((-1 + gas_gamma)*gas_gamma*
                                              ni*(-2*ei + ni*pow(ui,2))))));
        R(6,10)=0;
        R(6,11)=0;
        R(6,12)=0;
        R(6,13)=0;
        
        R(7,1)=0;
        R(7,2)=0;
        R(7,3)=0;
        R(7,4)=0;
        R(7,5)=1;
        R(7,6)=0;
        R(7,7)=0;
        R(7,8)=1;
        R(7,9)=1;
        R(7,10)=0;
        R(7,11)=0;
        R(7,12)=0;
        R(7,13)=0;
        
        R(8,1)=0;
        R(8,2)=1;
        R(8,3)=0;
        R(8,4)=0;
        R(8,5)=0;
        R(8,6)=0;
        R(8,7)=0;
        R(8,8)=0;
        R(8,9)=0;
        R(8,10)=0;
        R(8,11)=0;
        R(8,12)=0;
        R(8,13)=0;
        
        R(9,1)=0;
        R(9,2)=0;
        R(9,3)=0;
        R(9,4)=0;
        R(9,5)=0;
        R(9,6)=0;
        R(9,7)=0;
        R(9,8)=0;
        R(9,9)=0;
        R(9,10)=-(c0/ur);
        R(9,11)=0;
        R(9,12)=c0/ur;
        R(9,13)=0;

        R(10,1)=0;
        R(10,2)=0;
        R(10,3)=0;
        R(10,4)=0;
        R(10,5)=0;
        R(10,6)=0;
        R(10,7)=0;
        R(10,8)=0;
        R(10,9)=0;
        R(10,10)=0;
        R(10,11)=c0/ur;
        R(10,12)=0;
        R(10,13)=-(c0/ur);

        R(11,1)=1;
        R(11,2)=0;
        R(11,3)=0;
        R(11,4)=0;
        R(11,5)=0;
        R(11,6)=0;
        R(11,7)=0;
        R(11,8)=0;
        R(11,9)=0;
        R(11,10)=0;
        R(11,11)=0;
        R(11,12)=0;
        R(11,13)=0;
        
        R(12,1)=0;
        R(12,2)=0;
        R(12,3)=0;
        R(12,4)=0;
        R(12,5)=0;
        R(12,6)=0;
        R(12,7)=0;
        R(12,8)=0;
        R(12,9)=0;
        R(12,10)=0;
        R(12,11)=1;
        R(12,12)=0;
        R(12,13)=1;
        
        R(13,1)=0;
        R(13,2)=0;
        R(13,3)=0;
        R(13,4)=0;
        R(13,5)=0;
        R(13,6)=0;
        R(13,7)=0;
        R(13,8)=0;
        R(13,9)=0;
        R(13,10)=1;
        R(13,11)=0;
        R(13,12)=1;
        R(13,13)=0;


        for (int n1=1; n1<=mwave; n1++)
           delta(n1) = df(i,n1);

        // compute coefficients of the 13 eigenvectors
        for (int n1=1; n1<=mwave; n1++)
        {
            alpha(n1) = 0.;
            for (int n2=1; n2<=mwave; n2++)
                alpha(n1) = alpha(n1) + L(n1,n2)*delta(n2);
            alp = alpha(n1);
        }

        // compute waves

        for (int n1=1; n1<=mwave; n1++)
           for (int n2=1; n2<=mwave; n2++)
           {
               wave(i,n2,n1) = alpha(n1)*R(n2,n1);
               wa = wave(i,n2,n1);
           }

    }
    
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);
}
