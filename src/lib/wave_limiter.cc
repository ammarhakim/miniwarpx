#include "miniwarpx.h"
#include "rkdg_algo.h"
#include <math.h>

void
wave_limiter(Run_Data& rd, FArray<double>& wave, FArray<double>& s)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mwave = rd.mwave;

    double c,r,dotr,dotl,wnorm2,wlimitr;

    wlimitr = 0.0;

    // loop over each wave, applying limiters
    for(int mw=1; mw<=mwave; mw++)
    {
        if(rd.limiters[mw-1]!=0)
        {// limiters are to be applied to wave 'mw'

            dotr = 0.0;
            for(int i=0; i<=mx+1; i++)
            {
                wnorm2 = 0.0;
                dotl = dotr;
                dotr = 0.0;
                for(int m=1; m<=meqn; m++)
                {
                    wnorm2 += pow(wave(i,m,mw),2.);
                    dotr += wave(i,m,mw)*wave(i+1,m,mw);
                }
                
                if (i==0) goto loopend;
                if (wnorm2 == 0.0) goto loopend;

                if(s(i,mw) > 0)
                    r = dotl/wnorm2;
                else
                    r = dotr/wnorm2;

                switch(rd.limiters[mw-1])
                {
                    case 0:
                        // this can never happen
                        break;

                    case 1:
                        // minmod limiter
                        wlimitr = dmax(0.0, dmin(1.0,r));
                        break;

                    case 2:
                        // superbee
                        wlimitr = dmax(0.0, dmin(1.,2.*r), dmin(2,r));
                        break;

                    case 3:
                        // van Leer
                        wlimitr = (r+fabs(r))/(1.+fabs(r));
                        break;

                    case 4:
                        // monotized centered
                        c = (1.+r)/2.;
                        wlimitr = dmax(0.0, dmin(c,2.,2.*r));
                        break;

                    case 5:
                        // Beam-Warming
                        wlimitr = r;

                    default:
                        ;
                }

                // apply limiter to waves
                for(int m=1; m<=meqn; m++)
                    wave(i,m,mw) = wlimitr*wave(i,m,mw);

              loopend:
                ;
            }
        }
    }
}
