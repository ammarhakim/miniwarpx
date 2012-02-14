#include "miniwarpx.h"

extern 
void src_p(const Run_Data& rd, FArray<double>& sr, FArray<double>& q);

void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;
    FArray<double> sr_p(WxRange(1,26)), q_p(WxRange(1,26));

    // loop over gird points computing source terms
    for (int i=1; i<=mx; i++)
    {
        // copy conserved variables into q_p
        for(int j=1; j<=26; ++j)
            q_p(j) = q(i,j);

        // make call to compute source terms
        src_p(rd, sr_p, q_p);

        // copy sources into sr
        for(int j=1; j<=26; ++j)
            sr(i,j) = sr_p(j);
    }
}
