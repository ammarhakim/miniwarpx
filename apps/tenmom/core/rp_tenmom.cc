#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"

extern
void 
rp_10m(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
       FArray<double>& wave, FArray<double>& s, 
       FArray<double>& amdq, FArray<double>& apdq);

extern
void 
rp_maxwell(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
           FArray<double>& wave, FArray<double>& s, 
           FArray<double>& amdq, FArray<double>& apdq);

void
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    int start[4]; // for use in indicating starting indices of sliced arrays
    // arrays for use in slicing
    FArray<double> ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s;

    start[0] = 1-mbc; // first index of all sliced arrays is 1-mbc
    start[1] = 1; // second start index is 1
    start[2] = 1; // third start index is 1

    // slice arrays to send to ten moment ion fluid RP solver
    ql_s = ql.slice( Range(1-mbc,mx+mbc, 1,10), start );
    qr_s = qr.slice( Range(1-mbc,mx+mbc, 1,10), start );
    df_s = df.slice( Range(1-mbc,mx+mbc, 1,10), start );

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 1,10, 1,10), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 1,10), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 1,10), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 1,10), start );

    rd.meqn = 10;
    rd.mwave = 10;

    //
    // call RP solver for ten moment equations for ion fluid
    //
    rp_10m(rd,ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s);

    // slice arrays to send to ten moment electron fluid RP solver

    ql_s = ql.slice( Range(1-mbc,mx+mbc, 11,20), start );
    qr_s = qr.slice( Range(1-mbc,mx+mbc, 11,20), start );
    df_s = df.slice( Range(1-mbc,mx+mbc, 11,20), start );

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 11,20, 11,20), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 11,20), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 11,20), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 11,20), start );

    rd.meqn = 10;
    rd.mwave = 10;
    //
    // call RP solver for ten moment equations for electron fluid
    //
    rp_10m(rd,ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s);   

    // slice arrays to send to Maxwell RP solver
    ql_s = ql.slice( Range(1-mbc,mx+mbc, 21,26), start );
    qr_s = qr.slice( Range(1-mbc,mx+mbc, 21,26), start );
    df_s = df.slice( Range(1-mbc,mx+mbc, 21,26), start );

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 21,26, 21,23), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 21,26), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 21,26), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 21,23), start );

    rd.meqn = 6;
    rd.mwave = 3;
    //
    // call RP solver for Maxwell equations
    //
    rp_maxwell(rd,ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s);
    
    // restore values of meqn and mwave
    rd.meqn = 26;
    rd.mwave = 23;
         
}
