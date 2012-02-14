#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"
#include "cold.h"

void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;
    int mwave = rd.mwave;

    int start[4]; // for use in indicating starting indices of sliced arrays
    // arrays for use in slicing
    FArray<double> ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s;

    start[0] = 1-mbc; // first index of all sliced arrays is 1-mbc
    start[1] = 1; // second start index is 1
    start[2] = 1; // third start index is 1

    // slice arrays to send to cold-electron RP solver
    ql_s = ql.slice( Range(1-mbc,mx+mbc, 1,4), start );
    qr_s = qr.slice( Range(1-mbc,mx+mbc, 1,4), start );
    df_s = df.slice( Range(1-mbc,mx+mbc, 1,4), start );

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 1,4, 1,2), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 1,4), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 1,4), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 1,2), start );

    rd.meqn = 4;
    rd.mwave = 2;

    //
    // call RP solver for cold electron fluid
    //
    rp_electrons(rd,ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s);

    // slice arrays to send to Maxwell RP solver
    ql_s = ql.slice( Range(1-mbc,mx+mbc, 5,10), start );
    qr_s = qr.slice( Range(1-mbc,mx+mbc, 5,10), start );
    df_s = df.slice( Range(1-mbc,mx+mbc, 5,10), start );

    wave_s = wave.slice( Range(1-mbc,mx+mbc, 5,10, 3,5), start );
    amdq_s = amdq.slice( Range(1-mbc,mx+mbc, 5,10), start );
    apdq_s = apdq.slice( Range(1-mbc,mx+mbc, 5,10), start );

    s_s = s.slice( Range(1-mbc,mx+mbc, 3,5), start );

    rd.meqn = 6;
    rd.mwave = 3;
    //
    // call RP solver for Maxwell equations
    //
    rp_maxwell(rd,ql_s,qr_s,df_s,wave_s,s_s,amdq_s,apdq_s);
    
    // restore values of meqn and mwave
    rd.meqn = meqn;
    rd.mwave = mwave;
}
