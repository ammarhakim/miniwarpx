#include "wxrange.h"


WxRange::WxRange()
        : _rank(1), _start(new int[1]), _end(new int[1]) 
{
    _start[0] = 0; _end[0] = 0;
}


WxRange::WxRange(unsigned rank)
        : _rank(rank), _start(new int[rank]), _end(new int[rank]) 
{
    for(unsigned i=0; i<rank; ++i)
    {
        _start[i] = 0; _end[i] = 0;
    }
}


WxRange::WxRange(unsigned rank, int *shape)
        : _rank(rank), _start(new int[rank]), _end(new int[rank]) 
{
    for(unsigned i=0; i<rank; ++i)
    {
        _start[i] = 0; _end[i] = shape[i]-1;
    }
}


WxRange::WxRange(unsigned rank, int *start, int *end)
        : _rank(rank), _start(new int[rank]), _end(new int[rank]) 
{
    for(unsigned i=0; i<rank; ++i)
    {
        _start[i] = start[i]; _end[i] = end[i];
    }
}


WxRange::WxRange(int s1, int e1)
        : _rank(1), _start(new int[1]), _end(new int[1]) 
{
    _start[0] = s1; _end[0] = e1;
}


WxRange::WxRange(int s1, int e1, int s2, int e2)
        : _rank(2), _start(new int[2]), _end(new int[2]) 
{
    _start[0] = s1; _end[0] = e1;
    _start[1] = s2; _end[1] = e2;
}


WxRange::WxRange(int s1, int e1, int s2, int e2, int s3, int e3)
        : _rank(3), _start(new int[3]), _end(new int[3]) 
{
    _start[0] = s1; _end[0] = e1;
    _start[1] = s2; _end[1] = e2;
    _start[2] = s3; _end[2] = e3;
}


WxRange::WxRange(int s1, int e1, int s2, int e2, int s3, int e3, int s4, int e4)
        : _rank(4), _start(new int[4]), _end(new int[4]) 
{
    _start[0] = s1; _end[0] = e1;
    _start[1] = s2; _end[1] = e2;
    _start[2] = s3; _end[2] = e3;
    _start[3] = s4; _end[3] = e4;
}


WxRange::WxRange(const WxRange& r)
{
    _rank   = r._rank;
    _start  = new int[_rank];
    _end    = new int[_rank];
    for(unsigned i=0; i<_rank; ++i)
    {
        _start[i]  = r._start[i];
        _end[i]    = r._end[i];
    }
}



WxRange&
WxRange::operator=(const WxRange& r)
{
    if(this==&r)
        return *this;
    delete [] _start;
    delete [] _end;

    _rank   = r._rank;
    _start  = new int[_rank];
    _end    = new int[_rank];
    for(unsigned i=0; i<_rank; ++i)
    {
        _start[i]  = r._start[i];
        _end[i]    = r._end[i];
    }
    return *this;
}

WxRange::~WxRange()
{
    delete [] _start;
    delete [] _end;
}
