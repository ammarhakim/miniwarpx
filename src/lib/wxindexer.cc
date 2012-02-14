#include "wxindexer.h"

WxIndexer::WxIndexer(const WxRange &r)
        : _rank(r.rank()), _ai(new int[r.rank()+1]), _r(r)
{
    // set a_1, ... a_N
    _ai[1] = 1;
    for(unsigned i=2; i<=_rank; ++i) 
        _ai[i] = _ai[i-1]*(r.end(i-2)-r.start(i-2)+1);
    
    // set a_0
    int sum = 0;
    for(unsigned i=1; i<=_rank; ++i) sum += _ai[i]*r.start(i-1);
    _ai[0] = -sum;
}

WxIndexer::WxIndexer(const WxIndexer& idx)
        : _rank(idx._rank), _ai(new int[idx._rank+1]), _r(idx._r)
{
    for(unsigned i=0; i<=_rank; ++i) _ai[i] = idx._ai[i];
}

WxIndexer::WxIndexer(unsigned rank, int *ai, const WxRange& r)
        : _rank(rank), _ai(new int[rank+1]), _r(r)
{
    for(unsigned i=0; i<=_rank; ++i) _ai[i] = ai[i];
}

WxIndexer::~WxIndexer()
{
    delete [] _ai;
}

WxIndexer&
WxIndexer::operator=(const WxIndexer& idx)
{
    if(this==&idx)
        return *this;
    delete [] _ai;
    
    _rank = idx._rank;
    _ai = new int[_rank+1];
    for(unsigned i=0; i<=_rank; ++i) _ai[i] = idx._ai[i];
    _r = idx._r;
    
    return *this;
}
