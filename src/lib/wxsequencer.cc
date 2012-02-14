#include "wxsequencer.h"

WxSequencer::WxSequencer(const WxRange& r)
        : _range(r), _n(r.rank()), _len(1), _count(0), _indices(new int[r.rank()])
{
    int length;
    for(int i=0; i<_n; ++i)
    {
        length = r.end(i)-r.start(i)+1;
        _indices[i] = r.start(i);
        _len *= length;
    }
}

WxSequencer::~WxSequencer()
{
    delete [] _indices;
}

int
WxSequencer::step()
{
    int i,j;
    
    if(_count++==0)
        return 1;
    // further calls
    // check if no more indices
    if(_count > _len)
        return 0;

    _indices[0] += 1;
    if(_indices[0] > _range.end(0))
    {
        // run out of indices along 0th dimension:
        // check for closest indices to increase
        for(i = 1; i < _n; ++i)
        {
            _indices[i] += 1;
            if(!(_indices[i] > _range.end(i)))
            {
                // reset old
                for(j = 0; j < i; j++) _indices[j] = _range.start(j);
                return 1;
            }
        }
    }
    return 1;    
}

int*
WxSequencer::indices() const
{ return _indices; }
