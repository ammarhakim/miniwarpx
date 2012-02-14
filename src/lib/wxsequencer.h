#ifndef __wxsequencer__h__
#define __wxsequencer__h__

// WarpX includes
#include "wxrange.h"

class WxSequencer
{
public:
    WxSequencer(const WxRange& r);
    ~WxSequencer();

    int step();
    int* indices() const;

private:
    WxRange _range;
    int _n, _len, _count;
    int *_indices;

    // no copying allowed
    WxSequencer(const WxSequencer&);
    WxSequencer& operator=(const WxSequencer&);
};

#endif // __wxsequencer__h__
