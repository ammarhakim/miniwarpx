#ifndef __wxindexer__h__
#define __wxindexer__h__

#include <cassert>
#include "wxrange.h"

template<typename T, typename WXINDEXER> class WxArray;

/**
   WxIndexer provides a way to map an element of an n-dimensional
   index set into a single integer. This particular indexer maps an
   index set using a column major order.
 */
class WxIndexer
{
public:
    // this chummy-ness is presently required to make array slicing
    // work. MUST FIX THIS ASAP.
    template <typename T, typename WXINDEXER> friend class WxArray;

    WxIndexer(const WxRange& r);
    WxIndexer(unsigned rank, int *ai, const WxRange& r);
    ~WxIndexer();

    WxIndexer(const WxIndexer& idx);
    WxIndexer& operator=(const WxIndexer& idx);

    //! rank of indexer
    unsigned rank() const {
        return _rank;
    }

    //! range of indexer
    const WxRange& range() const {
        return _r;
    }

    //! Index 1D array
    int index(int k1) const {
#ifdef _DO_RANGE_CHECK_
        assert(_rank==1);
        assert((k1>=_r.start(0)) && (k1<=_r.end(0)));
#endif 
        return _ai[0]+_ai[1]*k1;
    }

    //! Index 2D array
    int index(int k1, int k2) const {
#ifdef _DO_RANGE_CHECK_
        assert(_rank==2);
        assert((k1>=_r.start(0)) && (k1<=_r.end(0)));
        assert((k2>=_r.start(1)) && (k2<=_r.end(1)));
#endif 
        return _ai[0]+_ai[1]*k1+_ai[2]*k2;         
    }

    //! Index 3D array
    int index(int k1, int k2, int k3) const {
#ifdef _DO_RANGE_CHECK_
        assert(_rank==3);
        assert((k1>=_r.start(0)) && (k1<=_r.end(0)));
        assert((k2>=_r.start(1)) && (k2<=_r.end(1)));
        assert((k3>=_r.start(2)) && (k3<=_r.end(2)));
#endif 
        return _ai[0]+_ai[1]*k1+_ai[2]*k2+_ai[3]*k3;         
    }

    //! Index 4D array
    int index(int k1, int k2, int k3, int k4) const {
#ifdef _DO_RANGE_CHECK_
        assert(_rank==4);
        assert((k1>=_r.start(0)) && (k1<=_r.end(0)));
        assert((k2>=_r.start(1)) && (k2<=_r.end(1)));
        assert((k3>=_r.start(2)) && (k3<=_r.end(2)));
        assert((k4>=_r.start(3)) && (k4<=_r.end(3)));
#endif  
        return _ai[0]+_ai[1]*k1+_ai[2]*k2+_ai[3]*k3+_ai[4]*k4;        
    }

    //! Index arbitrary dimensional array
    int index(int *k) const {
#ifdef _DO_RANGE_CHECK_
    for(unsigned i=0; i<rank(); ++i)
        assert(k[i]>=_r.start(i) && k[i]<=_r.end(i));
#endif
        int sum=_ai[0];
        for(unsigned i=1; i<=_rank; ++i) sum += _ai[i]*k[i-1];
        return sum;
    }

protected:
    unsigned _rank;
    int *_ai;
    WxRange _r;
};

#endif // __wxindexer__h__
