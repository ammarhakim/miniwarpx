#ifndef __wxrange__h__
#define __wxrange__h__

#include <cassert>

/**
   WxRange represents a slice of an n-dimensional space of integers.
*/
class WxRange
{
public:
    /**
       Default ctor is meant to be used as a place holder: beware of
       using WxRange object returned by this constructor
    */
    WxRange();

    //! 'rank' is rank of array 
    WxRange(unsigned rank);;

    /**
       'rank' is rank of array. All start indices are assumed to be 0
       while kth end index is shape[k]-1, k=0..rank-1
    */
    WxRange(unsigned rank, int *shape);
    
    /**
       'rank' is rank of array, 'start[k]', 'end[k]', are start index,
       end index along kth dimension, for k=0..rank-1
    */
    WxRange(unsigned rank, int *start, int *end);

    // constructors for ranks 1 to 4. In the following sk and ek are
    // start, and end for kth dimension
    
    //! rank-1 ctor
    WxRange(int s1, int e1);

    //! rank-2 ctor
    WxRange(int s1, int e1, int s2, int e2);

    //! rank-3 ctor
    WxRange(int s1, int e1, int s2, int e2, int s3, int e3);

    //! rank-4 ctor
    WxRange(int s1, int e1, int s2, int e2, int s3, int e3, int s4, int e4);

    WxRange(const WxRange& r);
    WxRange& operator=(const WxRange& r);

    ~WxRange();

    //! rank
    unsigned  rank() const { 
        return _rank; 
    }
    //! start index along dimension k
    int  start(unsigned k) const {
#ifdef _DO_RANGE_CHECK_
        assert(k<_rank);
#endif
        return _start[k]; 
    }

    //! end index along dimension k
    int  end(unsigned k) const {
#ifdef _DO_RANGE_CHECK_
        assert(k<_rank);
#endif
        return _end[k]; 
    }

    //! elements along dimension k
    int shape(unsigned k) const {
#ifdef _DO_RANGE_CHECK_
        assert(k<_rank);
#endif 
        return _end[k]-_start[k]+1;
    }

    //! total no of elements in array
    unsigned size() const {
        unsigned _size = 1;
        for(unsigned i=0; i<_rank; ++i)
            _size *= _end[i]-_start[i]+1;
        return _size;
    }

private:
    unsigned _rank;
    int *_start, *_end;
};

#endif // __wxrange__h__
