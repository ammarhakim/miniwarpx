#ifndef __wxarray__h__
#define __wxarray__h__

// WarpX includes
#include "wxrange.h"
#include "wxindexer.h"
#include "wxsequencer.h"

// std includes
#include <iostream>
#include <cassert>

// following macros are used to manipulate WxArray.traits
static unsigned int masks[] = 
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// is array contiguous?
#define CONTIGUOUS masks[0]
#define SET_CONTIGUOUS(bit) (bit) |= CONTIGUOUS
#define CLEAR_CONTIGUOUS(bit) (bit) &= ~CONTIGUOUS
#define IS_CONTIGUOUS(bit) (bit) & CONTIGUOUS

// was array data allocated with a malloc?
#define ALLOC masks[1]
#define SET_ALLOC(bit) (bit) |= ALLOC
#define CLEAR_ALLOC(bit) (bit) &= ~ALLOC
#define IS_ALLOC(bit) (bit) & ALLOC

/**
   WxArray provides a simple yet powerful N-dimensional array
   facility. The array start index can be arbitrary. The arrays are
   reference counted and can hence be passed and returned from
   functions without worrying about data copying. The arrays can be
   sliced or reshaped and can be created from raw C pointers. By
   deafult the data is stored contiguously in Fortran style, i.e in
   column major order. Thus, in most cases the arrays can be passed to
   Fortran routines. However, in certain cases when arrays are sliced
   its elements may no longer be stored contiguously and the array
   data can not be then passed to Fortran.
*/
template<typename T,
         typename WXINDEXER = WxIndexer>
class WxArray
{
public:

    /** 
        default ctor is meant to be used as a place holder: beware of
        using WxArray object returned by this constructor
    */
    WxArray();

    //! 'r' is range object which specifies array dimensions.
    WxArray(const WxRange& r);

    /**
       'r' is range object which specifies array dimensions. 'data' is
       initial value assigned to all array elements
    */
    WxArray(const WxRange& r, const T& data);

    /**
       'r' is range object which specifies array dimensions. Here, and
       in the following ctors 'cdata' is a pointer to a raw array of
       type T. Such constructor do not allocate any memory but simply
       use the space provided by 'cdata'. In general this is not safe
       as there is no way of knowing if the size of the 'cdata' array
       is at least as large as r.size(). However, this feature is very
       useful to "wrap" a raw C array so that it behaves like a
       WxArray.
    */
    WxArray(const WxRange& r, T* cdata);

    //! 'rank' is rank of array, 'dims[k]' is size along kth dimension
    WxArray(unsigned rank, int *dims);

    /**
       'rank' is rank of array, 'dims[k]' is size along kth dimension
       and 'data' is initial value assigned to all array elements
    */
    WxArray(unsigned rank, int *dims, const T& data);
    WxArray(unsigned rank, int *dims, T* cdata);

    /**
       Copy construtor does not allocate new memmory. The created
       array shares the data with the original one.
    */
    WxArray(const WxArray& f);
    /**
       Assignment opertor does not allocate new memmory. The created
       array shares the data with the original one.
    */
    WxArray& operator=(const WxArray& f);

    ~WxArray();

    //! range object of array
    const WxRange& range() const; 
    //! rank
    unsigned rank() const; // rank
    //! size along kth dimension
    int dims(unsigned k) const;
    //! start index along kth dimension
    int start(unsigned k) const;
    //! end index along kth dimension
    int end(unsigned k) const;
    //! total no of elements in array
    int size() const;

    /**
       actual pointer to data: modifying this pointer directly can be
       rather dangerous
    */
    T* data();
    T const* data() const;

    /**
       general indexing routine. 'indices[k]', k=0.._r.rank()-1 is the
       index into the kth dimension.
    */
    T operator()(int *indices) const;
    /**
       general indexing routine. 'indices[k]', k=0.._r.rank()-1 is the
       index into the kth dimension.
    */
    T& operator()(int *indices);

    T index(int *indices) const;
    T& index(int *indices);

    /**
     * Linear index into array
     */
    int linloc(int *indices) const;

    /**
       indexers for ranks 1 to 4. In the following kp is the index
       into the pth dimension. The index() functions are provided as
       operator overloading syntax for pointers is rather ugly
       constrast f->operator()(i,j); with f->index(i,j);
    */

    //! rank-1 indexer
    T operator()(int k1) const;
    //! rank-1 indexer
    T& operator()(int k1);
    //! rank-1 indexer
    T index(int k1) const;
    //! rank-1 indexer
    T& index(int k1);

    //! rank-2 indexer
    T operator()(int k1, int k2) const;
    //! rank-2 indexer
    T& operator()(int k1, int k2);
    //! rank-2 indexer
    T index(int k1, int k2) const;
    //! rank-2 indexer
    T& index(int k1, int k2);

    //! rank-3 indexer
    T operator()(int k1, int k2, int k3) const;
    //! rank-3 indexer
    T& operator()(int k1, int k2, int k3);
    //! rank-3 indexer
    T index(int k1, int k2, int k3) const;
    //! rank-3 indexer
    T& index(int k1, int k2, int k3);

    //! rank-4 indexer
    T operator()(int k1, int k2, int k3, int k4) const;
    //! rank-4 indexer
    T& operator()(int k1, int k2, int k3, int k4);
    //! rank-4 indexer
    T index(int k1, int k2, int k3, int k4) const;
    //! rank-4 indexer
    T& index(int k1, int k2, int k3, int k4);


    /**
       Slices array. The returned array shares data with the
       original. The returned array is not usually contiguous so the
       raw C array returned by its 'data()' function should NOT be
       used (without verifing with is_contiguous()) UNDER ANY
       CIRCUMSTANCES. For example, after slicing an WxArray the
       returned slice can not be usually passed to a Fortran
       routine. You have been warned!
    */
    WxArray<T,WXINDEXER> slice(const WxRange& r);

    /**
       This version of slice returns an array with starting indices
       start[0]...start[_range-1]
    */
    WxArray<T,WXINDEXER> slice(const WxRange&r, int start[]);
    
    // functions for querying array traits

    //! has array allocated any data?
    bool is_alloc() const;
    //! is array pointing to contiguous piece of memory?
    bool is_contiguous() const;

private:
    WxRange _range; // rank of array
    T *_data; // actual data stored
    unsigned *_count; // no of references to this array
    WXINDEXER _indexer; // indexing object
    unsigned _traits; // array traits
};

template<typename T, typename WXINDEXER>
inline
const WxRange&
WxArray<T,WXINDEXER>::range() const { return _range; }

template<typename T, typename WXINDEXER>
inline
unsigned
WxArray<T,WXINDEXER>::rank() const { return _range.rank(); }

template<typename T, typename WXINDEXER>
inline
int
WxArray<T,WXINDEXER>::dims(unsigned k) const { return _range.shape(k); }

template<typename T, typename WXINDEXER>
inline
int
WxArray<T,WXINDEXER>::start(unsigned k) const { return _range.start(k); }

template<typename T, typename WXINDEXER>
inline
int
WxArray<T,WXINDEXER>::end(unsigned k) const { return _range.end(k); }

template<typename T, typename WXINDEXER>
inline
int
WxArray<T,WXINDEXER>::size() const { return _range.size(); }

template<typename T, typename WXINDEXER>
inline
T*
WxArray<T,WXINDEXER>::data() { return _data; }

template<typename T, typename WXINDEXER>
inline
T const*
WxArray<T,WXINDEXER>::data() const { return _data; }

//
// ctors
//

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray()
        : _data(0), _count(new unsigned(0)), _indexer(WxRange()), _traits(0)
{
    CLEAR_ALLOC(_traits);
    CLEAR_CONTIGUOUS(_traits);
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(const WxRange& r)
        : _range(r), _data(new T[r.size()]), _count(new unsigned(1)), 
          _indexer(r), _traits(0)
{
    SET_ALLOC(_traits);
    SET_CONTIGUOUS(_traits);
    for(unsigned i=0; i<_range.size(); ++i) _data[i] = 0;
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(const WxRange& r, const T& data)
        : _range(r), _data(new T[r.size()]), _count(new unsigned(1)),
          _indexer(r), _traits(0)
{
    SET_ALLOC(_traits);
    SET_CONTIGUOUS(_traits);
    for(unsigned i=0; i<_range.size(); ++i) _data[i] = data;
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(const WxRange& r, T* cdata)
        : _range(r), _data(cdata), _count(new unsigned(1)), 
          _indexer(r), _traits(0)
{
    CLEAR_ALLOC(_traits);
    SET_CONTIGUOUS(_traits);
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(unsigned rank, int *dims)
        : _range(rank, dims), _count(new unsigned(1)), 
          _indexer(_range), _traits(0)
{
    SET_ALLOC(_traits);
    SET_CONTIGUOUS(_traits);
    _data = new T[_range.size()];
    for(unsigned i=0; i<_range.size(); ++i) _data[i] = 0;
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(unsigned rank, int *dims, const T& data)
        : _range(rank, dims), _count(new unsigned(1)), 
          _indexer(_range), _traits(0)
{
    SET_ALLOC(_traits);
    SET_CONTIGUOUS(_traits);
    _data = new T[_range.size()];
    for(unsigned i=0; i<_range.size(); ++i) _data[i] = data;
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(unsigned rank, int *dims, T* cdata)
        : _range(rank, dims), _data(cdata), _count(new unsigned(1)),
          _indexer(_range), _traits(0)
{
    CLEAR_ALLOC(_traits);
    SET_CONTIGUOUS(_traits);
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::WxArray(const WxArray& f)
        : _range(f._range), _data(f._data), _indexer(f._indexer), _traits(f._traits)
{
    CLEAR_ALLOC(_traits); // memory not allocated
    ++*f._count;
    _count = f._count;
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>&
WxArray<T,WXINDEXER>::operator=(const WxArray &f)
{
    ++*f._count; // increase use count of f
    if(--*_count == 0) // reduce own use count
    { // no other array was assigned to this array
        if(IS_ALLOC(_traits))
        { // delete data only if it was allocated
            delete [] _data;
            delete _count;
        }
    }
    _range = f._range;
    _indexer = f._indexer;
    _data = f._data;
    _count = f._count;
    _traits = f._traits;
    CLEAR_ALLOC(_traits); // memory not allocated
    return *this;
}

// dtor: reduce count by one and delete if needed
template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER>::~WxArray()
{
    if((--*_count == 0) && IS_ALLOC(_traits))
    {
        delete [] _data;
        delete _count;
    }
}

template<typename T, typename WXINDEXER>
inline
T
WxArray<T,WXINDEXER>::operator()(int *indices) const
{
    return _data[_indexer.index(indices)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::operator()(int *indices)
{
    return _data[_indexer.index(indices)]; 
}

template<typename T, typename WXINDEXER>
inline
T
WxArray<T,WXINDEXER>::operator()(int k1) const
{
    return _data[_indexer.index(k1)]; 
}

template<typename T, typename WXINDEXER>
inline
T&
WxArray<T,WXINDEXER>::operator()(int k1)
{ 
    return _data[_indexer.index(k1)]; 
}

template<typename T, typename WXINDEXER>
inline
T 
WxArray<T,WXINDEXER>::operator()(int k1, int k2) const
{
    return _data[_indexer.index(k1,k2)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::operator()(int k1, int k2)
{
    return _data[_indexer.index(k1,k2)]; 
}

template<typename T, typename WXINDEXER>
inline
T 
WxArray<T,WXINDEXER>::operator()(int k1, int k2, int k3) const
{
    return _data[_indexer.index(k1,k2,k3)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::operator()(int k1, int k2, int k3)
{
    return _data[_indexer.index(k1,k2,k3)]; 
}

template<typename T, typename WXINDEXER>
inline
T 
WxArray<T,WXINDEXER>::operator()(int k1, int k2, int k3, int k4) const
{
    return _data[_indexer.index(k1,k2,k3,k4)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::operator()(int k1, int k2, int k3, int k4)
{
    return _data[_indexer.index(k1,k2,k3,k4)]; 
}

template<typename T, typename WXINDEXER>
inline
T
WxArray<T,WXINDEXER>::index(int *indices) const
{
    return _data[_indexer.index(indices)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::index(int *indices)
{
    return _data[_indexer.index(indices)]; 
}

template<typename T, typename WXINDEXER>
inline
T
WxArray<T,WXINDEXER>::index(int k1) const
{
    return _data[_indexer.index(k1)]; 
}

template<typename T, typename WXINDEXER>
inline
T&
WxArray<T,WXINDEXER>::index(int k1)
{ 
    return _data[_indexer.index(k1)]; 
}

template<typename T, typename WXINDEXER>
inline
T 
WxArray<T,WXINDEXER>::index(int k1, int k2) const
{
    return _data[_indexer.index(k1,k2)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::index(int k1, int k2)
{
    return _data[_indexer.index(k1,k2)]; 
}

template<typename T, typename WXINDEXER>
inline
T 
WxArray<T,WXINDEXER>::index(int k1, int k2, int k3) const
{
    return _data[_indexer.index(k1,k2,k3)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::index(int k1, int k2, int k3)
{
    return _data[_indexer.index(k1,k2,k3)]; 
}

template<typename T, typename WXINDEXER>
inline
T 
WxArray<T,WXINDEXER>::index(int k1, int k2, int k3, int k4) const
{
    return _data[_indexer.index(k1,k2,k3,k4)]; 
}

template<typename T, typename WXINDEXER>
inline
T& 
WxArray<T,WXINDEXER>::index(int k1, int k2, int k3, int k4)
{
    return _data[_indexer.index(k1,k2,k3,k4)]; 
}

template<typename T, typename WXINDEXER>
inline
int
WxArray<T,WXINDEXER>::linloc(int *indices) const
{
    return _indexer.index(indices); 
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER> 
WxArray<T,WXINDEXER>::slice(const WxRange& r)
{
#ifdef _DO_RANGE_CHECK_
    assert(r.rank()==rank());
    for(unsigned i=0; i<r.rank(); ++i)
        assert(r.start(i)>=start(i) && r.end(i)<=end(i));
#endif
    int ndims = rank();
    int *tmp = new int[ndims+1];
    // compute dimensions of the array
    for(unsigned i=0; i<ndims; ++i)
        tmp[i] = r.end(i)-r.start(i)+1;

    WxArray<T,WXINDEXER> f; // allocate empty array

    *f._count = 1; // set use count
    f._range = WxRange(ndims, tmp); // set new range object

    // compute starting index in slice
    for(unsigned i=0; i<ndims; ++i)
        tmp[i] = r.start(i); 
    int loc = _indexer.index(tmp); // get its location in _data

    // compute indices
    tmp[0] = 0;
    for(unsigned i=1; i<ndims+1; ++i)
        tmp[i] = _indexer._ai[i];
    f._indexer = WxIndexer(ndims, tmp, f._range); // set new indexer object
    f._data = _data + loc; // set data array

    // set array traits
    CLEAR_ALLOC(f._traits);
    CLEAR_CONTIGUOUS(f._traits);

    delete [] tmp;

    return f;
    
}

template<typename T, typename WXINDEXER>
WxArray<T,WXINDEXER> 
WxArray<T,WXINDEXER>::slice(const WxRange& r, int beg[])
{
#ifdef _DO_RANGE_CHECK_
    assert(r.rank()==rank());
    for(unsigned i=0; i<r.rank(); ++i)
        assert(r.start(i)>=start(i) && r.end(i)<=end(i));
#endif
    unsigned ndims = rank();
    int *tmp = new int[ndims+1];
    // compute end indices of the array
    for(unsigned i=0; i<ndims; ++i)
        tmp[i] = beg[i]+r.end(i)-r.start(i);

    WxArray<T,WXINDEXER> f; // allocate empty array

    *f._count = 1; // set use count
    f._range = WxRange(ndims, beg, tmp); // set new range object

    // compute starting index in slice
    for(unsigned i=0; i<ndims; ++i)
        tmp[i] = r.start(i); 
    int loc = _indexer.index(tmp); // get its location in _data

    // compute indices
    tmp[0] = 0;
    for(unsigned i=1; i<ndims+1; ++i)
        tmp[i] = _indexer._ai[i];
    f._indexer = WxIndexer(ndims, tmp, f._range); // set new indexer object
    // compute starting location in f._data
    for(unsigned i=0; i<ndims; ++i)
        tmp[i] = beg[i];
    int locf = f._indexer.index(tmp);
    f._data = _data + loc - locf; // set data array

    // set array traits
    CLEAR_ALLOC(f._traits);
    CLEAR_CONTIGUOUS(f._traits);

    delete [] tmp;

    return f;
    
}

template<typename T, typename WXINDEXER>
bool 
WxArray<T,WXINDEXER>::is_alloc() const { return IS_ALLOC(_traits); }

template<typename T, typename WXINDEXER>
bool 
WxArray<T,WXINDEXER>::is_contiguous() const { return IS_CONTIGUOUS(_traits); }

#endif // __wxarray__h__
