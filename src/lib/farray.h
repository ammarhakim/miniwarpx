#ifndef __farray_cpp_2004__
#define __farray_cpp_2004__

// This file is for backward compatibility only. The FArray set of
// classes have now been replaced by new Wx... classed which should be
// used in the future.

// WarpX includes
#include "wxrange.h"
#include "wxindexer.h"
#include "wxarray.h"

// defines for compatibility: defines are used instead of typedefs as
// we can't have templated typedefs.
#define Range WxRange
#define Indexer WxIndexer
#define FArray WxArray

#endif // __farray_cpp_2004__
