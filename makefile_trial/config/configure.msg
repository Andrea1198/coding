--------------------------------------------------------------------
ESPRESSO can take advantage of several optimized numerical libraries
(essl, fftw, mkl...).  This configure script attempts to find them,
but may fail if they have been installed in non-standard locations.
If a required library is not found, the local copy will be compiled.

The following libraries have been found:
  BLAS_LIBS=$(TOPDIR)/extlibs/blas/lib/libblas.a 
  LAPACK_LIBS=$(TOPDIR)/extlibs/lapack/lib/liblapack.a
  @delete@
  @delete@
  FFT_LIBS=
  @delete@
  
  @delete@
Please check if this is what you expect.

If any libraries are missing, you may specify a list of directories
to search and retry, as follows:
  ./configure LIBDIRS="list of directories, separated by spaces"

Parallel environment not detected \(is this a parallel machine?\).\
Configured for compilation of serial executables.

For more info, read the ESPRESSO User's Guide (Doc/users-guide.tex).
--------------------------------------------------------------------
