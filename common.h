/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * Copyright (c) 2004-2019 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

//#include "ompi_config.h"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

//#include "ompi/datatype/ompi_datatype.h"

//#include "ompi/mpi/c/bindings.h"
//#include "ompi/runtime/params.h"
//#include "ompi/communicator/communicator.h"
//#include "ompi/errhandler/errhandler.h"
//#include "opal/datatype/opal_convertor.h"
//#include "ompi/memchecker.h"

#if 0 && OPEN_MPI
extern void ompi_datatype_dump( MPI_Datatype ddt );
#define MPI_DDT_DUMP(ddt) ompi_datatype_dump( (ddt) )
#else
#define MPI_DDT_DUMP(ddt)
#endif  /* OPEN_MPI */

/* Find an approximation for the different levels of cache size. */
#if defined(__APPLE__)
#define L1size (64*1024)
#define L2size (4096*1024)
#define L3size 0
#define Lcachelinesize 64
#else
#define L1size sysconf(_SC_LEVEL1_DCACHE_SIZE)
#define L2size sysconf(_SC_LEVEL2_CACHE_SIZE)
#define L3size sysconf(_SC_LEVEL3_CACHE_SIZE)
#define Lcachelinesize 64
#endif

