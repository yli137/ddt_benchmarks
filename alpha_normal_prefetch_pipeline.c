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

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ompi/include/mpi.h"
#include "ompi/datatype/ompi_datatype.h"

#include "ompi_config.h"
#include "ompi/mpi/c/bindings.h"
#include "ompi/runtime/params.h"
#include "ompi/communicator/communicator.h"
#include "ompi/errhandler/errhandler.h"
#include "opal/datatype/opal_convertor.h"
#include "ompi/memchecker.h"

#if 0 && OPEN_MPI
extern void ompi_datatype_dump( MPI_Datatype ddt );
#define MPI_DDT_DUMP(ddt) ompi_datatype_dump( (ddt) )
#else
#define MPI_DDT_DUMP(ddt)
#endif  /* OPEN_MPI */

#define L1size sysconf(_SC_LEVEL1_DCACHE_SIZE)
#define L2size sysconf(_SC_LEVEL2_CACHE_SIZE)
#define L3size sysconf(_SC_LEVEL3_CACHE_SIZE)
void cache_flush(){
    char *cache = (char*)calloc(L1size+L2size+L3size, sizeof(char));
    free(cache);
}

static int trash_tlb_n = 8;
static int trash_disp[8] = {};
static int trash_len[8] = {8,8,8,8,8,8,8,8};
static int trash_dst_disp[8] = {0,8,16,24,32,40,48,56};

static int three_next_n = 2;
static int three_disp[2] = {0, 64};
static int three_len[2] = {24, 8};
static int three_dst_disp[2] = {0, 24};

static MPI_Datatype
create_random_indexed( int count, int blen, int stride, MPI_Datatype dt, int seed )
{
    MPI_Datatype ddt;
    int indices[count], block[count];

    size_t extent, lb;
    MPI_Type_get_extent( dt, &lb, &extent );

    srand(seed);
    indices[0] = 0;
    block[0] = blen;
    for( int i = 1; i < count; i++ ){

        indices[i] = i * stride * extent + rand() % ( (stride - blen) * extent );
        block[i] = blen;

    }

    MPI_Type_indexed( count, block, indices, dt, &ddt );
    MPI_Type_commit( &ddt );

    return ddt;
}

static MPI_Datatype
create_diagonal( int count )
{
    MPI_Datatype ddt;
    int indices[count], block[count];

    for( int i = 0; i < count; i++ ){
        indices[i] = i + i * count;
        block[i] = 1;
    }

    MPI_Type_indexed( count, block, indices, MPI_DOUBLE, &ddt );
    MPI_Type_commit( &ddt );

    return ddt;
}

static MPI_Datatype
create_upper_triangle( int count )
{
    MPI_Datatype ddt;
    int indices[count], block[count];

    for( int i = 0; i < count; i++ ){
        indices[i] = i + i * count;
        block[i] = count - i;
    }

    MPI_Type_indexed( count, block, indices, MPI_DOUBLE, &ddt );
    MPI_Type_commit( &ddt );

    return ddt;
}

static MPI_Datatype
create_lower_triangle( int count )
{
    MPI_Datatype ddt;
    int indices[count], block[count];

    for( int i = 0; i < count; i++ ){
        indices[i] = i * count;
        block[i] = i;
    }

    MPI_Type_indexed( count, block, indices, MPI_DOUBLE, &ddt );
    MPI_Type_commit( &ddt );

    return ddt;
}

static MPI_Datatype
create_merged_contig_with_gaps(int count)  /* count of the basic datatype */
{
    int array_of_blocklengths[] = {1, 1, 1};
    MPI_Aint array_of_displacements[] = {0, 8, 16};
    MPI_Datatype array_of_types[] = {MPI_DOUBLE, MPI_LONG, MPI_CHAR};
    MPI_Datatype type;

    MPI_Type_create_struct(3, array_of_blocklengths,
                           array_of_displacements, array_of_types,
                           &type);
    if( 1 < count ) {
        MPI_Datatype temp = type;
        MPI_Type_contiguous(count, temp, &type);
    }
    MPI_Type_commit(&type);
    MPI_DDT_DUMP( type );
    return type;
}

/* Create a non-contiguous resized datatype */
struct structure {
    double not_transfered;
    double transfered_1;
    double transfered_2;
};

static MPI_Datatype
create_struct_constant_gap_resized_ddt( int number,  /* IGNORED: number of repetitions */
                                        int contig_size,  /* IGNORED: number of elements in a contiguous chunk */
                                        int gap_size )    /* IGNORED: number of elements in a gap */
{
    struct structure data[1];
    MPI_Datatype struct_type, temp_type;
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    int blocklens[2] = {1, 1};
    MPI_Aint disps[3];

    MPI_Get_address(&data[0].transfered_1, &disps[0]);
    MPI_Get_address(&data[0].transfered_2, &disps[1]);
    MPI_Get_address(&data[0], &disps[2]);
    disps[1] -= disps[2]; /*  8 */
    disps[0] -= disps[2]; /* 16 */

    MPI_Type_create_struct(2, blocklens, disps, types, &temp_type);
    MPI_Type_create_resized(temp_type, 0, sizeof(data[0]), &struct_type);
    MPI_Type_commit(&struct_type);
    MPI_Type_free(&temp_type);
    MPI_DDT_DUMP( struct_type );

    return struct_type;
}

/* Create a datatype similar to the one use by HPL */
static MPI_Datatype
create_indexed_constant_gap_ddt( int number,  /* number of repetitions */
                                 int contig_size,  /* number of elements in a contiguous chunk */
                                 int gap_size )    /* number of elements in a gap */
{
    MPI_Datatype dt, *types;
    int i, *bLength;
    MPI_Aint* displ;

    types = (MPI_Datatype*)malloc( sizeof(MPI_Datatype) * number );
    bLength = (int*)malloc( sizeof(int) * number );
    displ = (MPI_Aint*)malloc( sizeof(MPI_Aint) * number );

    types[0] = MPI_DOUBLE;
    bLength[0] = contig_size;
    displ[0] = 0;
    for( i = 1; i < number; i++ ) {
        types[i] = MPI_DOUBLE;
        bLength[i] = contig_size;
        displ[i] = displ[i-1] + sizeof(double) * (contig_size + gap_size);
    }
    MPI_Type_create_struct( number, bLength, displ, types, &dt );
    MPI_DDT_DUMP( dt );
    free(types);
    free(bLength);
    free(displ);
    MPI_Type_commit( &dt );
    return dt;
}

static MPI_Datatype
create_optimized_indexed_constant_gap_ddt( int number,  /* number of repetitions */
                                           int contig_size,  /* number of elements in a contiguous chunk */
                                           int gap_size )    /* number of elements in a gap */
{
    MPI_Datatype dt;

    MPI_Type_vector( number, contig_size, (contig_size + gap_size), MPI_DOUBLE, &dt );
    MPI_Type_commit( &dt );
    MPI_DDT_DUMP( dt );
    return dt;
}

typedef struct {
   int i[2];
   float f;
} internal_struct;
typedef struct {
   int v1;
   int gap1;
   internal_struct is[3];
} ddt_gap;

static MPI_Datatype
create_indexed_gap_ddt( void )
{
    ddt_gap dt[2];
    MPI_Datatype dt1, dt2, dt3;
    int bLength[2] = { 2, 1 };
    MPI_Datatype types[2] = { MPI_INT, MPI_FLOAT };
    MPI_Aint displ[2];

    MPI_Get_address( &(dt[0].is[0].i[0]), &(displ[0]) );
    MPI_Get_address( &(dt[0].is[0].f), &(displ[1]) );
    displ[1] -= displ[0];
    displ[0] -= displ[0];
    MPI_Type_create_struct( 2, bLength, displ, types, &dt1 );
    /*MPI_DDT_DUMP( dt1 );*/
    MPI_Type_contiguous( 3, dt1, &dt2 );
    /*MPI_DDT_DUMP( dt2 );*/
    bLength[0] = 1;
    bLength[1] = 1;
    MPI_Get_address( &(dt[0].v1), &(displ[0]) );
    MPI_Get_address( &(dt[0].is[0]), &(displ[1]) );
    displ[1] -= displ[0];
    displ[0] -= displ[0];
    types[0] = MPI_INT;
    types[1] = dt2;
    MPI_Type_create_struct( 2, bLength, displ, types, &dt3 );
    /*MPI_DDT_DUMP( dt3 );*/
    MPI_Type_free( &dt1 );
    MPI_Type_free( &dt2 );
    MPI_Type_contiguous( 10, dt3, &dt1 );
    MPI_DDT_DUMP( dt1 );
    MPI_Type_free( &dt3 );
    MPI_Type_commit( &dt1 );
    return dt1;
}

static MPI_Datatype
create_indexed_gap_optimized_ddt( void )
{
    MPI_Datatype dt1, dt2, dt3;
    int bLength[3];
    MPI_Datatype types[3];
    MPI_Aint displ[3];

    MPI_Type_contiguous( 40, MPI_BYTE, &dt1 );
    MPI_Type_create_resized( dt1, 0, 44, &dt2 );

    bLength[0] = 4;
    bLength[1] = 9;
    bLength[2] = 36;

    types[0] = MPI_BYTE;
    types[1] = dt2;
    types[2] = MPI_BYTE;

    displ[0] = 0;
    displ[1] = 8;
    displ[2] = 44 * 9 + 8;

    MPI_Type_create_struct( 3, bLength, displ, types, &dt3 );

    MPI_Type_free( &dt1 );
    MPI_Type_free( &dt2 );
    MPI_DDT_DUMP( dt3 );
    MPI_Type_commit( &dt3 );
    return dt3;
}


/********************************************************************
 *******************************************************************/

#define DO_CONTIG                         0x00000001
#define DO_CONSTANT_GAP                   0x00000002
#define DO_INDEXED_GAP                    0x00000004
#define DO_OPTIMIZED_INDEXED_GAP          0x00000008
#define DO_STRUCT_CONSTANT_GAP_RESIZED    0x00000010
#define DO_STRUCT_MERGED_WITH_GAP_RESIZED 0x00000020

#define DO_PACK                         0x01000000
#define DO_UNPACK                       0x02000000
#define DO_ISEND_RECV                   0x04000000
#define DO_ISEND_IRECV                  0x08000000
#define DO_IRECV_SEND                   0x10000000
#define DO_IRECV_ISEND                  0x20000000

#define MIN_LENGTH   1024
#define MAX_LENGTH   1024*(1024*1024)

static int cycles  = 10;
static int trials  = 10;
static int warmups = 0;

static void print_result( size_t length, int trials, double* timers )
{
    double bandwidth, clock_prec, temp;
    double min_time, max_time, average, std_dev = 0.0;
    double ordered[trials];
    int t, pos, quartile_start, quartile_end;

    for( t = 0; t < trials; ordered[t] = timers[t], t++ );
    for( t = 0; t < trials-1; t++ ) {
        temp = ordered[t];
        pos = t;
        for( int i = t+1; i < trials; i++ ) {
            if( temp > ordered[i] ) {
                temp = ordered[i];
                pos = i;
            }
        }
        if( pos != t ) {
            temp = ordered[t];
            ordered[t] = ordered[pos];
            ordered[pos] = temp;
        }
    }
    quartile_start = trials - (3 * trials) / 4;
    quartile_end   = trials - (1 * trials) / 4;
    clock_prec = MPI_Wtick();
    min_time = ordered[quartile_start];
    max_time = ordered[quartile_start];
    average = ordered[quartile_start];
    for( t = quartile_start + 1; t < quartile_end; t++ ) {
        if( min_time > ordered[t] ) min_time = ordered[t];
        if( max_time < ordered[t] ) max_time = ordered[t];
        average += ordered[t];
    }
    average /= (quartile_end - quartile_start);
    for( t = quartile_start; t < quartile_end; t++ ) {
        std_dev += (ordered[t] - average) * (ordered[t] - average);
    }
    std_dev = sqrt( std_dev/(quartile_end - quartile_start) );
    
    bandwidth = (length * clock_prec) / (1024.0 * 1024.0) / (average * clock_prec);
    printf( "%zu\t%15g\t%10.4f MB/s [min %10g max %10g std %2.2f%%]\n", length, average, bandwidth,
            min_time, max_time, (100.0 * std_dev) / average );
}

static int pack( int cycles,
                 MPI_Datatype sdt, int scount, void* sbuf,
                 void* packed_buf )
{
    int position, myself, c, t, outsize;
    double timers[trials];

    MPI_Type_size( sdt, &outsize );
    outsize *= scount;

    MPI_Comm_rank( MPI_COMM_WORLD, &myself );

    for( t = 0; t < trials; t++ ) {
        cache_flush();
        timers[t] = MPI_Wtime();
        for( c = 0; c < cycles; c++ ) {
            position = 0;
            MPI_Pack(sbuf, scount, sdt, packed_buf, outsize, &position, MPI_COMM_WORLD);
        }
        timers[t] = (MPI_Wtime() - timers[t]) / cycles;
    }
    print_result( outsize, trials, timers );
    return 0;
}

static int do_pipeline_pack( const void *inbuf, int incount, MPI_Datatype datatype,
                             void *outbuf, int outsize, int pipe_size,
                             int num, int *disp, int *len, int *dst_disp )
{
    int position, myself, c, t, i, ddt_size;
    double timers[trials];
    size_t extent, lb;

    int hold_incount = incount;
    int pref_count = incount;
    
    int pc, check;
    size_t pref = 0;
    int j = 0, current_pos = 0, keep_pos = 0;
    int k;

    MPI_Type_get_extent( datatype, &lb, &extent );
    MPI_Type_size( datatype, &ddt_size );

    MPI_Comm_rank( MPI_COMM_WORLD, &myself );
    
    for( t = 0; t < trials; t++ ) {
        cache_flush();
        timers[t] = MPI_Wtime();

        for( c = 0; c < cycles; c++ ) {
            current_pos = 0;
            i = 0;

            for( ; i < pipe_size; i++ ){
                __builtin_prefetch( inbuf + i * extent );
            }

            for( j = 0 ; j < incount - pipe_size; j++, i++, current_pos++ ){
                __builtin_prefetch( inbuf + i * extent );
                for( k = 0; k < num; k++ ){
                    memcpy( outbuf + current_pos * ddt_size + dst_disp[k],
                            inbuf + current_pos * extent + disp[k],
                            len[k] );
                }
            }

            for( j = 0; j < pipe_size; j++, current_pos++ ){
                for( k = 0; k < num; k++ ){
                    memcpy( outbuf + current_pos * ddt_size + dst_disp[k],
                            inbuf + current_pos * extent + disp[k],
                            len[k] );
                }
            }

        }
        timers[t] = (MPI_Wtime() - timers[t]) / cycles;
    }
    print_result( outsize, trials, timers );
    return 0;
}

static int pack_pipeline( int cycles,
        MPI_Datatype sdt, int scount, void* sbuf,
        void* packed_buf, int pipe_size,
        int num, int *disp, int *len, int *dst_disp )
{
    int position, myself, c, t, outsize;
    double timers[trials];

    MPI_Type_size( sdt, &outsize );
    outsize *= scount;

    do_pipeline_pack(sbuf, scount, sdt, packed_buf, outsize, pipe_size,
                     num, disp, len, dst_disp);
    return 0;
}

static int unpack( int cycles,
        void* packed_buf,
        MPI_Datatype rdt, int rcount, void* rbuf )
{
    int position, myself, c, t, insize;
    double timers[trials];

    MPI_Type_size( rdt, &insize );
    insize *= rcount;

    MPI_Comm_rank( MPI_COMM_WORLD, &myself );

    for( t = 0; t < warmups; t++ ) {
        for( c = 0; c < cycles; c++ ) {
            position = 0;
            MPI_Unpack(packed_buf, insize, &position, rbuf, rcount, rdt, MPI_COMM_WORLD);
        }
    }

    for( t = 0; t < trials; t++ ) {
        cache_flush();
        timers[t] = MPI_Wtime();
        for( c = 0; c < cycles; c++ ) {
            position = 0;
            MPI_Unpack(packed_buf, insize, &position, rbuf, rcount, rdt, MPI_COMM_WORLD);
        }
        timers[t] = (MPI_Wtime() - timers[t]) / cycles;
    }
    print_result( insize, trials, timers );
    return 0;
}

static int do_test_for_ddt( int doop, MPI_Datatype sddt, MPI_Datatype rddt, size_t length )
{
    MPI_Aint lb, extent;
    char *sbuf, *rbuf;
    int i;
    int max_length;

    MPI_Type_get_extent( sddt, &lb, &extent );
    MPI_Type_size( sddt, &max_length );

    length = 100000000. / max_length * extent;

    sbuf = (char*)malloc( length );
    rbuf = (char*)malloc( length );


    if( doop & DO_PACK ) {
        printf("# Pack (max length %zu)\n", length);
        for( i = 1; i < (length / extent); i*=2  ) {
            pack( cycles, sddt, i, sbuf, rbuf );
        }
    }

    free( sbuf );
    free( rbuf );
    return 0;
}

static int do_pipeline_test_for_ddt( int doop, MPI_Datatype sddt, MPI_Datatype rddt, size_t length,
                                     int num, int *disp, int *len, int *dst_disp )
{
    MPI_Aint lb, extent;
    char *sbuf, *rbuf;
    int i;
    int max_length;

    MPI_Type_get_extent( sddt, &lb, &extent );
    MPI_Type_size( sddt, &max_length );

    length = 100000000. / max_length * extent;

    sbuf = (char*)malloc( length );
    rbuf = (char*)malloc( length );

    printf("\n# Pack (max length %zu) Pipeline FULL ddt per segment\n",
            length,
            length / extent);

    for( i = 1; i < (length / extent); i*=2 ){
        pack_pipeline( cycles, sddt, i, sbuf, rbuf, i,
                num, disp, len, dst_disp );
    }

    for( int j = 1; j < 513; j *= 2 ){
        printf("\n# Pack (max length %zu) Pipeline %d ddt per segment\n", 
                length,
                j);
        for( i = j; i < (length / extent); i*=2  ) {
            pack_pipeline( cycles, sddt, i, sbuf, rbuf, j,
                           num, disp, len, dst_disp );
        }
    }

    free( sbuf );
    free( rbuf );
    return 0;
}

static int matrix_pack( int doop, MPI_Datatype sddt, MPI_Datatype rddt, size_t length )
{
    MPI_Aint lb, extent;
    char *sbuf, *rbuf;
    int i;
    int max_length;

    MPI_Type_get_extent( sddt, &lb, &extent );
    MPI_Type_size( sddt, &max_length );

    if( extent > max_length )
        length = extent;
    else
        length = max_length;

    sbuf = (char*)malloc( length );
    rbuf = (char*)malloc( length );

    pack( cycles, sddt, 1, sbuf, rbuf );

    free( sbuf );
    free( rbuf );
    return 0;
}

static int matrix_unpack( int doop, MPI_Datatype sddt, MPI_Datatype rddt, size_t length )
{
    MPI_Aint lb, extent;
    char *sbuf, *rbuf;
    int i;
    int max_length;

    MPI_Type_get_extent( sddt, &lb, &extent );
    MPI_Type_size( sddt, &max_length );

    if( extent > max_length )
        length = extent;
    else
        length = max_length;

    sbuf = (char*)malloc( length );
    rbuf = (char*)malloc( length );

    unpack( cycles, sbuf, rddt, 1, rbuf );

    free( sbuf );
    free( rbuf );
    return 0;
}

int main( int argc, char* argv[] )
{
    int run_tests = 0xffff;  /* do all datatype tests by default */
    int rank, size;
    MPI_Datatype ddt;

    run_tests |= DO_PACK | DO_UNPACK;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    if( rank != 0 ) {
        MPI_Finalize();
        exit(0);
    }

    printf("\n! 3 doubles 1 next cache line\n\n");
    int blen2[] = { 3, 1 };
    int disp2[] = { 0, 8 };
    MPI_Type_indexed( 2, blen2, disp2, MPI_DOUBLE, &ddt );
    MPI_Type_create_resized( ddt, 0, 64 * 8, &ddt );
    MPI_Type_commit( &ddt );

    do_test_for_ddt( run_tests, ddt, ddt, MAX_LENGTH );
    do_pipeline_test_for_ddt( run_tests, ddt, ddt, MAX_LENGTH,
                              three_next_n, &(three_disp[0]),
                              &(three_len[0]), &(three_dst_disp[0]) );
    MPI_Type_free( &ddt );

    printf("\n! Trashing TLB datatype\n\n");
    int disp5[] = { 0, 4096, 8, 4104, 16, 4112, 24, 4120 };
    int blen5[] = { 1, 1, 1, 1, 1, 1, 1, 1 };
    MPI_Type_indexed( 8, blen5, disp5, MPI_DOUBLE, &ddt );
    MPI_Type_create_resized( ddt, 0, 8192, &ddt );
    MPI_Type_commit( &ddt );

    do_test_for_ddt( run_tests, ddt, ddt, MAX_LENGTH );
    do_pipeline_test_for_ddt( run_tests, ddt, ddt, MAX_LENGTH,
                              trash_tlb_n, &(trash_disp[0]),
                              &(trash_len[0]), &(trash_dst_disp[0]) );
    MPI_Type_free( &ddt );

    MPI_Finalize();
    exit(0);
}

