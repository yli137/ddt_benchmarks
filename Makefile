TARGETS=alpha_normal_no_prefetch_pipeline alpha_normal_prefetch_pipeline alpha_prefetch_pipeline \
	alpha_pipeline_prefetch_pipeline alpha_pipeline_no_prefetch_pipeline

OMPI_SRC = /Users/bosilca/unstable/ompi/ompi/ompi/
OMPI_BIN = /Users/bosilca/unstable/ompi/ompi/build/fast
CC = mpicc
#CFLAGS = -O3 -I$(OMPI_SRC) -I$(OMPI_SRC)/opal/include -I$(OMPI_SRC)/ompi/include -I$(OMPI_BIN) -I$(OMPI_BIN)/opal/include
CFLAGS = -O3
LDFLAGS = -O3

all:$(TARGETS)

% : %.c common.h
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm -f $(TARGETS)
