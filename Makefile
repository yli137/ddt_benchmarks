TARGETS=alpha_normal_prefetch_pipeline
#TARGETS+=alpha_normal_no_prefetch_pipeline  alpha_prefetch_pipeline \
	alpha_pipeline_prefetch_pipeline alpha_pipeline_no_prefetch_pipeline

INTERNAL_TARGETS=$(addsuffix ${MPI},$(TARGETS))

ifeq ($(MPI),)
MPI:=
endif

ifeq ($(ARGS),)
ARGS:=
endif

ifeq "$(origin CC)" "default"
CC := mpicc
endif

CFLAGS = -O3
LDFLAGS = -O3

all:$(INTERNAL_TARGETS)

%${MPI} : %.c common.h
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

check: all
	$(foreach target,$(INTERNAL_TARGETS),./$(target) $(ARGS);)

clean:
	rm -rf $(INTERNAL_TARGETS) *.o $(addsuffix .dSYM,$(TARGETS))
