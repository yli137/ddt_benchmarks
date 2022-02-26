#!/opt/local/bin/bash

declare -A mpis=( ["master"]="/Users/bosilca/opt/ompi/master/fast"
                  ["4.x"]="/Users/bosilca/opt/ompi/releases/4.1/fast"
                  ["5.x"]="/Users/bosilca/opt/ompi/releases/5.0/fast"
                  ["iovec"]="/Users/bosilca/opt/ompi/sandbox/yicheng/iovec/fast"
                  ["iovec_gather"]="/Users/bosilca/opt/ompi/sandbox/yicheng/iovec_gather/fast")
ARGS=""

for key in "${!mpis[@]}" ; do
  compiler="${mpis[$key]}/bin/mpicc"
  [[ -f "$compiler" && -x "$compiler" ]] && echo $compiler
  eval make check MPI=.${key} CC=${compiler} ARGS=${ARGS}
done