## Before build
    Use vscode-"Remote-Containers" or directly use .devcontainer/Dockerfile to build the dev docker.
    Then use src/autodiff/derivative.hpp to replace /usr/local/include/autodiff/forward/utils/derivative.hpp.
## Build Instructions:
    mkdir -p build && cd build
    cmake ../
    make

## Test Case:
    An incomplete case is in examples/test1.cpp(just show derivate functions). 