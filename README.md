## Prequires
    First install docker, see https://docs.docker.com/engine/install/ubuntu/ .
    Next run "sudo usermod -aG docker <username>" to let normal user use docker.
    Then install vscode, see https://code.visualstudio.com/docs/setup/linux .
    And install vscode-"Remote-Containers" plugins and build dev docker.
    
## Build Instructions:
    sudo cp src/autodiff/derivative.hpp /usr/local/include/autodiff/forward/utils/derivative.hpp
    mkdir -p build && cd build
    cmake ../
    make

## Test Case:
    An incomplete case is in examples/test1.cpp(just show derivate functions). 