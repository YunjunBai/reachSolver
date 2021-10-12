# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.194.3/containers/ubuntu/.devcontainer/base.Dockerfile

# [Choice] Ubuntu version: hirsute, bionic, focal
ARG VARIANT="focal"
FROM mcr.microsoft.com/vscode/devcontainers/base:0-${VARIANT}

# [Optional] Uncomment this section to install additional OS packages.
# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends vim ca-certificates

RUN cp /etc/apt/sources.list /etc/apt/sources.list.bak \
    && sed -i "s/archive.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g" /etc/apt/sources.list \
    && sed -i "s/security.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g" /etc/apt/sources.list

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends build-essential cmake libgtest-dev libeigen3-dev \
    && ln -s /usr/include/eigen3/Eigen /usr/include/Eigen

RUN git clone https://github.com/catchorg/Catch2.git /tmp/Catch2\
    && cd /tmp/Catch2 && git checkout v2.x \
    && cmake -Bbuild -H. -DBUILD_TESTING=OFF \
    && cmake --build build/ --target install \
    && rm -r /tmp/Catch2

RUN git clone https://github.com/autodiff/autodiff /tmp/autodiff \
    && cd /tmp/autodiff && mkdir .build && cd .build \
    && cmake .. -DAUTODIFF_BUILD_PYTHON=0 \
    && cmake --build . --target install \
    && rm -r /tmp/autodiff
