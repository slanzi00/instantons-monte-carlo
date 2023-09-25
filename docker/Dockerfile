FROM gcc:13

RUN apt update && export DEBIAN_FRONTEND=noninteractive \
    && apt install -y --no-install-recommends cppcheck \
    cmake make git ca-certificates unzip curl libboost-dev \
    libfmt-dev
