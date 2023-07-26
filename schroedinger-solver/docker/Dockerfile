FROM ubuntu:20.04

RUN apt update && export DEBIAN_FRONTEND=noninteractive \
    && apt install -y --no-install-recommends g++ cppcheck \
    cmake make git ca-certificates unzip curl
