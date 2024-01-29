# Monte Carlo simulation for Instantons

## Compilation Instructions

Follow these steps to compile the project:

1. Install libraries and run the code with CMake:

```cppcheck``` and ```libboost-dev``` (Boost C++ libraries) are required. Choose which of the two subproject
you want build  (```schroediger-solver``` or ```lattice-solution```) and build:
 ```bash
 cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
 build/scrouedinger-solver(or build/lattice-solutions)
```

2. Using **Docker**: 

build the docker image:
 ```bash
docker build -t instantons_docker docker
```
open the image and include the directories (example of command):
 ```bash
docker run -it -v $(pwd):/workspace instantons_docker docker
# insert the full name of the container
```
then choose project and build with CMake
 ```bash
 cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```
