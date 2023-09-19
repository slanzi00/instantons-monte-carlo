# Monte Carlo simulation for Instantons

## Compilation Instructions

Follow these steps to compile the project:

1. **Docker**: If you haven't already, clone this repository to your local machine using Git:

  ```bash
  docker build - < docker/Dockerfile
  ```

choose the schroediger-solver or lattice-solution directory and build with cmake:

 ```bash
 cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```
