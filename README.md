# Monte Carlo simulation for Instantons

## Compilation Instructions

Follow these steps to compile the project:

1. **Docker**:

  ```bash
  docker build - < docker/Dockerfile
  ```

choose the schroediger-solver or lattice-solution directory and build with cmake:

 ```bash
 cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```
