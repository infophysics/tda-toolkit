# bottleneck

C++ code for calculating bottleneck distance between two [peristence diagrams](https://en.wikipedia.org/wiki/Topological_data_analysis) with mex code for Matlab interop

Based on the [Hopcroft-Karp algorithm](https://en.wikipedia.org/wiki/Hopcroft%E2%80%93Karp_algorithm)


### Compiling

Depends on `g++` compiler

```sh
make
```

#### Building for Matlab

From the Matlab command line, run

```sh
cd <THIS_DIRECTORY>
mex bottleneck_distance.cpp
```

Then from your Matlab code you can call

```matlab
bottleneck_distance(mat1, mat2)
```

where `mat1` and `mat2` are two `Mx2` matrices representing persistence diagrams.
The function outputs the distance between the two matrices.

#### Running on the command line

You can also run from the command line if you have saved the persistence diagrams to files

```sh
./bottleneck pers_dia1.txt pers_dia2.txt
```


### Testing

run `test_bottleneck.sh`

```sh
./test_bottleneck.sh
```

### Authors

Written by Miro Kramer (2013), modified/optimized by Kelvin Abrokwa (2016)

