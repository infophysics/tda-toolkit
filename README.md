
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1436034.svg)](https://doi.org/10.5281/zenodo.1436034)


# TDA
Topological Data Analysis Toolkit


Implementation of Persistent Homology algorithms, mostly by wrapping existing algorithms in python.  The following algorithms are wrapped;

- RIPSER; https://github.com/Ripser
U. Bauer (2016)

- Cubical RIPSER; https://github.com/CubicalRipser
Takeki Sudo and Kazushi Ahara (2017)

- PERSEUS; http://people.maths.ox.ac.uk/nanda/perseus/index.html
Software - V. Nanda
Paper - "Morse Theory for Filtrations and Efficient Computation of Persistent Homology", V. Nanda and K. Mischaikow (http://people.maths.ox.ac.uk/nanda/source/DMPers-Final.pdf)


Implementation of an algorithm to compute the Persistent Homology Dimension associated to a barcode diagram based on the paper; "Measuring Shape with Topology" by Robert MacPherson and Benjamin Schweinhart - arXiv:1011.2258 


Implementation of a Bottleneck Distance algorithm written by Miro Kramer (2013), modified/optimized by Kelvin Abrokwa (2016); https://github.com/kelvinabrokwa/bottleneck




## Installing from source

Requirements: You must have CMake>=2.8.12 and a C++11 compatible compiler (GCC>=4.8) to build.

```
  git clone https://github.com/infophysics/TDA.git
  cd TDA
  sudo python3 setup.py install
```
## Implementation
```python
    grid = [[1,1,1,1,1],[1,0,0,0,1],[1,0,0,0,1],[1,0,0,0,1],[1,1,1,1,1]]
    with open("square.csv", 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows(grid)
    #   try 2D von neumann filter
    filt = Filter2D()
    filt.loadBinaryFromFile("square.csv")
    filt.filterBinaryVonNeumann(10)
    filt.saveBinaryFiltration("square2.csv")
    cube2D = CubicalRipser2D()
    convert_csv_to_dipha("square2.csv", "square_dipha.csv")
    cube2D.ComputeBarcode("square_dipha.csv", "test.csv", "DIPHA", "LINKFIND", 10, True)
    barcode = cube2D.getBarcode()
    plot_persistence_diagram(barcode)

    #   Now try with RIPSER directly
    rips = Ripser()
    # save as point cloud format
    save_binary_cells_to_point_cloud(grid, "square_cloud.csv")

    # Run ripser on this
    rips.ComputeBarcode("square_cloud.csv", 2, 10, 1, "point-cloud", 1)

    # Plot the barcode
    barcode2 = rips.getBarcode()
    plot_persistence_diagram(barcode2)
```

For more examples on possible calls, please see the Jupyter Notebook example.



### Support

* Bugs: Please report bugs to the [issue tracker on Github](https://github.com/infophysics/TDA/issues) such that we can keep track of them and eventually fix them.  Please explain how to reproduce the issue (including code) and which system you are running on.
* Help: Help can be provided also via the issue tracker by tagging your issue with 'question'
* Contributing:  Please fork this repository then make a pull request.  In this pull request, explain the details of your change and include tests.

## Technical implementation

This package is a [pybind11](https://pybind11.readthedocs.io/en/stable/intro.html) wrapper of several persistent homology algorithm implementations as well as general purpose plotting and a computation of the [Persistent Homology Dimension](https://people.math.osu.edu/schweinhart.2/MeasuringShapeWithTopology.pdf)

* Help from Chris Tunnel, which included a great [binding tutorial](https://indico.cern.ch/event/694818/contributions/2985778/attachments/1682465/2703470/PyHEPTalk.pdf)
* Implementation also based on [this](http://www.benjack.io/2018/02/02/python-cpp-revisited.html)

See AUTHORS.md for information on the developers.

## Citation

When you use `TDA`, please say so in your slides or publications (for publications, see Zenodo link above).  You can mention this in addition to how you cite TDA.  This is important for us being able to get funding to support this project.
