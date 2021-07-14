## Marching algorithms
These are methods used to visualize volumetric data, by finding iso-contours of the underlying function that the data is supposedly generated from. The 2D and 3D analogues of the algorithm are the marching cubes and squares, respectively.  

### Useful resources
- The [original paper](https://dl.acm.org/doi/10.1145/37401.37422) by Lorenson and Cline.
- [Isosurfaces: Geometry, Topology, and Algorithms](https://books.google.co.in/books/about/Isosurfaces.html?id=62DOBQAAQBAJ&redir_esc=y)

This was done in Winter 2020, as an attempt to learn basic fortran. The fortran code calculates vertices of the triangles, and produces a .vtk file which is then visualized in [ParaView](https://www.paraview.org/). There is a single python file to manage vtk formatting, as I could not wrangle with fortran to get rid of some whitespace.
