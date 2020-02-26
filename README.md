3D Monte Carlo Tomography 

**Authors:**
 - Xin Zhang x.zhang2@ed.ac.uk
 - Andrew Curtis

## MCTomo

3D Monte Carlo tomography package using reversible jump McMC method and 3D Voronoi tessellation.

Compilation and installation
============================

## Requirements

* A C compiler that supports C++11 and a Fortran compiler which supports Fortran 2003.

* CGAL 4.8 or later ( for 3d delaunay and Voronoi support) CGAL can be downloaded from (http://www.cgal.org/) and can be built and installed following [CGAL installation manual] (http://doc.cgal.org/latest/Manual/installation.html). 

If you are using MAC OS X, CGAL can be installed by following way:
     
     sudo port install cgal
  
  On debian/Ubuntu, 
     
     sudo apt-get install libcgal-dev

* NetCDF4 or later ( support for storaging the samples in a protable way ). If you do not have a netcdf4 library, you can install it through a package management program, such as rpm, yum, homebrew, macports, adept, and others.  Alternatively, you can download it from (http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html).
  This is not necessary. If you do not want to use netcdf4 library, just leave the NETCDF variable undefined in the makefile.

## Algorithms and Dependencies

** Normal Modes **

* The surface wave modal approximation code is from COMPUTER PROGRAM IN SEISMOLOGY wrote by Herrmann, R.B. (http://www.eas.slu.edu/eqc/eqccps.html).

** Fast Marching **

* The code used the fast marching method to calculate the first arrival times. The 2D fast marching code is from Nick Rawlinson (http://rses.anu.edu.au/seismology/soft/fmmcode/). For the 3D body wave tomography, we used a general N-dimensional fast marching code by Gomez, J.V. and Pardeiro, J. (https://github.com/jvgomez/fast_methods). Both of them are revised for our usage.

** 3D Voronoi **

* The 3D Voronoi tessellation is calculated based on CGAL (see above). However, we still used the KD-tree method to convert the Voronoi diagram to a regular grid based representation (https://github.com/jmhodges/kdtree2).

## Compilation

* Add your cgal include path (e.g. MYCGAL/include) and link path in the makefile
* If available and necessary, add your netcdf include path and link path in the makefile
* Then make

## Examples
Examples are in the examples directory

## References

CGAL, Computational Geometry Algorithms Library, http://www.cgal.org.

Chen, X. (1993). A systematic and efficient method of computing normal modes for multilayered half-space. Geophysical Journal International, 115(2), 391-409.

Herrmann, R. B. (2013). Computer programs in seismology: An evolving tool for instruction and research. Seismological Research Letters, 84(6), 1081-1088.

Kennel, M. B. (2004). KDTREE 2: Fortran 95 and C++ software to efficiently search for near neighbors in a multi-dimensional Euclidean space. arXiv preprint physics/0408067.

Rawlinson, N., & Sambridge, M. (2004). Wave front evolution in strongly heterogeneous layered media using the fast marching method. Geophysical Journal International, 156(3), 631-647.

Valero-Gomez, A., Gomez, J. V., Garrido, S., & Moreno, L. (2013). The path to efficiency: fast marching method for safer, more efficient mobile robot trajectories. IEEE Robotics & Automation Magazine, 20(4), 111-120.

Wu, B., & Chen, X. (2016). Stable, accurate and efficient computation of normal modes for horizontal stratified models. Geophysical Journal International, 206(2), 1281-1300.

Zhang, X., Curtis, A., Galetti, E., & de Ridder, S., 2018. 3-D Monte Carlo surface wave tomography, Geophysical Journal International, 215(3), 1644–1658.

Zhang, X., Hansteen, F., Curtis, A., & de Ridder, S., 2019a. 1D, 2D and 3D Monte Carlo ambient noise tomography using a dense passive seismic array installed on the North Sea seabed, Journal of Geophysical
Research: Solid Earth.

Zhang, X., Roy, C., Curtis, A., Nowacki, A., & Baptie, B., 2019b. 3D Tomographic Monte Carlo joint inversion of earthquake body wave travel times and ambient noise surface wave dispersion – imaging the
subsurface using induced seismicity and ambient noise, submitted.




