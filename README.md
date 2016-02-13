# LLG-Code

1.installation:

To compile this code, Eigen 3 and GNU Scientific Library (GSL) are  required to be installed.

For more information:

http://eigen.tuxfamily.org/index.php?title=Main_Page

http://www.gnu.org/software/gsl/

The magnon binary has been build under MAC OSX 10.11. To modify the source code, you're required to rebuild the makefile with CMake. 

2.About the program:

Current version of Magnon can be used to solve the Landau–Lifshitz–Gilbert (LLG) equation in two-dimensional square/ hexagonal magnetic superlattice with circle or ellipse cross-section. Based on the eigenvalues and eigenvectors from LLG, associated Green functions and time response of spin wave pulse can be achieved. With the solved band structure, density of states of magnon and vector plot of total magnetic/electric field can be calculated. For more theoretical derivation, please refer to the following paper:

http://arxiv.org/abs/1111.2506
