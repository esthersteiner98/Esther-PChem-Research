# PChem-Research
The calculate deltas program inputs molecule geometry (x,y,z) coordinates from each input file (txt format).  The input files are assumed to contain a list of atoms with coordinates.  The atoms must be in the same order in all input files.

Calculates:
1) bond distances for each input file.  Distances greater than 1.9 A between atoms is assumed to be non-bonding.
2) deltas - for input files 2, 3, 4, ... calculate the difference between the reference bond length and the corresponding bond length of the input file.    Input file 1 is assumed to be the reference.

Here's an example of the input format:

 C          1.22172        0.00000       -0.69374
 C          0.00000       -0.00000        1.42470
 C          2.43721        0.00000        1.42164
 C          2.45857        0.00000        0.01186
 C          0.00000        0.00000        0.00000
 H          1.22897        0.00000       -1.77455
 C          1.23785       -0.00000        2.09459
 H          1.24346       -0.00000        3.17557
 C         -1.22172        0.00000       -0.69374
 C         -2.45857        0.00000        0.01186
 C         -2.43721        0.00000        1.42164
 H         -1.22897        0.00000       -1.77455
