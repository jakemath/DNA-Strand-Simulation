# DNA-Strand-Simulation

A randomized simulation of large DNA/RNA strands, from which statistics are computed and reported regarding the composition of the strand.

The statistics are primarily an analysis of the proportions of nitrogen bases, codons, and amino acids present within the strand. The program assumes a uniform probability for each base in the strand, but this can easily be changed in the strand.cpp file.

To run a small driver program: place all the files in the same directory, and build with
  
    ./build.sh
    
Then execute with an integer argument for the size of the DNA strand:

    ./a.out <size>
