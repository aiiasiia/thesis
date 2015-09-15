# EMS-GT: Exact Motif Search

A solution to the (l,d) planted motif problem for any arbitrary instance, l <= 17.

- `/src/EMS-GT.java` is plain 32-bit EMS-GT.
- `/src/EMS-GT_32.java` is 32-bit EMS-GT with the block bit-masking optimization.
- `/src/EMS-GT_64.java` is 64-bit EMS-GT with the block bit-masking optimization.

EMS-GT was benchmarked against state-of-the-art algorithms [PMS8](http://engr.uconn.edu/~man09004/PMS8/) and [qPMS9](https://github.com/mariusmni/qpms9/tree/master/qpms9).

## Prerequisites: 
- unix or Cygwin environment
- Java 7+ for EMS-GT
- openmpi library (mpi, mpi_cxx, mpiCC) for PMS8 and qPMS9
- gcc, g++ and make to compile PMS8 and qPMS9 from source

## To benchmark:
Run `./experiments.sh`. Synthetic EMS-GT datasets will be generated in `/datasets`, with FASTA versions in `/datasets/FASTA` for PMS8 and qPMS9. Summaries of experimental runs will be found in `/results`with the naming convention `[ProgramName]-[l,d]`.