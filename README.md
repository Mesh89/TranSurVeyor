# TranSurVeyor

## Compiling

To be compiled, TranSurveyor requires g++ 4.7.2 or higher and htslib 1.4 or higher (https://github.com/samtools/htslib). It was tested using g++ 4.7.2 and htslib 1.6, so those versions are recommended.
Make sure that htslib is in your library search path, and run 
```
cmake . && make
```

## Running 

Once the c++ code is compiled, TranSurveyor can be run. Python and libraries NumPy (http://www.numpy.org/), PyFaidx (https://github.com/mdshw5/pyfaidx) and PySam (https://github.com/pysam-developers/pysam) are required. Python 2.7, NumPy 1.10, PyFaidx 0.4 and PySam 0.12 are the recommended versions.

The bare minimum command for running TranSurveyor is 
```
python surveyor.py /path/to/bamfile /an/empty/working/directory /path/to/reference/fasta
```

Note also that the reference fasta file should be indexed by both bwa and samtools. For example, assuming the file is hg19.fa, you should run
```
bwa index hg19.fa
samtools faidx hg19.fa
```

Other parameters which may be important are the number of threads and the location of the bwa and samtools executable.
```
python surveyor.py /path/to/bamfile /an/empty/working/directory	/path/to/reference/fasta --threads 40 --samtools /path/to/samtools --bwa /path/to/bwa
```

After TranSurveyor has been successfully run, the calls can be retrieved with the command
```
./filter /path/to/working/directory
```
