//Parameters for the coalescence simulation program : fsimcoal2.exe
1 samples to simulate :
//Population effective sizes (number of genes)
NCurrent$
//Samples sizes and samples age 
100
//Growth rates	: negative growth implies population expansion
0
//Number of migration matrices : 0 implies no migration between demes
0
///hist event: time, source, sink, migrants, bot intensity, growth rate, migr. mat
1 historical event
Time$ 0 0 0 botIntensity$ 0 0 instbot
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   1.0e-8 OUTEXP
