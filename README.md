Make Augustus hints file from GFF3
==================================

The gene prediction tool [Augustus](http://bioinf.uni-greifswald.de/augustus/) can use hints to supplement its ab-initio gene prediction
algorithm. At SANBI, to generate these hints we use [exonerate](https://www.ebi.ac.uk/~guy/exonerate/) with either transcript
or protein sequence as input and wrap this in a tool that generates GFF3.

The hints format for Augustus is largely similar to GFF3, and the package provides a tool (exonerate2hints.pl) to parse exonerate output
and produce this hints format. This script emulates the behaviour of the tool, including the default trimming of 15 base pairs off 
the start and end of each CDS. (The usage string for exonerate2hints.pl mentions 9 bases, but the code uses 15 by default.) In exonerate
output, cds and exon lines appear to be interchangeable, and thus this tool optionally converts exon features to CDS.