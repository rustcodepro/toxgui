# toxannotator

- a toxdb annotator for genomic comparison for comparative genomics. 
- ready to use tables for comparative analysis. 
- genomic comparison and annotation plotters. 

```
cargo build
```

```
_____                                                  _             _                  
|_   _|   ___   __  __   __ _   _ __    _ __     ___   | |_    __ _  | |_    ___    _ __ 
 | |    / _ \  \ \/ /  / _` | | '_ \  | '_ \   / _ \  | __|  / _` | | __|  / _ \  | '__|
 | |   | (_) |  >  <  | (_| | | | | | | | | | | (_) | | |_  | (_| | | |_  | (_) | | |   
 |_|    \___/  /_/\_\  \__,_| |_| |_| |_| |_|  \___/   \__|  \__,_|  \__|  \___/  |_|   
                                                                                        

A toxodb annotator.
     ************************************************
     Gaurav Sablok,
     Email: codeprog@icloud.com
    ************************************************

Usage: toxannotator <COMMAND>

Commands:
protein-compare              Only compare protein coding
protein-compare-seq-command  Compare protein coding coordinates and sequences also
protein-plotter              Plot the protein coding regions
protein-tensor               Prepare the protein tensor for the machine and deep learning
extract-seq                  plot the specific ids information
comparem-rna                 compare mRNA annotations
compare-fasta-gff            compare fasta and gff for the same and different ids and sequences
help                         Print this message or the help of the given subcommand(s)

Options:
-h, --help     Print help
-V, --version  Print version


```

```
Only compare protein coding

Usage: toxannotator protein-compare <GFFFILE1> <GFFFILE2>

Arguments:
  <GFFFILE1>  path to the first gff file
  <GFFFILE2>  path to the second gff file

Options:
  -h, --help  Print help
  
toxannotator protein-compare ./testfiles/a1.gff ./testfiles/b1.gff  4

commonggenes-same-strand.txt
TGME49_200010	TGME49_200010	2245476	2245476	2249210	2248187	-	-

GeneName	Start1	Start2	End1	End2	Strand1	Strand2	Start Difference	End Difference
TGME49_200010	2245476	2245476	2249210	2248187	-	-	0	1023

```

```

Compare protein coding coordinates and sequences also

Usage: toxannotator protein-compare-seq-command <GFFFILE1> <GFFFILE2> <FASTAFILE_1> <FASTAFILE_2> <THREADS>

Arguments:
  <GFFFILE1>     path to the first gff file
  <GFFFILE2>     path to the second gff file
  <FASTAFILE_1>  fasta1 file
  <FASTAFILE_2>  fasta2 file
  <THREADS>      threads

Options:
  -h, --help  Print help
```

```
Plot the protein coding regions

Usage: toxannotator protein-plotter <INPUTFILE1> <INPUTFILE2> <THREADS>

Arguments:
  <INPUTFILE1>  input file 1
  <INPUTFILE2>  input file 2
  <THREADS>     threads

Options:
  -h, --help  Print help
  
toxannotator protein-plotter ./testfiles/a1.gff ./testfiles/b1.gff 4
```

```
Prepare the protein tensor for the machine and deep learning

Usage: toxannotator protein-tensor <INPUTFILE> <THREADS>

Arguments:
  <INPUTFILE>  input file for the tensor
  <THREADS>    threads

Options:
  -h, --help  Print help

```

```
plot the specific ids information

Usage: toxannotator extract-seq <ANNOTATIONFILE> <IDSFILE> <THREADS>

Arguments:
  <ANNOTATIONFILE>  file to the annotation
  <IDSFILE>         idsfile
  <THREADS>         threads

Options:
  -h, --help  Print help
  
toxannotator extract-seq ./testfiles/a1.gff ./testfiles/id.test 4
```

```
compare mRNA annotations

Usage: toxannotator comparem-rna <GFF_1> <GFF_2> <THREADS>

Arguments:
  <GFF_1>    path to the first gff
  <GFF_2>    path to the second gff
  <THREADS>  threads

Options:
  -h, --help  Print help
```
```
compare fasta and gff for the same and different ids and sequences

Usage: toxannotator compare-fasta-gff <GFF_1_INPUT> <GFF_2_INPUT> <FASTA_1_INPUT> <FASTA_2_INPUT> <THREADS>

Arguments:
  <GFF_1_INPUT>    path to the first gff
  <GFF_2_INPUT>    path to the second gff
  <FASTA_1_INPUT>  path to the first fasta
  <FASTA_2_INPUT>  path to the second fasta
  <THREADS>        threads

Options:
  -h, --help  Print help
```

- To install windows version:

```
rustup component add llvm-tools
rustup target add x86_64-pc-windows-msvc
git clone https://github.com/IBCHgenomic/ensemblcov.git
cd ensemblcov
cargo xwin build --target x86_64-pc-windows-msvc

Gaurav Sablok \
codeprog@icloud.com
