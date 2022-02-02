# CircRNA analysis suite of tools
## Circulator

Circulator is a multi-threaded software that scans for CDSs in circRNA suquences. 
Provide an input FASTA file that contains circRNA sequence records and the software will scan and output predicted protein-formatted fasta records.


## Comandline help
```
usage: Circulator.py [-h] [-o Output] [-m MIN] [-t THREADS] [-l] [-r] [-v]
                     Input

       .__                     .__                       __   
  ____ |__|______   ____  __ __|  | _____ _______  _____/  |_ 
_/ ___\|  \_  __ \_/ ___\|  |  \  | \__  \_  __ \/  _ \   __\
\  \___|  ||  | \/\  \___|  |  /  |__/ __ \|  | \(  <_> )  |  
 \___  >__||__|    \___  >____/|____(____  /__|   \____/|__|  
     \/                \/                \/                   
Translate potential CDSs in circRNAs.
By Amin Mahpour.

positional arguments:
  Input                 Input FASTA file.

optional arguments:
  -h, --help            show this help message and exit
  -o Output, --out Output
                        Output protein FASTA file.
  -m MIN, --min MIN     Cutoff CDS size. default: 20.
  -t THREADS, --threads THREADS
                        Requested processor numbers. Default(0) uses all
                        available cores. example: -t 4
  -l, --longest         Only keep longest ORFs.
  -r, --rna             Input fasta file is RNA sequence.
  -v, --version         show program's version number and exit

Examples:
Circulator.py input.fa > output.fa #Simple usage
Circulator.py -t 4 input.fa -o output.fa #Add 4 threads
Circulator.py -t 4 -l input.fa -o output.fa #Add 4 threads and keep the longest CDS.
Circulator.py -t 4 -l -m 25 input.fa > output.fa #Add 4 threads and keep the longest CDS and set a cutoff of 25 amino acids.

```
