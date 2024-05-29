# <p name="h1">jasmine genome</p> 
Some scripts for plot in jasmine genome

## <a name="C1">Installation </a>
```sh
git clone https://github.com/HansongYan666/jasmine_genome.git
```


## <a name="C2">Requirement</a>
- [**Biopython**](https://github.com/biopython/biopython) Needs to be installed with pip.
- [**pyfaidx**](https://github.com/mdshw5/pyfaidx.git) Needs to be installed with pip.
- [**pyMSAviz**](https://github.com/moshi4/pyMSAviz/) Needs to be installed with pip.
- [**muscle**](https://github.com/rcedgar/muscle.git) Needs to be installed in PATH(for script plotmsa.py).


## <a name="C3">Options and usage</a>

### plot base and AA

```shell
python plotbase.py -h 
usage: plotbase.py [-h] [-i INFILE] [-o OUTFILE]

plot base and aa sequences

options:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        the input fasta file
  -o OUTFILE, --output OUTFILE
                        the output png file
```


### plot mas alignment base and AA

```shell
python plotmas.py -h  
usage: plotmas.py [-h] [-i INFILE] [-o OUTFILE]

plot msa align and aa sequences, muscle should be in PATH

options:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        the input fasta file
  -o OUTFILE, --output OUTFILE
                        the output png file
```


