# <p name="h1">jasmine genome</p>

Some scripts for plot in jasmine genome

## <a name="C1">Installation </a>

### <a name="C1">pyMSAviz must be installed </a>
```sh
pip uninstall pyMSAviz
git clone https://github.com/HansongYan666/pyMSAviz.git
cd pyMSAviz && pip install ./
```

## <a name="C2">Requirement</a>

- [**Biopython**](https://github.com/biopython/biopython) Needs to be installed with pip.
- [**pyfaidx**](https://github.com/mdshw5/pyfaidx.git) Needs to be installed with pip.
- [**muscle**](https://github.com/rcedgar/muscle.git) Needs to be installed in PATH(for script plotmsa.py).
- **[vmatch](https://github.com/uwb-linux/vmatch)** Needs to be installed in PATH(for script RepFind_stat.py).

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

### find long repeat

This script is used for long repeat identification and statistical analysis of the results in chloroplast and mitochondrial genomes. Before using the script, please ensure that mkvtree and vmatch are installed and have been run to generate preliminary results.

```shell
python PATH_to_jasmine_genome/find_long_repeat/RepFind_stat.py --help
usage: RepFind_stat.py [-h] [-l L] [-seedlength SEEDLENGTH] fasta_list prefix

Script for get chloroplast statistics

positional arguments:
  fasta_list            input fasta file list,separate with ,
  prefix                output prefix

options:
  -h, --help            show this help message and exit
  -l L                  search minimal length. [default: 30]
  -seedlength SEEDLENGTH
                        Specify the seed length. [default: 7]
```
