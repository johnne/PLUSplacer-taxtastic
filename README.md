

------------------------
Summary
------------------------

This repository extends the [pplacer-SCAMPP repository](https://github.com/chry04/PLUSplacer), which llows pplacer and EPA-ng to run on ultra-large reference trees. This repository contains code for pplacer-SCAMPP-FastTree. 

The main conclusion is that using FastTree instead of RAxML to estimate numeric parameters on the backbone tree allows pplacer to place query sequences into much larger backbone trees, with similar scalability and better accuracy than comparably scalable
methods. For further details, please refer to the [paper](https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad008/7009227).

They are both python programs that can be run on **Linux and MacOS**

------------------------
Requirements
------------------------
1. Python: Version >= 3.0
2. Treeswift
3. numpy
4. pplacer: Version 1.1.alpha19 (for pplacer-SCAMPP)
   or 
   EPA-ng: Version 0.3.8 (for EPA-ng-SCAMPP)
5. taxtastic (can be found [here](https://github.com/fhcrc/taxtastic)). Please install from the git repository with the following steps, inside this PLUSplacer-taxtastic repository.

```
git clone https://github.com/fhcrc/taxtastic.git
cd taxtastic
python3 -m venv taxtastic-env
source taxtastic-env/bin/activate
pip install .
```

----------------------------------
Input & Output Specification
----------------------------------

The input parameters for pplacer-tax-SCAMPP are as follows:
    
    Required arguments: 
    -i, --info : statistics file path produced by FastTree containing the substitution rates
    -t, --tree : backbone tree file path (tree expected in newick format) (parameter T)
    -d, --outdir : directory to where output file pplacer-SCAMPP.jplace is written
    -a, --alignment : fasta format file path containing the multiple sequence alignment (MSA) of reference and query sequences 
    -r, --refaln : fasta format file path containing the multiple sequence alignment (MSA) of reference sequences

    Optional arguments:
    -m, --model : DNA substitution model such as GTR for nucleotides
    -q, --qalignment : file containing a list of query sequences aligned to the reference MSA in fasta format (needed when not contained in file with reference MSA)
    -b, --subtreesize : maximum size of the subtree for placement (with pplacer 2000 is recommended, and with EPA-ng 10,000 is recommended) (parameter B) 
    -s, --subtreetype : options for collecting nodes for the subtree - d (default) for edge weighted distances, n for node distances, h for Hamming distances
    -n, --tmpfilenbr : number for working directories in the current file path (to allow to run multiple instances concurently)
    -q, --qalignment : path to query sequence alignment in fasta format (if reference alignment is not provided)
    -f, --fragmentflag : boolean, True if queries contain fragmentary sequences, and you wish to mask leading and trailing gaps in each query sequence when finding closest sister taxon.
    -v, --version : show the version number
  
----------------------------------
Usage
----------------------------------
Please run from a directory containing the respective phylogenetic placement method. This would be pplacer (available at https://github.com/matsen/pplacer) for pplacer-SCAMPP or EPA-ng (available at https://github.com/Pbdas/epa-ng) for EPA-ng-SCAMPP.

python3 pplacer-tax-SCAMPP.py -i INFO -t TREE -d OUTDIR -a ALIGNMENT
