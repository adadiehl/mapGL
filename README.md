# mapGL
#### Cross-species mapping of DNA elements combined with gain/loss prediction for non-mapping features based on a phylogenetic maximum parsimony algorithm.

Label input regions as orthologous, gained in the query species, or lost in
the target species. Chain alignment files are used to map features from query
to target and one or more outgroup species. Features that map directly from
query to target are labeled as orthologs, and ortholgous coordinates in the
target species are given in the output. Non-mapping features are assigned as
gains or losses based on a maximum-parsimony algorithm predicting
presence/absence in the most-recent common ancestor. Features that, when
mapped, span multiple chains or multiple chromosomes are silently filtered out
unless the -k argument is given.

Based on bnMapper.py, by Ogert Denas (James Taylor lab):
  * https://github.com/bxlab/bx-python/blob/master/scripts/bnMapper.py
  * https://travis-ci.org/bxlab/bx-python


## Usage

```mapGL.py [-h] [-o FILE] [-t FLOAT] [-g GAP] [-v {info,debug,silent}] [-k] input tree qname tname alignments [alignments ...] ```

## Required Arguments

  | input | Input regions to process. Should be in standard bed format. Only the first four bed fields will be used. |
  | tree | Phylogenetic tree describing relationships of query and target species to outgroups. Must be in standard Newick format. Branch lengths are optional, and will be ignored. |
  | qname | Name of the query species. Regions from this species will be mapped to target species coordinates. |
  | tname | Name of the target species. Regions from the query species will be mapped to coordinates from this species. |
  | alignments | Alignment files (.chain or .pkl): One for the target species and one per outgroup species. Files should be named according to the convention: qname.tname[...].chain.gz, where qname is the query species name and tname is the name of the target/outgroup species. Names used for qname and tname must match names used in the newick tree. |

## Options

  * -h, --help            show this help message and exit
  * -o FILE, --output FILE
                          Output file. Default stdout. (default: stdout)
  * -t FLOAT, --threshold FLOAT
                          Mapping threshold i.e., |elem| * threshold <=
                          |mapped_elem| (default: 0.0)
  * -g GAP, --gap GAP     Ignore elements with an insertion/deletion of this or
                          bigger size. (default: -1)
  * -v {info,debug,silent}, --verbose {info,debug,silent}
                          Verbosity level (default: info)
  * -k, --keep_split      If elements span multiple chains, report the segment
                          with the longest overlap instead of silently dropping
                          them. (This is the default behavior for liftOver.)
                          (default: False)

Copyright 2018, Adam Diehl (adadiehl@umich.edu), Boyle Lab, University of Michigan
