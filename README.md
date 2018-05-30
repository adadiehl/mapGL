# mapGL
#### Cross-species mapping of DNA elements combined with gain/loss prediction for non-mapping features based on a phylogenetic maximum parsimony algorithm.

Label input regions as orthologous, gained in the query species, or lost in
the target species, based on inferred presence/absence in the most-recent
common ancestor (MRCA). Chained alignment files are used to map features from
query to target and one or more outgroup species. Features that map directly from
query to target are labeled as orthologs, and ortholgous coordinates in the
target species are given in the output. Non-mapping features are assigned as
gains or losses based on a maximum-parsimony algorithm predicting presence
or absence in the MRCA. Features that, when mapped, span multiple chains or
multiple chromosomes are silently filtered out unless the -k argument is given.

Based on bnMapper.py, by Ogert Denas (James Taylor lab):
  * https://github.com/bxlab/bx-python/blob/master/scripts/bnMapper.py
  * https://travis-ci.org/bxlab/bx-python


## Usage

```mapGL.py [-h] [-o FILE] [-t FLOAT] [-g GAP] [-v {info,debug,silent}] [-k] input tree qname tname alignments [alignments ...] ```

## Required Arguments

  | Argument | Description |
  |---|---|
  | input | Input regions to process. Should be in standard bed format. Only the first four bed fields will be used. |
  | tree | Phylogenetic tree describing relationships of query and target species to outgroups. Must be in standard Newick format. Branch lengths are optional, and will be ignored. |
  | qname | Name of the query species. Regions from this species will be mapped to target species coordinates. |
  | tname | Name of the target species. Regions from the query species will be mapped to coordinates from this species. |
  | alignments | Alignment files (.chain or .pkl): One for the target species and one per outgroup species. Files should be named according to the convention: qname.tname[...].chain.gz, where qname is the query species name and tname is the name of the target/outgroup species. Names used for qname and tname must match names used in the phylogenetic tree. |

## Options

  | Option | Description |
  |---|---|
  | -h, --help | Show help message and exit. |
  | -o FILE, --output FILE | Output file. (default: stdout) |
  | -t FLOAT, --threshold FLOAT | Mapping threshold i.e., |elem| * threshold <= |mapped_elem| (default: 0.0) |
  | -g GAP, --gap GAP | Ignore elements with an insertion/deletion of this or bigger size. (default: -1) |
  | -v {info,debug,silent}, --verbose {info,debug,silent} | Verbosity level (default: info) |
  | -k, --keep_split | Use liftOver mapping convention for split alignments: eep elements that span multiple chains, reporting the longest aligned segment, instead of silently dropping them. (default: False) |

## Output

Predictions are reported in tab-delimited format with the first four columns following the BED4 convention. The predicted evolutionary history (i.e., ortholog, gain in query, or loss in target) is reported in the "status" column. The final three columns contain the mapped location, in target coordinates, of mapped (ortholog) elements.

| Column | Description |
|---|---|
| chrom | Chromosome on which the query element is located. |
| start | Start position on query chromosome. |
| end | End position on query chromosome. |
| name | Element name or ID. |
| status | Predicted phylogenetic history: __ortholog__, __gain_qname__, or __loss_tname__ |
| mapped chrom | For mapped (ortholog) elements, the chromosome on which the mapped element is located, in target coordinates. |
| mapped start | For mapped (ortholog) elements, the start position on the target chromosome on which the mapped element is located. |
| mapped end | For mapped (ortholog) elements, the end position on the target chromosome on which the mapped element is located. |

Copyright 2018, Adam Diehl (adadiehl@umich.edu), Boyle Lab, University of Michigan
