#!/usr/bin/env python

"""
Label input regions as orthologous, gained in the query species, or lost
in the target species. Chain alignment files are used to map features from query
to target and one or more outgroup species. Features that map directly from query
to target are labeled as orthologs, and ortholgous coordinates in the target 
species are given in the output. Non-mapping features are assigned as gains or
losses based on a maximum-parsimony algorithm predicting presence/absence in the 
most-recent common ancestor.

Based on bnMapper.py, by Ogert Denas (James Taylor lab)
(https://github.com/bxlab/bx-python/blob/master/scripts/bnMapper.py)
"""
import logging
import os
import sys
from itertools import groupby
from operator import attrgetter, concat, itemgetter
import numpy as np
from six.moves import reduce

import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)

from lib.phylo import newick
from lib.bx.align import epo
from lib.bx.align.epo import bed_union as elem_u
from lib.bx.cookbook import argparse
from lib.bx.intervals.intersection import IntervalTree, Interval

elem_t = np.dtype([('chrom', np.str_, 30), ('start', np.int64), ('end', np.int64), ('id', np.str_, 100)])
narrowPeak_t = np.dtype([('chrom', np.str_, 30), ('start', np.int64), ('end', np.int64), ('id', np.str_, 100),
                         ('score', np.int64), ('strand', np.str_, 1), ('signalValue', np.float),
                         ('pValue', np.float), ('qValue', np.float), ('peak', np.int64)])
LOG_LEVELS = {"info" : logging.INFO, "debug" : logging.DEBUG, "silent" : logging.ERROR}

logging.basicConfig()
log = logging.getLogger()

class GIntervalTree( IntervalTree ):
    """a set of IntervalTrees that is indexed by chromosomes"""

    def __init__(self, data=[]):
        self._trees = {}

    def add(self, chrom, element):
        """insert an element. use this method as the IntervalTree one.
        this will simply call the IntervalTree.add method on the right tree

        :param chrom: chromosome
        :param element: the argument of IntervalTree.insert_interval
        :return: None
        """

        self._trees.setdefault(chrom, IntervalTree()).insert_interval( element )

    def find(self, chrom, start, end):
        """find the intersecting elements

        :param chrom: chromosome
        :param start: start
        :param end: end
        :return: a list of intersecting elements"""

        tree = self._trees.get( chrom, None )
        if tree:
            return tree.find( start, end )
        #return always a list
        return []

def transform(elem, chain_CT_CQ, max_gap):
    """transform the coordinates of this elem into the other species.

    elem intersects this chain's ginterval.
    :return: a list of the type [(to_chr, start, end, elem[id]) ... ]"""
    (chain, CT, CQ) = chain_CT_CQ
    start, end = max(elem['start'], chain.tStart) - chain.tStart, min(elem['end'], chain.tEnd) - chain.tStart

    assert np.all( (CT[:,1] - CT[:,0]) == (CQ[:,1] - CQ[:,0]) )
    to_chrom = chain.qName
    to_gab_start = chain.qStart

    start_idx = np.where( CT[:,1] > start )[0][0]
    end_idx = np.where( CT[:,0] < end )[0][-1]

    if start_idx > end_idx: #maps to a gap region on the other species
        return []

    ## apply the gap threshold
    if max_gap >= 0 and start_idx < end_idx - 1:
        if np.max(CT[(start_idx+1):end_idx,0] - CT[start_idx:(end_idx-1),1]) > max_gap or np.max(CQ[(start_idx+1):end_idx,0] - CQ[start_idx:(end_idx-1),1]) > max_gap:
            return []

    assert start < CT[start_idx, 1]
    assert  CT[end_idx, 0] < end
    to_start = CQ[start_idx, 0] + max(0, start - CT[start_idx,0]) # correct if on middle of interval
    to_end = CQ[end_idx, 1] - max(0, CT[end_idx, 1] - end)        # idem

    if start_idx == end_idx: #elem falls in a single run of matches
        slices = [(to_start, to_end)]
    else:
        slices = [(to_start, CQ[start_idx,1])]
        slices += [(CQ[i,0], CQ[i,1]) for i in range(start_idx+1, end_idx)]
        slices.append( (CQ[end_idx,0], to_end) )
    if chain.qStrand == '-':
        Sz = chain.qEnd - chain.qStart
        slices =  [(Sz-t[1], Sz-t[0]) for t in slices]
    return [(to_chrom, to_gab_start + t[0], to_gab_start + t[1], elem['id']) for t in slices]

def union_elements(elements):
    """elements = [(chr, s, e, id), ...], this is to join elements that have a
    deletion in the 'to' species
    """

    if len(elements) < 2: return elements
    assert set( [e[3] for e in elements] ) == set( [elements[0][3]] ), "more than one id"
    el_id = elements[0][3]

    unioned_elements = []
    for ch, chgrp in groupby(elements, key=itemgetter(0)):
        for (s, e) in elem_u( np.array([itemgetter(1, 2)(_) for _ in chgrp], dtype=np.uint) ):
            if (s < e):
                unioned_elements.append( (ch, s, e, el_id) )
    assert len(unioned_elements) <= len(elements)
    return unioned_elements

def transform_by_chrom(all_epo, from_elem, tree, chrom, opt):
    """ 
    AGD: Dropped loop over multiple elements used in bnMapper,
    so this now only maps a single element.
    """
    elems_mapped = []

    mapped_elem_count = 0
    matching_block_ids = [attrgetter("value")(_) for _ in tree.find(chrom, from_elem['start'], from_elem['end'])]

    # do the actual mapping
    to_elem_slices = [_ for _ in (transform(from_elem, all_epo[i], opt.gap) for i in matching_block_ids) if _]
    max_elem_idx = 0
    if len(to_elem_slices) == 0:
        log.debug("%s: no match in target: discarding." % (str(from_elem)))
        return elems_mapped
    elif len(to_elem_slices) > 1:
        log.debug("%s spans multiple chains/chromosomes. Using longest alignment." % (str(from_elem)))
        max_elem_len = 0
        for i in xrange(len(to_elem_slices)):
            elem_len = to_elem_slices[i][-1][2] - to_elem_slices[i][0][2]
            if elem_len > max_elem_len:
                max_elem_len = elem_len
                max_elem_idx = i
    elif len(to_elem_slices) > 1 and opt.drop_split:
        log.debug("%s spans multiple chains/chromosomes: discarding." % (str(from_elem)))
    to_elem_slices = to_elem_slices[max_elem_idx]

    # apply threshold
    if (from_elem[2] - from_elem[1]) * opt.threshold > reduce(lambda b,a: a[2]-a[1] + b, to_elem_slices, 0):
        log.debug("%s did not pass threshold" % (str(from_elem)))

    # if to_species had insertions you can join elements
    to_elem_list = sorted(union_elements(to_elem_slices), key=lambda a: a[1])
    if to_elem_list:
        mapped_elem_count += 1
        log.debug("\tjoined to %d elements" % (len(to_elem_list)))
        # Don't store the whole to_elem_list -- just the chrom,
        # start, and end of the mapped region as a whole.
        start = to_elem_list[0][1]
        end = to_elem_list[-1][2]
        peak = int((start + end)/2) - start
        if opt.in_format == "narrowPeak":
            # Map the peak location
            matching_block_ids = [attrgetter("value")(_) for _ in tree.find(chrom, from_elem['peak'], from_elem['peak'])]
            p_elem_slices = [_ for _ in (transform( np.array((chrom, from_elem['peak'], from_elem['peak'], '.'), dtype=elem_t), all_epo[i], opt.gap) for i in matching_block_ids) if _]
            if len(p_elem_slices) >= 1:
                #sys.stderr.write("{}\n".format(p_elem_slices))
                # Make sure the peak is between the start and end positions
                if p_elem_slices[0][0][1] >= start and p_elem_slices[0][0][1] <= end:
                    peak = p_elem_slices[0][0][1] - start
                    
        elems_mapped.append([to_elem_list[0][0], start, end, peak, from_elem['id']])
    return elems_mapped

def transform_file(ELEMS, ofname, TREES, leaves, phylo_full, phylo_pruned, opt):
    """
    Transform each element of the input file to the target. If the
    element maps, print the labelled element, with mapped coordinates,
    to the output file on 'ofname'. If it does not, map to each of the
    outgroup species and store a label to indicate positive/negative
    mapping at each leaf, then perform a post-order tree traversal to
    predict presence/absence of the element at the tree root (most-recent
    common ancestor). According to maximum parsimony, elements present in
    the MRCA are most likely lost in the target species, while absence in
    the MRCA most likely inidcates gain in the query species. Element
    labels: ortholog, gain-(qname), loss-(tname), ambiguous. In cases
    where an unambiguous gain/loss call cannot be made (which should be
    rare!), the most-terminal outgroup will be pruned from the tree,
    which normally results in an unambigous call. If this still does not
    work, the "ambiguous" label will be applied.
    """
    OUT_FRM = "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\n"


    log.info("Mapping {} features...".format(len(ELEMS.flat)))
    with open(ofname, 'w') as out_fd:

        for elem in ELEMS.flat:
            peak = int((elem[2] - elem[1])/2)
            if opt.in_format == "narrowPeak":
                peak =  elem[9] - elem[1]
            
            # First try mapping to the target species.
            mapped_el = map_elem(elem, TREES[opt.tname]["EPO"],
                                 TREES[opt.tname]["TREE"], opt)

            if mapped_el:
                #sys.stderr.write("{}\n".format(mapped_el))
                out_fd.write(OUT_FRM % (elem[0], elem[1], elem[2], elem[3], peak,
                                        "ortholog", mapped_el[0], mapped_el[1], mapped_el[2], mapped_el[3]))
            
            else:
                # No match to target. Map to outgroups.
                mapped_spp = {opt.qname: 1,
                              opt.tname: 0} # Mapping labels for query and target pre-filled
                for spp in leaves:
                    if spp == opt.tname or spp == opt.qname:
                        continue
                    mapped_el = map_elem(elem, TREES[spp]["EPO"],
                                         TREES[spp]["TREE"], opt)
                    if mapped_el:
                        mapped_spp[spp] = 1
                    else:
                        mapped_spp[spp] = 0
                log.debug("Species mappings for elem {}: {}\n".format(elem, mapped_spp))

                # Now do a post-order traversal of the tree, finding
                # the most-likely lable at each ancestral node, down
                # to the root node.
                in_mrca = infer_ancestral(phylo_full, mapped_spp)

                orth_str = "ambiguous"
                if in_mrca == 0:
                    orth_str = "gain_{}".format(opt.qname)
                elif in_mrca == 1:
                    orth_str = "loss_{}".format(opt.tname)
                elif in_mrca == "N":
                    # If the full tree gives an ambiguous result, try the
                    # pruned tree.
                    in_mrca = infer_ancestral(phylo_pruned, mapped_spp)
                    if in_mrca == 0:
                        orth_str = "gain_{}".format(opt.qname)
                    elif in_mrca == 1:
                        orth_str = "loss_{}".format(opt.tname)

                out_fd.write(OUT_FRM % (elem[0], elem[1], elem[2], elem[3], peak,
                                        orth_str, ".", -1, -1, -1))
    log.info("DONE!")


def map_elem(elem, EPO, TREE, opt):
    """
    map a single elem
    """
    mapped_els = transform_by_chrom(EPO, elem, TREE,
                                    elem['chrom'], opt)
    if mapped_els:
        if len(mapped_els) > 1:
            log.info("Warning: %s maps to multiple locations. Returning only the first mapping.\n" % elem_qry)
        return(mapped_els[0])
    return


def infer_ancestral(phylo, mapped_spp):
    """
    Infer the ancestral state of the most-recent common ancestor
    using maximum parsimony. Returns 1 for present in MRCA, 0 for
    absent in MRCA, -1 for ambiguous.
    """
    # Set mapping labels in leaves
    phylo.visit(lambda n: setattr(n, 'label', mapped_spp[n.name]),
                lambda n: n.is_leaf,
                mode = "postorder")

    # Infer ancestral labels
    phylo.visit(lambda n: n.infer_labels(),
                lambda n: not n.is_leaf,
                mode = "postorder")

    # Return the label at the root    
    return(phylo.get_label_at_root())
    
    
def loadChains(path):
    "name says it."

    EPO = epo.Chain._parse_file(path, True)
    ## convert coordinates w.r.t the forward strand (into slices)
    ## compute cummulative intervals
    for i in range( len(EPO) ):
        ch, S, T, Q = EPO[i]
        if ch.tStrand == '-':
            ch = ch._replace(tEnd = ch.tSize - ch.tStart,
                    tStart = ch.tSize - ch.tEnd)
        if ch.qStrand == '-':
            ch = ch._replace(qEnd = ch.qSize - ch.qStart,
                    qStart = ch.qSize - ch.qEnd)
        EPO[i] = (ch,
                epo.cummulative_intervals(S, T),
                epo.cummulative_intervals(S, Q)
                )
    ##now each element of epo is (chain_header, target_intervals, query_intervals)
    assert all( t[0].tStrand == '+' for t in EPO ), "all target strands should be +"
    return EPO

def loadFeatures(path):
    "load bed4 features (all other columns are ignored)"

    log.info("loading from %s ..." % path)
    data = []
    with open(path) as fd:
        for line in fd:
            cols = line.split()
            data.append( (cols[0], int(cols[1]), int(cols[2]), cols[3]) )
    return np.array(data, dtype=elem_t)

def loadFeatures(path, opt):
    """
    Load features. For BED, only BED4 columns are loaded.
    For narrowPeak, all columns are loaded.
    """

    log.info("loading from %s ..." % path)
    data = []
    if opt.in_format == "BED":
        with open(path) as fd:
            for line in fd:
                cols = line.split()
                data.append( (cols[0], int(cols[1]), int(cols[2]), cols[3]) )
        data = np.array(data, dtype=elem_t)
    else:
        with open(path) as fd:
            for line in fd:
                cols = line.split()
                data.append( (cols[0], int(cols[1]), int(cols[2]), cols[3], int(cols[4]),
                              cols[5], float(cols[6]), float(cols[7]), float(cols[8]),
                              int(cols[-1])+int(cols[1])) )
        data = np.array(data, dtype=narrowPeak_t)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, epilog="Adam Diehl (Boyle Lab)",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input",
            help="Input regions to process. Should be in standard bed format. Only the first four bed fields will be used.")
    parser.add_argument("tree",
            help="Tree, in standard Newick format, with or without branch lengths, describing relationships of query and target species to outgroups.")
    parser.add_argument("qname",
            help="Name of the query species. Regions from this species will be mapped to target species coordinates.")
    parser.add_argument("tname",
            help="Name of the target species. Regions from the query species will be mapped to coordinates from this species.")
    parser.add_argument("alignments", nargs='+',
            help="Alignment files (.chain or .pkl): One for the target species and one per outgroup species. Files should be named according to the convention: qname.tname[...].chain.gz, where qname is the query species name and tname is the name of the target/outgroup species. Names used for qname and tname must match names used in the newick tree.")    
    parser.add_argument("-o", '--output', metavar="FILE", default='stdout',
            type=lambda s: ((s in ('stdout', '-') and "/dev/stdout") or s),
            help="Output file. Default stdout.")
    parser.add_argument("-t", '--threshold', metavar="FLOAT", default=0., type=float,
            help="Mapping threshold i.e., |elem| * threshold <= |mapped_elem|")
    parser.add_argument('-g', '--gap', type=int, default=-1,
            help="Ignore elements with an insertion/deletion of this or bigger size.")
    parser.add_argument('-v', '--verbose', type=str, choices=list(LOG_LEVELS.keys()), default='info',
            help='Verbosity level')
    parser.add_argument("-d", '--drop_split', default=False, action='store_true',
            help="If elements span multiple chains, report them as non-mapping. These will then be reported as gains or losses, according to the maximum-parsimony predictions. This is the default mapping behavior for bnMapper. By default, mapGL.pys will follow the mapping convention used by liftOver, whereas the longest mapped alignment is reported for split elements.")
    parser.add_argument("-i", "--in_format", choices=["BED", "narrowPeak"], default="BED",
            help="Input file format.")

    opt = parser.parse_args()
    log.setLevel(LOG_LEVELS[opt.verbose])

    # Load up the newick tree
    log.info("Parsing species tree: {}".format(opt.tree))
    phylo_full = newick.parse_node(opt.tree)
    log.debug("Full tree:\n{}".format(phylo_full.ascii_art(show_internal=False, strict=True)))

    # Prune the terminal outgroup (furthest from the query species)
    # to use in ambiguous cases. For now, this is assumed to be the
    # last species in the leaves list.
    phylo_pruned = newick.parse_node(opt.tree)
    leaves = phylo_pruned.get_leaves()
    leaves[-1].ancestor.descendants.remove(leaves[-1])
    phylo_pruned = newick.parse_node(leaves[-1].ancestor.newick)
    phylo_pruned.remove_redundant_nodes()
    log.debug("Pruned tree:\n{}".format(phylo_pruned.ascii_art(show_internal=False, strict=True)))
    leaves = phylo_full.get_leaf_names()
    
    if opt.qname not in leaves:
        sys.stderr.write("Query species name {} not present in tree: {}\n".format(opt.qname, phylo_full.newick))
        exit(1)
    if opt.tname not in leaves:
        sys.stderr.write("Target species name {} not present in tree: {}\n".format(opt.tname, phylo_full.newick))
        exit(1)

    # Load up alignment chains. Need reciprocal-best chains for the pair
    # of species to be compared, and for three outgroup species. (Four
    # chains in all). TREES is a dictionary, keyed according to the
    # names of the target and output species, containing EPO and TREE
    # for each species.
    TREES = dict()
    for chain in opt.alignments:
        # Get the target species name from the file name.
        cname_parts = chain.split("/")
        cname_parts = cname_parts[-1].split(".")
        if cname_parts[0] != opt.qname:
            sys.stderr.write("Chain {} does not appear to contain the correct query species!\n".format(chain))
            exit(1)
        if cname_parts[1] not in leaves:
            sys.stderr.write("Chain {} target species not present in tree {}\n".format(chain, phylo_full.newick))
            exit(1)
            
        #loading alignments from the chain/pkl file
        EPO = dict( (ch[0].id, ch) for ch in loadChains(chain) )

        ## create an interval tree based on chain headers (from_species side)
        ## for fast feature-to-chain_header searching
        log.info("indexing %d chains ..." % (len(EPO),))
        TREE = GIntervalTree()
        for gabid in EPO:
            chain, t, q = EPO[gabid]
            TREE.add(chain.tName, Interval(chain.tStart, chain.tEnd, chain.id))

        TREES[cname_parts[1]] = {}
        TREES[cname_parts[1]]["EPO"] = EPO
        TREES[cname_parts[1]]["TREE"] = TREE

    if len(TREES) < len(leaves)-1:
        sys.stderr.write("Not enough alignments for the given tree!\n")
        exit(1)

    # transform elements
    transform_file(loadFeatures( opt.input, opt ), opt.output, TREES, leaves, phylo_full, phylo_pruned, opt)
