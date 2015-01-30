#!/usr/bin/env python2.7
"""
CGquery handler

given arguments, creates a pickled dictionary mapping TCGA analysis IDs to full cghub bam slicer
query strings. These strings are used in a cghub bam slicer query to extract regions from the bam
they represent.
"""

import sys, os, argparse, urllib2, itertools
import xml.etree.cElementTree as ET
import cPickle as pickle

from lib.general_lib import FullPaths

#this hard-coded list stores which assemblies use '1' instead of 'chr1'
#because people suck
non_chr_names = ['GRCh37-lite','GRCh37']

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", "-g", nargs="+", required=True, 
            help="Genomes to send to cgquery.")
    parser.add_argument("--tissue_types", "-u", nargs="+", required=True, 
            help="Tissue types to look at.")
    parser.add_argument("--size", nargs=2, default=["100000000000","*"], 
            help="two integer values representing # of bytes we want to look at. Second value can be *")
    parser.add_argument("--study", default="phs000178", 
            help="Study to look at. default=phs000178.")
    parser.add_argument("--library_strategy", default="WGS", 
            help="library type. default=WGS")
    parser.add_argument("--target_range", "-t", nargs="+", default=["chr1:120309986-145388461"], 
            help="target range(s) of assembly. default = chr1:120309986-145388461. include the word chr.")
    parser.add_argument("--out", "-o", default="queries/queries.pickle", type=argparse.FileType("wb"), 
            action=FullPaths, 
            help="Location to write out the pickled queries to. Default = queries/queries.pickle")
    return parser.arse_args()


def search_metadata(genomes, tissue_types, size, study, library_strategy):
    """
    given arguments to program, query cghub for all xmls.
    returns a list of xmls representing each combination of queries.
    """
    size = '[' + "%20TO%20".join(size) + ']'
    xml_list = []
    base_path = "https://cghub.ucsc.edu/cghub/metadata/analysisDetail?"
    for genome, tissue in itertools.product(genomes, tissue_types):
        search = ["=".join([x,y]) for x,y in zip(["study","sample_type","library_strategy",
                "refassem_short_name","filesize"],[study,tissue,library_strategy,genome,size])]
        url = base_path+"&".join(search)
        xml_list.append(ET.parse(urllib2.urlopen(url)).getroot())
    #return a flat list of xml elements
    return [item for sublist in xml_list for item in sublist]


def parse_metadata_xml(xml_list):
    """
    Parses a metadata xml retrieved from cgquery.
    """
    analysis_dict = {}
    for item in xml_list:
        if item.find('analysis_id') is not None and item.find('state').text != 'suppressed':
            analysis = item.find('analysis_id').text
            refassem_short_name = item.find('refassem_short_name').text
            analysis_dict[analysis] = refassem_short_name
    return analysis_dict


def build_bamslicer_query(analysis, refassem_short_name, ranges):
    """
    Builds a query string for the CGhub bam slicer.
    """
    if refassem_short_name in non_chr_names:
        ranges = [x.replace("chr","") for x in ranges]
    formatted_ranges = "".join(["&range=", "&range=".join(ranges)]) 
    base_path = "https://slicer.cghub.ucsc.edu/analyses/{}/slices?ref={}&format=bam".format(analysis, refassem_short_name)
    query_string = "".join([base_path, formatted_ranges])
    return query_string


def main():
    args = parse_args(args)
    xml_list = search_metadata(args.genomes, args.tissue_types, args.size, args.study, args.library_strategy)
    analysis_dict = parse_metadata_xml(xml_list)
    
    query_dict = {}
    for analysis, refassem_short_name in analysis_dict.iteritems():
        query_string = build_bamslicer_query(analysis, refassem_short_name, args.target_range)

    pickle.dump(args.out, query_dict)


if __name__ == '__main__':
    sys.exit(main())
