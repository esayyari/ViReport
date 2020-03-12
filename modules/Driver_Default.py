#! /usr/bin/env python3
'''
Default implementation of the "Driver" module
'''
from Driver import Driver
import ViReport_GlobalContext as GC
from os import makedirs
from os.path import isfile
from sys import stdout

def print_message():
    '''
    Print author message
    '''
    title = "ViReport v%s" % GC.VIREPORT_VERSION
    devel = "Niema Moshiri 2020"
    l = max(len(title),len(devel))
    print("/-%s-\\" % ('-'*l))
    for e in (title,devel):
        lpad = int((l-len(e))/2)
        rpad = l - lpad - len(e)
        print("| %s%s%s |" % (lpad*' ',e,rpad*' '))
    print("\\-%s-/\n" % ('-'*l))

class Driver_Default(Driver):
    def init():
        pass

    def cite():
        return GC.CITATION_VIREPORT

    def run(seqs_filename, sample_times_filename, outgroups_filename, categories_filename):
        # organize citations
        GC.CITATIONS = set()
        for m in GC.SELECTED:
            cite = GC.SELECTED[m].cite()
            if isinstance(cite,str):
                GC.CITATIONS.add(cite.strip())
            elif isinstance(cite,set) or isinstance(cite,list):
                for c in cite:
                    GC.CITATIONS.add(c.strip())
        GC.CITATIONS = sorted(GC.CITATIONS)

        # print starting messages
        print_message()
        print("========================   Workflow Process   ========================")
        print("ViReport was run as follows: %s" % GC.VIREPORT_COMMAND)
        print("Output directory: %s" % GC.OUT_DIR_PRINT)
        print("Starting viral analysis workflow...")

        # check input files
        if not isfile(seqs_filename):
            raise ValueError("Invalid sequence file: %s" % seqs_filename)
        GC.INPUT_SEQS = seqs_filename
        if not isfile(sample_times_filename):
            raise ValueError("Invalid sample times file: %s" % sample_times_filename)
        GC.INPUT_TIMES = sample_times_filename
        if outgroups_filename is not None and not isfile(outgroups_filename):
            raise ValueError("Invalid outgroups list file: %s" % outgroups_filename)
        GC.INPUT_OUTGROUPS = outgroups_filename
        if categories_filename is not None and not isfile(categories_filename):
            raise ValueError("Invalid sample categories file: %s" % categories_filename)
        GC.INPUT_CATEGORIES = categories_filename

        # set up output and intermediate folders
        GC.OUT_DIR_OUTFILES = "%s/output_files" % GC.OUT_DIR
        makedirs(GC.OUT_DIR_OUTFILES, exist_ok=True)
        GC.OUT_DIR_TMPFILES = "%s/intermediate_files" % GC.OUT_DIR
        makedirs(GC.OUT_DIR_TMPFILES, exist_ok=True)
        GC.OUT_DIR_REPORTFILES = "%s/report_files" % GC.OUT_DIR
        makedirs(GC.OUT_DIR_REPORTFILES, exist_ok=True)
        GC.OUT_DIR_REPORTFIGS = '%s/figs' % GC.OUT_DIR_REPORTFILES
        makedirs(GC.OUT_DIR_REPORTFIGS, exist_ok=True)

        # initialize module implementations
        print("\nInitializing module implementations...", end=' '); stdout.flush()
        for k in GC.SELECTED:
            GC.SELECTED[k].init()
        print("done")

        # run preprocessing
        print("\nRunning '%s'..." % GC.SELECTED['Preprocessing'].__name__); stdout.flush()
        GC.PROCESSED_SEQS, GC.PROCESSED_TIMES, GC.PROCESSED_OUTGROUPS, GC.PROCESSED_CATEGORIES = GC.SELECTED['Preprocessing'].preprocess(GC.INPUT_SEQS, GC.INPUT_TIMES, GC.INPUT_OUTGROUPS, GC.INPUT_CATEGORIES)
        GC.SEQ_TYPE = GC.predict_seq_type(GC.PROCESSED_SEQS)
        print("Preprocessed sequences output to: %s" % GC.PROCESSED_SEQS)
        print("Preprocessed sample times output to: %s" % GC.PROCESSED_TIMES)
        if GC.PROCESSED_OUTGROUPS is not None:
            print("Preprocessed outgroups list output to: %s" % GC.PROCESSED_OUTGROUPS)
        if GC.PROCESSED_CATEGORIES is not None:
            print("Preprocessed sample categories output to: %s" % GC.PROCESSED_CATEGORIES)

        # align the preprocessed sequences
        print("\nRunning '%s'..." % GC.SELECTED['MultipleSequenceAlignment'].__name__); stdout.flush()
        GC.ALIGNMENT_WITH_OUTGROUP = GC.SELECTED['MultipleSequenceAlignment'].align(GC.PROCESSED_SEQS)
        GC.ALIGNMENT = GC.remove_outgroups_fasta(GC.ALIGNMENT_WITH_OUTGROUP, GC.PROCESSED_OUTGROUPS)
        print("Multiple sequence alignment output to: %s" % GC.ALIGNMENT)

        # compute pairwise sequence distances
        print("\nRunning '%s'..." % GC.SELECTED['PairwiseDistancesSequence'].__name__); stdout.flush()
        GC.PAIRWISE_DISTS_SEQS = GC.SELECTED['PairwiseDistancesSequence'].pairwise_distances(GC.ALIGNMENT)
        print("Pairwise sequence distances output to: %s" % GC.PAIRWISE_DISTS_SEQS)

        # infer a phylogeny
        print("\nRunning '%s'..." % GC.SELECTED['PhylogeneticInference'].__name__); stdout.flush()
        GC.TREE_UNROOTED_WITH_OUTGROUP = GC.SELECTED['PhylogeneticInference'].infer_phylogeny(GC.ALIGNMENT_WITH_OUTGROUP)
        GC.TREE_UNROOTED = GC.remove_outgroups_newick(GC.TREE_UNROOTED_WITH_OUTGROUP, GC.PROCESSED_OUTGROUPS)
        print("Inferred (unrooted) phylogeny output to: %s" % GC.TREE_UNROOTED)

        # compute pairwise phylogenetic distances
        print("\nRunning '%s'..." % GC.SELECTED['PairwiseDistancesTree'].__name__); stdout.flush()
        GC.PAIRWISE_DISTS_TREE = GC.SELECTED['PairwiseDistancesTree'].pairwise_distances(GC.TREE_UNROOTED)
        print("Pairwise phylogenetic distances output to: %s" % GC.PAIRWISE_DISTS_TREE)

        # root the phylogeny
        print("\nRunning '%s'..." % GC.SELECTED['Rooting'].__name__); stdout.flush()
        GC.TREE_ROOTED_WITH_OUTGROUP = GC.SELECTED['Rooting'].root(GC.TREE_UNROOTED_WITH_OUTGROUP)
        GC.TREE_ROOTED = GC.remove_outgroups_newick(GC.TREE_ROOTED_WITH_OUTGROUP, GC.PROCESSED_OUTGROUPS)
        print("Rooted phylogeny output to: %s" % GC.TREE_ROOTED)

        # date the rooted phylogeny
        print("\nRunning '%s'..." % GC.SELECTED['Dating'].__name__); stdout.flush()
        GC.TREE_DATED = GC.SELECTED['Dating'].date(GC.TREE_ROOTED, GC.PROCESSED_TIMES)
        print("Dated phylogeny output to: %s" % GC.TREE_DATED)

        # infer ancestral sequence(s)
        print("\nRunning '%s'..." % GC.SELECTED['AncestralSequenceReconstruction'].__name__); stdout.flush()
        GC.ANCESTRAL_SEQS = GC.SELECTED['AncestralSequenceReconstruction'].reconstruct(GC.TREE_ROOTED, GC.ALIGNMENT)
        print("Ancestral sequence(s) output to: %s" % GC.ANCESTRAL_SEQS)

        # perform transmission clustering
        print("\nRunning '%s'..." % GC.SELECTED['TransmissionClustering'].__name__); stdout.flush()
        GC.TRANSMISSION_CLUSTERS = GC.SELECTED['TransmissionClustering'].infer_transmission_clusters()
        print("Transmission clusters output to: %s" % GC.TRANSMISSION_CLUSTERS)

        # write the report
        print("\nWriting report using '%s'..." % GC.SELECTED['WriteReport'].__name__); stdout.flush()
        GC.REPORT = GC.SELECTED['WriteReport'].write_report()
        print("Report written to: %s" % GC.REPORT)

        # print citations
        print("\n\n===========================   Citations   ============================")
        for cite in GC.CITATIONS:
            print(cite)
