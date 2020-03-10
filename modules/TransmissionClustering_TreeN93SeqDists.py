#! /usr/bin/env python3
'''
Implementation of the "TransmissionClustering" module using TreeN93 using pairwise sequence distances
'''
from TransmissionClustering import TransmissionClustering
import ViReport_GlobalContext as GC
from os import makedirs
from subprocess import call

class TransmissionClustering_TreeN93SeqDists(TransmissionClustering):
    def init():
        pass

    def finalize():
        pass

    def cite():
        return GC.CITATION_TREEN93

    def blurb():
        return "Transmission clustering was performed using TreeN93 (Moshiri, 2018) using pairwise sequence distances."

    def infer_transmission_clusters():
        out_filename = '%s/transmission_clusters.tsv' % GC.OUT_DIR_OUTFILES
        treen93_dir = '%s/TreeN93' % GC.OUT_DIR_TMPFILES
        makedirs(treen93_dir)
        treen93_tree_filename = '%s/treen93.tre' % treen93_dir
        log1 = open('%s/log1.txt' % treen93_dir, 'w'); log2 = open('%s/log2.txt' % treen93_dir, 'w')
        command1 = ['TreeN93.py', '-v', '-i', GC.PAIRWISE_DISTS_SEQS, '-o', treen93_tree_filename]
        command2 = ['TreeN93_cluster.py', '-v', '-i', treen93_tree_filename, '-o', out_filename]
        f = open('%s/command.txt' % treen93_dir, 'w')
        f.write('%s\n' % ' '.join(command1))
        f.write('%s\n' % ' '.join(command2))
        f.close()
        call(command1, stderr=log1); call(command2, stderr=log2)
        log1.close(); log2.close()
        return out_filename