#! /usr/bin/env python3
'''
Implementation of the "PhylogeneticInference" module using RAxML-NG
'''
from PhylogeneticInference import PhylogeneticInference
import ViReport_GlobalContext as GC
from glob import glob
from os import makedirs
from os.path import isfile
from shutil import move
from subprocess import check_output
MODEL = {
    'DNA': 'General Time-Reversible (GTR) model (Tavare, 1986)',
    'AA': 'LG model (Le & Gascuel, 2008)',
}

class PhylogeneticInference_RAxMLNG(PhylogeneticInference):
    def init():
        pass

    def finalize():
        pass

    def cite():
        return GC.CITATION_RAXML_NG

    def blurb():
        return "A maximum-likelihood phylogeny was inferred under the %s using RAxML-NG (Kozlov et al., 2019) with GAMMA among-site rate heterogeneity (+G) and potential invariant sites (+I)." % MODEL[GC.SEQ_TYPE]

    def infer_phylogeny(aln_filename):
        if not isfile(aln_filename):
            raise ValueError("Invalid alignment file: %s" % aln_filename)
        raxmlng_dir = '%s/RAxML-NG' % GC.OUT_DIR_TMPFILES
        out_filename = '%s/unrooted.tre' % GC.OUT_DIR_OUTFILES
        if GC.GZIP_OUTPUT:
            out_filename += '.gz'
        if isfile(out_filename) or isfile('%s.gz' % out_filename):
            GC.SELECTED['Logging'].writeln("Inferred phylogeny exists. Skipping recomputation.")
        else:
            makedirs(raxmlng_dir, exist_ok=True)
            if aln_filename.lower().endswith('.gz'):
                unzipped_filename = '%s/aln_unzipped.fas' % raxmlng_dir
                GC.write_file('\n'.join(GC.read_file(aln_filename)), unzipped_filename)
                aln_filename = unzipped_filename
            command = ['raxml-ng', '--force', '--msa', aln_filename, '--model']
            if GC.SEQ_TYPE == 'DNA':
                command.append('GTR+I+G')
            elif GC.SEQ_TYPE == 'AA':
                command.append('LG+I+G')
            else:
                raise ValueError("Invalid sequence type: %s" % GC.SEQ_TYPE)
            if GC.NUM_THREADS is not None:
                command += ['--threads', str(GC.NUM_THREADS)]
            f = open('%s/command.txt' % raxmlng_dir, 'w'); f.write('%s\n' % ' '.join(command)); f.close()
            check_output(command)
            GC.write_file('\n'.join(GC.read_file('%s.raxml.bestTree' % aln_filename)), out_filename)
            for f in glob('%s.*' % aln_filename):
                move(f, '%s/%s' % (raxmlng_dir, f.split('/')[-1]))
        return out_filename
