import argparse
from multiprocessing import cpu_count
from os import environ, makedirs, getenv
from os.path import abspath, expanduser, isdir,isfile, join
from shutil import copyfile, rmtree
from sys import argv, path, modules
WS_HOME = getenv('WS_HOME')
path.insert(0, join(WS_HOME, 'pasta/'))
path.insert(0, join(WS_HOME, 'sepp/'))
import pandas as pd
from pasta import *
import os
from pasta.alignment import *
from glob import glob
import numpy as np
from tempfile import TemporaryDirectory, TemporaryFile
import dendropy
from pasta.tree import PhylogeneticTree
import subprocess
from pasta.new_decomposition import *
from pasta.pastaalignerjob import bisect_tree
from subprocess import Popen, PIPE
import multiprocessing
from copy import deepcopy
from time import sleep

def decompose_phylogeny(phy, max_size, min_size):
    trees_map = []
    tree_list = [phy]
    while len(tree_list) > 0:
        tmp_phy = tree_list.pop()
        t1, t2 = bisect_tree(tree=tmp_phy, breaking_edge_style='midpoint', max_size=max_size)
        if t1.count_leaves() > min_size:
            tree_list.append(deepcopy(t1))
        else:
            trees_map.append(deepcopy(t1))
        if t2.count_leaves() > min_size:
            tree_list.append(deepcopy(t2))
        else:
            trees_map.append(deepcopy(t2))
    return trees_map


def run_process(command):
    p = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    input_fp = command[-1]
    output_fp = command[-1].replace('fa', 'aln')
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    rc = p.returncode
    with open(output_fp, 'wb') as f:
        f.write(output)
    f.close()

    aligned_core_ca = CompactAlignment()
    aligned_core_ca.read_filepath(output_fp)
    aligned_core_ca.write_filepath(output_fp)
    return rc, err, output_fp


def run_alignment(process_name, tasks, results):
    command = tasks.get()
    while True:
        if command == -1:
            print('[%s] evaluation routine quits' % process_name)

            # Indicate finished
            results.put(-1)
            break
        else:
            rc, err, output_fp = run_process(command)
            results.put((rc, err, output_fp))
    return


if __name__ == '__main__':
    # main_dir = '/Users/esayyari/UCSD/oasis/viruses/all/corona_overall.April7/'
    main_dir = argv[1]
    num_processes = int(argv[2])
    max_size = int(argv[3])
    min_size = int(argv[4])

    ca = CompactAlignment()
    ca.read_filepath(join(main_dir, 'dna-sequences.fasta'))
    print("The total number of taxa is", ca.get_num_taxa())
    tree = dendropy.Tree.get(path=join(main_dir, 'sequences', 'fastme_tree.nwk'), schema="newick")
    phy = PhylogeneticTree(tree)
    orig_phy = deepcopy(phy)
    trees_map = decompose_phylogeny(phy, max_size=max_size, min_size=min_size)
    core_ca = CompactAlignment()
    core_ca.read_filepath(join(main_dir, 'dna-sequences.core.fasta'))
    IDs = list(core_ca.keys())
    tmp_dir = join(main_dir, 'sub_process')
    i = 0
    commands = []
    for tmp_tre in trees_map:
        i += 1
        command = ['mafft', '--reorder', '--nomemsave', '--thread', '1', '--auto',
                   join(tmp_dir, 'dna-sequences.' + str(i) + '.fa')]
        commands.append(command)
    core_ca.clear()
    commands_to_write = []
    for command in commands:
        commands_to_write.append(" ".join(command))

    with open(join(main_dir, 'subalignment_jobs.txt'), 'w') as g:
        g.write("\n".join(commands_to_write) + "\n")
    print("The alignment jobs are written to file", join(main_dir, 'subalignment_jobs.txt'))
    print("All jobs are set up to run!")
    print("The temporary directory is", tmp_dir)
    # Define IPC manager
    manager = multiprocessing.Manager()

    # Define a list (queue) for tasks and computation results
    tasks = manager.Queue()
    results = manager.Queue()
    # Create process pool with four processes
    pool = multiprocessing.Pool(processes=num_processes)
    processes = []
    # Initiate the worker processes
    for i in range(num_processes):
        # Set process name
        process_name = 'P%i' % i

        # Create the process, and connect it to the worker function
        new_process = multiprocessing.Process(target=run_alignment, args=(process_name, tasks, results))

        # Add new process to the list of processes
        processes.append(new_process)

        # Start the process
        new_process.start()
    # Fill task queue
    all_tasks = []
    for i, command in enumerate(commands):
        if i == num_processes:
            break
        all_tasks.append(command)
    sleep(5)
    for i in range(num_processes):
        all_tasks.append(-1)

    for single_task in all_tasks:
        tasks.put(single_task)
    print("The total number of tasks is:", len(all_tasks))
    print("The number of subprocesses is", num_processes)

    num_finished_processes = 0

    while True:
        # Read result
        new_result = results.get()

        # Have a look at the results
        if new_result == -1:
            # Process has finished
            num_finished_processes += 1

            if num_finished_processes == num_processes:
                break
        else:
            rc, err, output_fp = new_result
            print("Output file path", output_fp, "generated with the output code", rc)
            if rc != 0:
                print("Something went wrong running", output_fp)
                print(err)
