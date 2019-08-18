#!/usr/bin/env python3

import os
import json
import itertools
from shutil import copyfile

try:
    SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()

def main():

    data_src_dir        = os.path.join(SCRIPT_DIR, "data")
    sampled_dist_dir    = os.path.join(SCRIPT_DIR, "sampled_data")
    contigs_src_dir     = os.path.join(data_src_dir, "contigs")
    contigs_dist_dir    = os.path.join(sampled_dist_dir, "contigs")
    experiment_src_dir  = os.path.join(data_src_dir, "experiment")
    experiment_dist_dir = os.path.join(sampled_dist_dir, "experiment")

    def mkdir(path):
        if not os.path.exists(path):
            os.mkdir(path)
    mkdir(sampled_dist_dir)
    mkdir(contigs_dist_dir)
    mkdir(experiment_dist_dir)

    contigs_with_experiment = json.load(
        open(os.path.join(SCRIPT_DIR, "data/contigs_with_experiment.json")))

    # First 2000 contigs with experiment data
    sampled_contigs = set(contigs_with_experiment[:2000])

    ## Copy index.json
    src_index = json.load(open(os.path.join(data_src_dir, "index.json")))
    new_index = []
    for item in src_index:
        if item['contig_id'] in sampled_contigs:
            new_index.append(item)
    json.dump(new_index, open(os.path.join(sampled_dist_dir, "index.json"), "w"))
    json.dump(new_index, open(os.path.join(sampled_dist_dir, "top_abundance_index.json"), "w"))

    ## Copy experiment_info.json
    copyfile(os.path.join(data_src_dir, "experiment_info.json"),
             os.path.join(sampled_dist_dir, "experiment_info.json"))

    ## Copy contigs and experiment
    for contig_id in sampled_contigs:
        copyfile(os.path.join(contigs_src_dir, contig_id + ".json"),
                 os.path.join(contigs_dist_dir, contig_id + ".json"))
    for contig_id in sampled_contigs:
        copyfile(os.path.join(experiment_src_dir, contig_id + ".json"),
                 os.path.join(experiment_dist_dir, contig_id + ".json"))

if __name__ == "__main__":
    main()


