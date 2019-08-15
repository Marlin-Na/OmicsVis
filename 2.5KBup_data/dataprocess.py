#!/usr/bin/env python3

import os
import json
import csv
import itertools
from collections import namedtuple
from Bio import SeqIO

try:
    SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()

GENE_FASTA_FILE = os.path.join(SCRIPT_DIR, "raw/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline_nucl.fasta")
CONTIG_FILE = os.path.join(SCRIPT_DIR, "raw/KansasContig_length_abundance_sorted.txt")
WD_PROTEIN_FILE = os.path.join(SCRIPT_DIR, "raw/Kansas_WD_2.5KBup_QValue0p01_2019Aug01_Protein_SpectraCount_crosstab.csv")

GENE_RECORD_FIELDS = [
    "contig_id","gene_start","gene_end","gene_strand",
    "gene_id","gene_partial","gene_start_type","gene_rbs_motif","gene_rbs_spacer",
    "gene_gc_cont"
]

CONTIG_RECORD_FIELDS = ["contig_id", "contig_length", "contig_abundance", "genes"]

def fixeddict(fields):
    def ans(**args):
        userans = dict()
        for key, value in args.items():
            if key in fields:
                userans[key] = value
            else:
                raise Exception("key {} does not match".format(key))
        return userans
    return ans

def create_contig_map(file):
    ContigRecord = fixeddict(CONTIG_RECORD_FIELDS)
    ans = dict()
    with open(file, "r") as f:
        for line in f:
            contig_id, contig_length, contig_abundance = line.rstrip().split()
            contig_id = contig_id.split("_")[1]
            contig_length = int(contig_length)
            contig_abundance = float(contig_abundance)
            ans[contig_id] = ContigRecord(
                contig_id = contig_id,
                contig_length = contig_length,
                contig_abundance = contig_abundance,
                genes = []
            )
    return ans

def read_gene_records(fasta):
    GeneRecord = fixeddict(GENE_RECORD_FIELDS)
    fasta = SeqIO.parse(fasta, "fasta")
    for record in fasta:
        words = record.description.split()
        ans = GeneRecord(
            contig_id    = words[0].split("_")[0],
            gene_start   = int(words[2]),
            gene_end     = int(words[4]),
            gene_strand  = "+" if words[6] == "1" else "-",
            #gene_id         = words[8].split(";")[0].split("=")[1],
            gene_id         = words[0],
            gene_partial    = words[8].split(";")[1].split("=")[1],
            gene_start_type = words[8].split(";")[2].split("=")[1],
            gene_rbs_motif  = words[8].split(";")[3].split("=")[1],
            gene_rbs_spacer = words[8].split(";")[4].split("=")[1],
            gene_gc_cont    = float(words[8].split(";")[5].split("=")[1]),
        )
        yield ans


## Write all the needed files for visualization
def write_contig_map(contig_map):
    data_dir = os.path.join(SCRIPT_DIR, "data")
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    if not os.path.exists(os.path.join(data_dir, "contigs")):
        os.makedirs(os.path.join(data_dir, "contigs"))

    ## Write index.json and top_abundance_index.json
    indexobj = []
    for contig in contig_map.values():
        indexobj.append({"contig_id"        : contig['contig_id'],
                         "contig_length"    : contig['contig_length'],
                         "contig_abundance" : contig['contig_abundance']})
    with open(os.path.join(data_dir, "index.json"), "w") as f:
        json.dump(indexobj, f, indent = 1)
    with open(os.path.join(data_dir, "top_abundance_index.json"), "w") as f:
        indexobj_sorted = sorted(indexobj,
                                 key = lambda x: x["contig_abundance"], reverse = True)
        json.dump(indexobj_sorted[:10000], f, indent = 1)

    ## Write json for each contig
    for contig in contig_map.values():
        contig_id = contig['contig_id']
        file_name = os.path.join(data_dir, "contigs", contig_id + ".json")
        with open(file_name, "w") as f:
            json.dump(contig, f, indent = 1)

## Put one gene record into the contig_map
def put_gene_record(contig_map, gene_record):
    if gene_record['contig_id'] not in contig_map:
        print("Contig {} not found for gene {}".format(
            gene_record['contig_id'], gene_record['gene_id']))
        return
    contig_map[gene_record['contig_id']]['genes'].append(gene_record)
    return

def process_wd_protein_data():
    data_dir = os.path.join(SCRIPT_DIR, "data")
    def read_wd_protein_data():
        with open(WD_PROTEIN_FILE) as csvfile:
            reader = csv.reader(csvfile)
            first_row = True
            for row in reader:
                if first_row:
                    first_row = False
                    continue
                gene_id = row[0]
                contig_id = gene_id.split("_")[0]
                wd_data = row[3:30]
                for i in range(len(wd_data)):
                    value = None if wd_data[i] == "" else int(float(wd_data[i]))
                    wd_data[i] = value
                yield {
                    "contig_id": contig_id,
                    "gene_id": gene_id,
                    "experiment_data": wd_data
                }

    # Write each json file
    if not os.path.exists(os.path.join(data_dir, "experiment")):
        os.makedirs(os.path.join(data_dir, "experiment"))

    grouped_by_contig = itertools.groupby(
        sorted(read_wd_protein_data(), key=lambda x: x['contig_id']),
        lambda x: x['contig_id']
    )

    for contig_id, experiments in grouped_by_contig:
        with open(os.path.join(data_dir, "experiment", contig_id + ".json"), "w") as dist_file:
            json.dump({"contig_id": contig_id, "experimentsbygene": list(experiments)},
                      dist_file, indent = 1)

    # Write experiment_info.json
    with open(os.path.join(data_dir, "experiment_info.json"), "w") as dist_file:
        ## TODO: use correct field name
        field_names = [ "metric" + str(i+1) for i in range(0,27) ]
        obj = {"field_names": field_names}
        json.dump(obj, dist_file, indent = 1)


def main():
    if __name__ != "__main__":
        return
    contig_map = create_contig_map(CONTIG_FILE)
    for record in read_gene_records(GENE_FASTA_FILE):
        put_gene_record(contig_map, record)
    write_contig_map(contig_map)
    process_wd_protein_data()

main()



# ‘2.5KBup’ data:
# Contig files:
# /pic/dtn2/wuru978/Kansas_contig_gene_forJialin/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline.fasta          (no need)
# Gene calling file:
# /pic/dtn2/wuru978/Kansas_contig_gene_forJialin/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline_prodigal_out   (no need)
# /pic/dtn2/wuru978/Kansas_contig_gene_forJialin/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline_score          (no need?)
# /pic/dtn2/wuru978/Kansas_contig_gene_forJialin/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline_nucl.fasta
# /pic/dtn2/wuru978/Kansas_contig_gene_forJialin/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline_prot.fasta     (no need)

#/pic/dtn2/wuru978/Kansas_contig_gene_forJialin/LengthAppend_toHeader_seqlength_seqname_seq_KS_ctg_2500_oneline_nucl.fasta

