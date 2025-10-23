#!/usr/bin/env python3
"""Go over several TOGA directories and arrange GLP table.

Table looks like:
gene species1_status species2_status etc.
"""
import argparse
import os
import sys
from collections import defaultdict

__author__ = "Bogdan Kirilenko, 2022; Anastasiya Stepanenko, 2025"

LOSS_SUMM_DATA = "loss_summary.tsv"
GENE_NAMES = "orthology_classification.tsv"
SPASS_DUMP = "/projects/hillerlab/genome/data/spass_dump.tsv"


def errprint(msg):
    sys.stderr.write(f"{msg}\n")


def parse_args():
    """Command line arguments parser."""
    app = argparse.ArgumentParser("Output goes to STDOUT")
    app.add_argument("toga_dirs_list",
                     help="Text file containing paths to "
                          "TOGA directories of interest."
                    )
    app.add_argument("--genes_list",
                     "-g",
                     help="(non-mandatory) text file containing a list of genes "
                          "of interest -> in case the full list of genes is not needed.")
    app.add_argument("--scientific_names",
                     "-s",
                     action="store_true",
                     dest="scientific_names", 
                     help="Output scientific species names."
                          "(Works in Hillelab environment only)"
                          )
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()

    if args.scientific_names and not os.path.isfile(SPASS_DUMP):
        errprint(f"Error! File {SPASS_DUMP} does not exist")
        errprint(f"Cannot execute the script with --scientific_names parameter")
        sys.exit(1)
    return args


def rm_vs_prefix(elem):
    """Remove vs_ prefix if present."""
    if elem.startswith("vs_"):
        return elem[3:]
    return elem


def read_genes_list(genes_list_arg):
    """Read genes list from file."""
    if genes_list_arg is None:
        return None
    f = open(genes_list_arg, "r")
    genes_all = [x.rstrip() for x in f.readlines()]
    ret = [x for x in genes_all if len(x) > 0]
    f.close()
    return ret


def get_sp_name_dict(sci_names_get):
    ret = {}
    if sci_names_get is False:
        # no need to open anythin
        # we don't use scientific names
        return ret
    f = open(SPASS_DUMP, "r")
    for line in f:
        ld = line.rstrip().split("\t")
        asm_name = ld[0]
        sci_name = ld[20].replace(" ", "_")
        ret[asm_name] = sci_name
    f.close()
    return ret
    

def get_toga_dirs(toga_dirs_list, sci_names=False):
    """Get TOGA dirs list."""
    f = open(toga_dirs_list, "r")
    any_non_exist = False
    any_loss_summ_data_is_absent = False
    sp_name_dict = get_sp_name_dict(sci_names)

    ret = []
    for line in f:
        path = line.rstrip()
        if not os.path.isdir(path):
            any_non_exist = True
            errprint(f"Fatal: TOGA directory {path} does not exist")
        loss_summ_path = os.path.join(path, LOSS_SUMM_DATA)
        gene_names_path = os.path.join(path, GENE_NAMES)
        if not (os.path.isfile(loss_summ_path) and os.path.isfile(gene_names_path)):
            any_loss_summ_data_is_absent = True
            errprint(f"Fatal: {LOSS_SUMM_DATA} not found in {path}")
        sp_name = rm_vs_prefix(os.path.basename(path))
        if sci_names is True:
            # weird construction but: get what is in the dictionary
            # if nothing -> return what it was
            _sp_name = sp_name_dict.get(sp_name, sp_name)
            sp_name = _sp_name
        item = (sp_name, loss_summ_path, gene_names_path)
        ret.append(item)
    if any_non_exist or any_loss_summ_data_is_absent:
        sys.exit(1)
    f.close()
    return ret


def extract_gen_stats(species_to_dir, genes_of_interest):
    ret = defaultdict(dict) # gene: sp: stat
    species_dirs = list(species_to_dir)
    gene_dict = {}
    for sp, path, genes in species_dirs:
        f = open(path, 'r')
        g = open(genes, 'r')
        for line in g:
            ld = line.rstrip().split("\t")
            ens_id = ld[0]
            gene_name = ld[1].split('#')
            if len(gene_name)>1 and ens_id not in gene_dict:
                gene_dict[ens_id] = ens_id+'#'+gene_name[1]
        for line in f:
            ld2 = line.rstrip().split("\t")
            if ld2[0] != "GENE":
                continue
            gene = ld2[1]
            stat = ld2[2]
            if genes_of_interest is None:
                if gene not in gene_dict or gene_dict[gene] is None:
                    errprint(f"Warning! No data for gene {gene}")
                    continue
                else:
                    gene_with_name = gene_dict[gene]
                    ret[gene_with_name][sp] = stat
            else:
                ret[gene][sp] = stat
        f.close()
        g.close()
    if genes_of_interest is None:
        return ret
    # ok, we have specified some genes of interest
    ret_spec = {g: ret.get(g, None) for g in genes_of_interest}
    ret_final = {}
    for k, v in ret_spec.items():
        if v is None:
            errprint(f"Warning! No data for gene {k}")
            continue
        if k in gene_dict:
            k = gene_dict[k]
        ret_final[k] = v
    return ret_final


def main():
    """Entry point."""
    args = parse_args()
    # None -> output all genes
    genes_of_interest = read_genes_list(args.genes_list)
    species_to_dir = get_toga_dirs(args.toga_dirs_list, sci_names=args.scientific_names)
    species_order = [x[0] for x in species_to_dir]
    gene_to_sp_and_stat = extract_gen_stats(species_to_dir, genes_of_interest)

    # output results
    spec_line = "\t".join(species_order)
    print(f"GENE\t{spec_line}")
    for gene, sp_stat in gene_to_sp_and_stat.items():
        stats_ordered = []
        for sp in species_order:
            stat = sp_stat.get(sp, "N")
            stats_ordered.append(stat)
            if stat == "N":
                errprint(f"Warning! Cannot find stat for gene {gene} in {sp}")
        stats_str = "\t".join(stats_ordered)
        print(f"{gene}\t{stats_str}")
    #print(gene_to_sp_and_stat)

if __name__ == '__main__':
    main()
