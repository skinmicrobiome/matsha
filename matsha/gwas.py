import logging
import warnings
from scipy import stats
import numpy as np
import os
from collections import defaultdict
from matsha import utils


def gwas_correlation(genotype_dct, phenotype_dct, min_maf, min_sample_size):
    stats_lst = []
    for v in genotype_dct:
        geno_lst = []
        pheno_lst = []
        for sample in phenotype_dct:
            if sample in genotype_dct[v]:
                if genotype_dct[v][sample] == -1:
                    continue
                else:
                    geno_lst.append(genotype_dct[v][sample])
                    pheno_lst.append(float(phenotype_dct[sample]))
        valid_counts = len(geno_lst)
        try:
            if min(Counter(np.array(geno_lst)>0).values()) >= valid_counts*min_maf and valid_counts >= min_sample_size:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    S_corr, S_p = stats.spearmanr(geno_lst, pheno_lst)
                stats_lst.append((v, 'correlation', S_corr, S_p))
        except:
            continue
    
    return stats_lst


def gwas_binary(genotype_dct, phenotype_dct, min_maf, min_sample_size):
    stats_lst = []
    uniq_phenotype = set(phenotype_dct.values())
    for v in genotype_dct:
        contingency_dct = {p:{'g0':0, 'g1':0} for p in uniq_phenotype}
        geno_count_dct = {'g0':0, 'g1':0}
        for sample in phenotype_dct:
            if sample in genotype_dct[v]:
                if genotype_dct[v][sample] == -1:
                    continue
                elif genotype_dct[v][sample] == 0:
                    contingency_dct[phenotype_dct[sample]]['g0'] += 1
                    geno_count_dct['g0'] += 1
                else:
                    contingency_dct[phenotype_dct[sample]]['g1'] += 1
                    geno_count_dct['g1'] += 1
        valid_counts = sum(geno_count_dct.values())
        try:
            if min(geno_count_dct.values()) >= valid_counts*min_maf and valid_counts >= min_sample_size:
                contingency_tb = [list(i.values()) for i in contingency_dct.values()]
                odds_ratio, p_value = stats.fisher_exact(contingency_tb, alternative='two-sided')
                stats_lst.append((v, 'enrichment', odds_ratio, p_value))
        except:
            continue
    return stats_lst


def merge_variants_by_gene(genotype_by_variant, genes_lst):
    genotype_by_gene = defaultdict(dict)
    for v, samples in genotype_by_variant.items():
        contig, pos = v.split(':')
        pos = int(pos) - 1
        for g in genes_lst:
            if g[0] == contig and g[1] < pos < g[2]:
                for sample, val in samples.items():
                    genotype_by_gene[g[3]][sample] = min(1, max(val, genotype_by_gene[g[3]].get(sample, 0)))
    return genotype_by_gene


def save_gwas_results(stats, output_file):
    """Save GWAS results."""
    with open(output_file, 'w') as wf:
        wf.write('\t'.join(['Phenotype', 'Hit', 'Test', 'Statistic', 'p', 'q'])+'\n')
        for l in stats:
            wf.write('\t'.join(str(x) for x in l) + '\n')


def run_GWAS(genotype_dct, phenotype_dct, genes_lst, mode, min_maf, min_sample_size, qcutoff, output_dir):
    """GWAS test."""
    for gwas_mode in ["gene", "variant"]:
        if mode not in [gwas_mode, "both"]:
            continue

        # Prepare genotypes
        if gwas_mode == "gene":
            if len(genes_lst) > 0:
                genotypes = merge_variants_by_gene(genotype_dct, genes_lst)
            else:
                logging.info(f"No gene annotation found; skipping GWAS ({gwas_mode})...")
                continue
        else:
            genotypes = {}
            # Attach gene info if available
            for var in genotype_dct:
                contig, pos = var.split(':')
                pos = int(pos) - 1
                for g in genes_lst:
                    if g[0] == contig and g[1] < pos < g[2]:
                        new_var = f"{var}|{g[3]}"
                        break
                else:
                    new_var = f"{var}|"
                genotypes[new_var] = genotype_dct[var]

        # Run statistical tests
        logging.info(f"Running GWAS ({gwas_mode})...")
        all_stats_lst = []
        for ip in phenotype_dct:
            cur_phenotype_dct = phenotype_dct[ip]
            uniq_phenotype = set(cur_phenotype_dct.values())
            if len(uniq_phenotype) > 2: # quantitative variable
                logging.info(f"Treating phenotype {ip} as quantitative variable.")
                stats_lst = gwas_correlation(genotypes, cur_phenotype_dct, min_maf, min_sample_size)
            elif len(uniq_phenotype) == 2: # binary variable
                logging.info(f"Treating phenotype {ip} as binary variable.")
                stats_lst = gwas_binary(genotypes, cur_phenotype_dct, min_maf, min_sample_size)
            else:
                logging.info(f"Skipping phenotype {ip} due to having fewer than two classes.")
                continue
            logging.info(f"GWAS completed. {len(stats_lst)} tests were done.")
            pval_lst = [s[-1] for s in stats_lst]
            qval_lst = utils.fdr(np.array(pval_lst))
            cur_stats_lst = [[ip]+[i for i in stat]+[q] for stat, q in zip(stats_lst, qval_lst)]
            all_stats_lst.extend(cur_stats_lst)

        # Save results
        save_gwas_results(all_stats_lst, os.path.join(output_dir, f"gwas_output.{gwas_mode}.tsv"))
        save_gwas_results([s for s in all_stats_lst if s[-1] < qcutoff], os.path.join(output_dir, f"gwas_output_significant.{gwas_mode}.tsv"))
        