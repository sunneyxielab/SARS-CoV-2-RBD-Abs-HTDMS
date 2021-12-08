import pandas as pd
from scipy.stats import binom_test, poisson
import numpy as np
import argparse
import sys, os
import Bio.SeqIO
import time
import dms_variants.globalepistasis
import dms_variants.codonvarianttable
import dms_variants.binarymap
import collections, itertools
from lets_plot import GGBunch, ylim, ggplot, ggsave, aes, geom_boxplot, geom_point, geom_text, scale_x_discrete,ggtitle

parser = argparse.ArgumentParser(description="Step 2: calculate escape score for every barcode.")

parser.add_argument('--reference_variant_counts', '-r', type=str)
parser.add_argument('--escape_variant_counts', '-e', type=str)
parser.add_argument('--pos_or_neg', '-d', type=str)
parser.add_argument('--experiment', '-exp', type=str, default='MACS')
parser.add_argument('--wildtype', '-w', type=str)
parser.add_argument('--table', '-t', type=str)
parser.add_argument('--frac_escape', type=float, default = None)
parser.add_argument('--cells_sorted', type=int, default = None)
parser.add_argument('--output_dir', '-o', type=str)
parser.add_argument('--primary_target', '-p', type=str, default='SARS-CoV-2')
parser.add_argument('--wtseqstart', '-S', type=int, default=330)
parser.add_argument('--mutbindexpr', type=str)
parser.add_argument('--variant_bind', type=str)
parser.add_argument('--variant_expr', type=str)
parser.add_argument('--exprmin', type=float, default=-1.0)
parser.add_argument('--bindmin', type=float, default=-100.0)
parser.add_argument("--score", type=str, default="ratio")
parser.add_argument('--site_average', type=str, default='naive') # naive, topN, blosum

def site_average(x, method, flag, plot_max):
    if method == 'naive':
        return min(plot_max, x.mean(skipna = True))
    elif method == 'topN':
        N = 10
        if sum(~pd.isna(x)) <= N:
            return x.mean(skipna = True)
        if flag == 'pos':
            tmp = sorted(x.fillna(100000))
            return min(plot_max, np.mean(tmp[0:N]))
        elif flag == 'neg':
            tmp = sorted(x.fillna(-100000), key=lambda k: -k)
            return min(plot_max, np.mean(tmp[0:N]))
        else:
            raise ValueError


def build_variant_table(table_path, wt_seq_path):
    wt_seqrecord = Bio.SeqIO.read(wt_seq_path, 'fasta')
    geneseq = str(wt_seqrecord.seq)
    primary_target = wt_seqrecord.name
    variants = dms_variants.codonvarianttable.CodonVariantTable(
               geneseq = geneseq,
               barcode_variant_file = table_path,
               substitutions_are_codon = True,
               substitutions_col = 'codon_substitutions',
               primary_target = primary_target)
    return variants, primary_target


def print_some_logs(args):
    sys.stdout.write('Reference: ' + args.reference_variant_counts + '\n')
    sys.stdout.write('Escape: ' + args.escape_variant_counts + '\n')
    sys.stdout.write('Experiment: ' + args.experiment + '\t' + args.pos_or_neg + '\n')
    if args.experiment == 'FACS':
        sys.stdout.write('FACS result: ' + str(args.cells_sorted) + ' cells , and ' + str(args.frac_escape) + ' escaped\n')

def filter_expr_bind(escape_scores, mutbindexpr, variant_bind, variant_expr, min_expr_mut, min_expr_variant, min_bind_mut, min_bind_variant, filter_on_variants=False):
    # filter on mutations
    mut_bind_expr = pd.read_csv(mutbindexpr)
    assert mut_bind_expr['mutation_RBD'].nunique() == len(mut_bind_expr)
    for prop in ['bind', 'expr']:
        muts_adequate = set(mut_bind_expr
                            .query(f"{prop}_avg >= {(min_expr_mut if prop == 'expr' else min_bind_mut)}")
                            ['mutation_RBD']
                            )
        sys.stdout.write(str(len(muts_adequate))+" of "+str(len(mut_bind_expr))+" mutations have adequate "+prop+"\n")
        escape_scores[f"muts_pass_{prop}_filter"] = (
            escape_scores
            ['aa_substitutions']
            .map(lambda s: set(s.split()).issubset(muts_adequate))
            ) 
    if filter_on_variants:    
        # filter on variants
        for prop, col in [('bind', 'delta_log10Ka'), ('expr', 'delta_ML_meanF')]:
            filter_name = f"variant_pass_{prop}_filter"
            variant_pass_df = (
                pd.read_csv((variant_bind if prop == 'bind' else variant_expr), keep_default_na=False, na_values=['NA'])
                .groupby(['library', 'target', 'barcode'])
                .aggregate(val=pd.NamedAgg(col, 'mean'))
                .reset_index()
                .assign(pass_filter=lambda x: x['val'] >= (min_expr_variant if prop == 'expr' else min_bind_variant))
                .rename(columns={'pass_filter': filter_name,
                                'val': f"variant_{prop}"})
                )
            escape_scores = (
                escape_scores
                .drop(columns=filter_name, errors='ignore')
                .merge(variant_pass_df,
                    how='left',
                    validate='many_to_one',
                    on=['library', 'target', 'barcode'],
                    )
                )
            assert escape_scores[filter_name].notnull().all()
    else:
        escape_scores = (
            escape_scores.assign(variant_pass_bind_filter=True)
            .assign(variant_pass_expr_filter=True)
        )
    # annotate as passing overall filter if passes all mutation and binding filters:
    escape_scores['pass_ACE2bind_expr_filter'] = (
            escape_scores['muts_pass_bind_filter'] &
            escape_scores['muts_pass_expr_filter'] &
            escape_scores['variant_pass_bind_filter'] &
            escape_scores['variant_pass_expr_filter']
            )
    
    return escape_scores.query('pass_ACE2bind_expr_filter == True')
    # return escape_scores

def main():
    bunch = GGBunch()
    plot_max = 1.0
    protein_chain = 'E'    
    args = parser.parse_args()
    SCORE = args.score
    
    assert (args.experiment == 'MACS') or (args.experiment == 'FACS' and (args.frac_escape is not None) and (args.cells_sorted is not None)), 'experiment argument error\n'

    ref = pd.read_csv(args.reference_variant_counts)
    escape = pd.read_csv(args.escape_variant_counts)

    ref_name = '_'.join(args.reference_variant_counts.split('/')[-1].split('_')[0:-2])
    escape_name = '_'.join(args.escape_variant_counts.split('/')[-1].split('_')[0:-2])

    lib = ref_name.split('_')[-1]
    assert(ref_name.split('_')[-1] == escape_name.split('_')[-1])

    primary_target = args.primary_target
    output_dir = os.path.join(args.output_dir, escape_name+'_over_'+ref_name)

    os.makedirs(output_dir, exist_ok=True)
    print_some_logs(args)

    ncounts_ref = np.sum(ref['count'].to_numpy())
    ncounts_escape = np.sum(escape['count'].to_numpy())
    frac_escape = args.frac_escape
    if frac_escape is None:
        frac_escape = 1.0

    variants, primary_target = build_variant_table(args.table, args.wildtype)

    ref = (ref.merge(escape[['sample','barcode','count']], on = ['barcode'])
              .rename(columns = {'count_x':'reference_count','count_y':'escape_count'}))
    if SCORE == "ratio":
        ref = ref.assign(escape_score = lambda x: (x.escape_count/ncounts_escape) / (x.reference_count/ncounts_ref) * frac_escape)
    # elif SCORE == "binom":
    #     _ref_ratio = ref['reference_count'] / ncounts_ref
    #     _scores = []
    #     for i in range(ref.shape[0]):
    #         cnt = ref['escape_count'][i]
    #         p = binom_test(cnt, n = ncounts_escape, p=_ref_ratio[i])
    #         if p < 1e-20:
    #             p = 20
    #         else:
    #             p = np.log10(p)
    #         _scores.append(p)
    #     ref = ref.assign(escape_score = _scores)
    # elif SCORE == "poisson":
    #     _ref_ratio = ref['reference_count'] / ncounts_ref
    #     _scores = []
    #     for i in range(ref.shape[0]):
    #         cnt = ref['escape_count'][i]
    #         p = 1-poisson.cdf(k=cnt, mu=ncounts_escape*_ref_ratio[i])
    #         if p < 1e-20:
    #             p = 20
    #         else:
    #             p = np.log10(p)
    #         _scores.append(p)
    #     ref = ref.assign(escape_score = _scores)
    else:
        raise NotImplementedError
    ref = ref[(ref['escape_score'] != np.inf) & (ref['escape_score'] >= 0) & (ref['reference_count'] > 5)].sort_values('escape_score')
    
    # num_use = int(len(ref)*0.99)
    # ref = ref[0:num_use]
    df = ref.assign(aa_sub_type = lambda x: (['>1' if i >1 else str(i) for i in x.n_aa_substitutions]))
    df = df[df['target'] == primary_target].fillna('')
    df = (df.assign(norm_score = (df['escape_score'])/np.std(df['escape_score']), condition = escape_name+'_over_'+ref_name)
            .drop(columns = ['sample_x','sample_y'])
            .pipe(variants.classifyVariants,
                        primary_target=primary_target,
                        syn_as_wt=False))

    mn = np.nanquantile(df['escape_score'], 0.01)
    mx = np.nanquantile(df['escape_score'], 0.99)

    df = df.assign(raw_escape_score = df['escape_score'])
    df['escape_score'] = [(i if (i > mn and i < mx) else (0 if (i < (mn+mx)/2) else mx))/mx for i in df['raw_escape_score']]

    p1 = ggplot(df, aes(x = 'variant_class', y = 'escape_score')) + geom_boxplot()
    bunch.add_plot(p1, 0, 0)

    
    # For RD libraries don't apply this filter
    df = filter_expr_bind(df, args.mutbindexpr, args.variant_bind, args.variant_expr, args.exprmin, args.exprmin, args.bindmin, args.bindmin).query('variant_class != "stop"')

    sys.stdout.write('Fitting epistasis model ...\n')
    binary_map = dms_variants.binarymap.BinaryMap(df, func_score_col='escape_score')
    model = dms_variants.globalepistasis.MonotonicSplineEpistasisGaussianLikelihood(binary_map)
    model.fit()
    sys.stdout.write('Done.\n')

    variant_counter = collections.Counter(model.binarymap.substitution_variants)
    muts_counter = collections.Counter(itertools.chain.from_iterable(s.split() for s in model.binarymap.substitution_variants))
    
    effects_df = (pd.DataFrame({'condition': escape_name+'_over_'+ref_name,
                    'library': lib,
                    'mutation': model.binarymap.all_subs})
        .assign(n_single_mut_measurements=lambda x: x['mutation'].map(variant_counter),
                n_any_mut_measurements=lambda x: x['mutation'].map(muts_counter),
                )
        .pipe(model.add_phenotypes_to_df, substitutions_col='mutation')
        .drop(columns='latent_phenotype')
        .rename(columns={'observed_phenotype': 'epistasis_model_score'})
        )
    sys.stdout.write("Removing mutations that do not have EITHER >= 1 single-mutant measurements or >= 3 any-mutant measurements.\n")
    effects_df = (
        effects_df
        .assign(sufficient_measurements=lambda x: (
                    (x['n_single_mut_measurements'] >= 1) |
                    (x['n_any_mut_measurements'] >= 3))
                )
        .query('sufficient_measurements == True')
        .drop(columns='sufficient_measurements')
    )
    raw_avg_single_mut_scores = (
        df.query('n_aa_substitutions == 1')
        .rename(columns={'aa_substitutions': 'mutation'})
        .groupby(['condition', 'library', 'mutation'])
        .aggregate(raw_single_mut_score=pd.NamedAgg('escape_score', 'mean'))
        .reset_index()
        )
    
    mn = np.nanquantile(effects_df['epistasis_model_score'], 0.01)
    mx = np.nanquantile(effects_df['epistasis_model_score'], 0.99)
    effects_df = effects_df.merge(raw_avg_single_mut_scores, how='outer', validate='one_to_one').assign(site=lambda x: x['mutation'].str[1: -1].astype(int),
            wildtype=lambda x: x['mutation'].str[0],
            mutation=lambda x: x['mutation'].str[-1],
            ).assign(
                mut_escape=lambda x: [(i if (i >= mn and i <= mx) else (mn if (i < (mn+mx)/2) else mx)) for i in x['epistasis_model_score']],
                protein_chain=protein_chain
            )

    mn = np.nanquantile(effects_df['raw_single_mut_score'], 0.01)
    mx = np.nanquantile(effects_df['raw_single_mut_score'], 0.99)

    effects_df['site'] = effects_df['site'] + args.wtseqstart
    effects_df = effects_df.assign(
        single_mut_escape=lambda x: [(i if (i >= mn and i <= mx) else (mn if (i < (mn+mx)/2) else mx)) for i in x['raw_single_mut_score']],
        protein_site=effects_df['site'],
        label_site=effects_df['wildtype']+effects_df['site'].astype(str)
    )

    site_effects_df = (
        effects_df
        .groupby('site')
        .aggregate(
            site_avg_escape_frac_epistasis_model=pd.NamedAgg('mut_escape',
                                                            lambda s: site_average(s, args.site_average, args.pos_or_neg, plot_max)),
            site_total_escape_frac_epistasis_model=pd.NamedAgg('mut_escape', 'sum'),
            site_avg_escape_frac_single_mut=pd.NamedAgg('single_mut_escape',
                                                        lambda s: site_average(s,args.site_average, args.pos_or_neg, plot_max)),
            site_total_escape_frac_single_mut=pd.NamedAgg('single_mut_escape', 'sum'),
            )
        .reset_index()
        )
    site_effects_df.index = site_effects_df['site']

    effects_df = effects_df.assign(
            site_total_escape=site_effects_df['site_total_escape_frac_epistasis_model'][effects_df['site']].to_numpy(),
            site_mean_escape=site_effects_df['site_avg_escape_frac_epistasis_model'][effects_df['site']].to_numpy(),
        )

    effects_df.to_csv(os.path.join(output_dir, "effects_df.csv"), index=False)
    site_effects_df.to_csv(os.path.join(output_dir, "site_effects_df.csv"), index=False)
    site_effects_df = site_effects_df.assign(targetsite = site_effects_df['site'])
    
    plots = []

    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_epistasis_model')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ggtitle(args.escape_variant_counts))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_single_mut')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_total_escape_frac_epistasis_model')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ggtitle(args.escape_variant_counts))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_total_escape_frac_single_mut')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_epistasis_model')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ggtitle(args.escape_variant_counts) + ylim(0,plot_max))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_single_mut')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ylim(0,plot_max))
    
    for i in range(len(plots)):
        plt = plots[i]
        bunch.add_plot(plt, 0, (i+1)*400,1500,300)

    ggsave(bunch, 'call_escape_score_plots.html', path = output_dir)

if __name__ == '__main__':
    main()
