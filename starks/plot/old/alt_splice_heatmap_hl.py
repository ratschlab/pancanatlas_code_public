import os
import sys
import numpy as np
import pandas as pd
import alt_splice_heatmap as ash
sns = ash.sns

BASEDIR = os.path.dirname(os.path.dirname(__file__))
sys.path.append(BASEDIR)
import compute.alt_splice as preproc
import config
import utils

DEBUG = False
if DEBUG: print "Debug mode active!"
MAX_EVENTS = 5000

RESET_DF_CACHE = False
RESET_LINKAGE_CACHE = False

ash.RESET_DF_CACHE = RESET_DF_CACHE
ash.RESET_LINKAGE_CACHE = RESET_LINKAGE_CACHE
CACHE_DIR = ash.CACHE_DIR

METHOD = 'ward'
METRIC = 'cosine'

PLOT_BASE = os.path.join(config.plot_dir, 'altsplice')
if not os.path.exists(PLOT_BASE): os.makedirs(PLOT_BASE)


def filter_df(df, nkeep):
    data , columns = ash.filter_to_high_var(df.values, df.columns, nkeep)
    df = pd.DataFrame(data, index=df.index, columns=columns)
    assert df.shape[1] == nkeep
    return df

def load_full_dataset(map_etype_to_file, max_events, desc, reset):
    cache_path = os.path.join(CACHE_DIR, '%s_df.tsv'%desc)
    if reset or not os.path.exists(cache_path):
        print "Calculating psi values from raw data"
        combined_df = ash.load_high_var_events_combined_df(
                        map_etype_to_file,
                        max_events=max_events,
                        condense=False)
        combined_df.to_csv(cache_path, sep='\t')
    else:
        print "Reading psi from %s" %cache_path
        combined_df = pd.read_csv(cache_path, sep='\t', index_col=0)
    return combined_df

def get_row_colors(metadata_df, subtype, cnc):
    # Add the subtype information
    subtype_series = metadata_df[subtype].loc[metadata_df['cnc'] == cnc]
    full_index = subtype_series.index
    subtype_series = subtype_series.dropna()
    if subtype_series.size == 0: return None, None
    subtype_series = subtype_series.loc[subtype_series.map(lambda x: not 'normal' in x.lower())]
    values = np.sort(np.append(subtype_series.unique(), 'Normal'))
    colors = sns.color_palette("Set1", n_colors=values.size, desat=.5).as_hex()
    color_lut = dict(zip(values, colors))
    row_colors = subtype_series.map(lambda x:color_lut[x])
    full_row_colors = pd.Series(color_lut['Normal'], index=full_index)
    full_row_colors.loc[row_colors.index] = row_colors.values
    return full_row_colors, color_lut

if __name__ == '__main__':
    # event_type to path of raw psi values
    map_etype_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                         'intron_retention': config.alt_splice_intron_retention_path,
                         'alt_3prime': config.alt_splce_alt_3prime_path,
                         'alt_5prime': config.alt_splce_alt_5prime_path}
    etype_desc = 'combined'
    col_cmap = np.array(sns.color_palette('Set2', len(map_etype_to_file.keys())))
    col_color_lut = dict(zip(map_etype_to_file.keys(), col_cmap))

    # Load the data
    combined_df_desc = '%s_%d_high_var_events_concat' %(etype_desc, MAX_EVENTS)
    combined_df = load_full_dataset(map_etype_to_file, MAX_EVENTS, combined_df_desc, RESET_DF_CACHE)
    metadata_df = utils.load_metadata_df(config.metadata_path, combined_df.index)
    metadata_df, subtype_names = utils.append_subtype(config.subtype_path, metadata_df)

    event_types = map_etype_to_file.keys() + ['concatenated']
    for event in event_types:
        if event != 'concatenated':
            event_mask = combined_df.columns.map(lambda x: x.startswith(event)).values
            event_df = combined_df.loc[:, event_mask]
        else:
            event_df = combined_df

        cancer_types = metadata_df['cnc'].unique()
        if DEBUG:
            cancer_types = ['BRCA']
            print "DEBUG: cancer_types = %s" %str(cancer_types)

        cnc_groups = metadata_df.groupby('cnc')
        for cnc in cancer_types:
            filtering_desc = '%s_%s_%d_high_var_events' %(cnc, etype_desc, MAX_EVENTS)
            cnc_index = cnc_groups.get_group(cnc).index
            cnc_data_df = filter_df(event_df.loc[cnc_index], MAX_EVENTS)

            linkage_desc = filtering_desc + '_' + METHOD + '_' + METRIC
            row_linkage, col_linkage = ash.get_linkage(cnc_data_df, linkage_desc,
                                                       method=METHOD, metric=METRIC,
                                                       reset=RESET_LINKAGE_CACHE)

            # coloring tasks
            for subtype in subtype_names:
                row_colors, row_color_lut = get_row_colors(metadata_df, subtype, cnc)
                if row_colors is None: continue
                col_colors = ash.map_col_colors(cnc_data_df.columns, col_color_lut)

                # Plot
                graph = sns.clustermap(cnc_data_df,
                                   row_colors=row_colors, col_colors=col_colors,
                                   row_linkage=row_linkage, col_linkage=col_linkage,
                                   xticklabels=False, yticklabels=False,
                                   linewidths=0)
                title = subtype.replace('Subtype_', '') + '_' + linkage_desc
                graph.ax_col_dendrogram.set_title("AltSplice %s Clustering" %title.replace('_', ' ').title())
                graph.ax_heatmap.set_xlabel("Events")
                graph.ax_heatmap.set_ylabel("Samples")
                graph.cax.set_title("psi")
                ash.add_legend(graph, row_color_lut)
                ash.add_col_legend(graph, col_color_lut)

                outpath = os.path.join(PLOT_BASE, 'concatenated', 'heatmaps', 'cancer_types', cnc, '%s_heatmap.png'%title.strip().lower().replace(' ', '_'))
                if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))
                print "Saving heatmap to: %s" %outpath
                sns.plt.savefig(outpath, bbox_inches='tight')


