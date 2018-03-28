import os
import sys
import re
import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import embedding_plotter as eplt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

DO_TSNE = True

DO_INDIVIDUAL_EVENTS = True
DO_COMBINED_EVENTS = True

TSNE_LR_LIST = 'all'
TSNE_PP_LIST = 'all'

PLOT_DIR = os.path.join(config.plot_dir, 'altsplice')
EMBED_DIR = os.path.join(config.embed_dir, 'altsplice')
EVENT_LIST = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'concatenated']

TSNE_PP_PLOT_SET  = 'all'
TSNE_LR_PLOT_SET  = 'all'

def load_cancers():
    metadata_df = utils.load_metadata_df(config.metadata_path)
    if DEBUG: return ['COAD']
    return np.sort(metadata_df['cnc'].unique())

def save_as_pdf(png_path):
    assert png_path.endswith('.png')
    pdf_out = re.sub('.png$', '.pdf', png_path)
    print "Writing %s" %pdf_out
    plt.savefig(pdf_out, bbox_inches='tight')
    return

def _assert_correct_index(metadata_df, tsne_embeds, cnc):
    index_with_meta = tsne_embeds.index.intersection(metadata_df.index)
    if index_with_meta.size < tsne_embeds.index.size:
        print "WARNING: %f size reduction in index" %(float(index_with_meta.size) / tsne_embeds.index.size)
    assert index_with_meta.size > 0
    assert metadata_df['cnc'].loc[index_with_meta].unique().size == 1
    assert metadata_df['cnc'].loc[index_with_meta].unique()[0] == cnc
    return index_with_meta


def run_subtype_embedding_procs():
    cancers = load_cancers()
    metadata_df = None
    for event in ['alt_3prime']:
        for cnc in cancers:
            cnc_event_embed_dir = os.path.join(EMBED_DIR, event, 'cancer_only', cnc)
            assert os.path.exists(cnc_event_embed_dir)
            pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(cnc_event_embed_dir)
            if metadata_df is None or not metadata_df.index.equals(pca_embeds.index):
                metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
                metadata_df, subtype_names = utils.append_subtype(config.subtype_path, metadata_df)
            for subtype in subtype_names:
                plotting_df, legend_kwargs, default_kwargs = \
                    eplt.load_subtype_color_tumor_marker_kwargs(metadata_df, subtype, cnc)
                if plotting_df is None: continue
                default_kwargs.update({'alpha': .6, 's':20})
                for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                    if not TSNE_PP_PLOT_SET == 'all' and not pp in TSNE_PP_PLOT_SET: continue
                    if not TSNE_LR_PLOT_SET == 'all' and not lr in TSNE_LR_PLOT_SET: continue
                    index_with_meta = _assert_correct_index(metadata_df, tsne_embeds, cnc)

                    tsne_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_subtype_cnc_subset', cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
                    if SKIP_EXISTING_PLOTS and os.path.exists(tsne_plot_out): continue
                    if not os.path.exists(os.path.dirname(tsne_plot_out)): os.makedirs(os.path.dirname(tsne_plot_out))

                    embeds = tsne_embeds.loc[index_with_meta]
                    fig, ax = eplt.plot_tsne_embeds(embeds, plotting_df, **default_kwargs)
                    eplt.add_legend(ax, legend_kwargs, loc='bottom', ncol=5)
                    fig.suptitle('AltSplice %s %s %s\ntSNE(perplexity=%d, learning_rate=%d)' %(subtype.title(), cnc, event, pp, lr))
                    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                    print "Writing %s" %tsne_plot_out
                    plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                    save_as_pdf(tsne_plot_out)
                    plt.close()

def run_staging_embedding_procs():
    cancers = load_cancers()
    metadata_df = None
    for event in EVENT_LIST:
        for cnc in cancers:
            cnc_event_embed_dir = os.path.join(EMBED_DIR, event, 'cancer_only', cnc)
            assert os.path.exists(cnc_event_embed_dir)
            pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(cnc_event_embed_dir)
            if metadata_df is None or not metadata_df.index.equals(pca_embeds.index):
                metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
                staging_df = utils.load_staging(config.staging_info, metadata_df.index)
                subtype_names = staging_df.columns
                metadata_df = metadata_df.loc[staging_df.index]
                metadata_df = pd.concat((metadata_df, staging_df), axis=1)

            for subtype in subtype_names:
                plotting_df, legend_kwargs, default_kwargs = \
                    eplt.load_subtype_color_tumor_marker_kwargs(metadata_df, subtype, cnc)
                if plotting_df is None: continue
                default_kwargs.update({'alpha': .6, 's':20})
                for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                    if not TSNE_PP_PLOT_SET == 'all' and not pp in TSNE_PP_PLOT_SET: continue
                    if not TSNE_LR_PLOT_SET == 'all' and not lr in TSNE_LR_PLOT_SET: continue
                    index_with_meta = _assert_correct_index(metadata_df, tsne_embeds, cnc)

                    tsne_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_subtype_cnc_subset', cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
                    if SKIP_EXISTING_PLOTS and os.path.exists(tsne_plot_out): continue
                    if not os.path.exists(os.path.dirname(tsne_plot_out)): os.makedirs(os.path.dirname(tsne_plot_out))

                    embeds = tsne_embeds.loc[index_with_meta]
                    fig, ax = eplt.plot_tsne_embeds(embeds, plotting_df, **default_kwargs)
                    eplt.add_legend(ax, legend_kwargs, loc='bottom', ncol=5)
                    fig.suptitle('AltSplice %s %s %s\ntSNE(perplexity=%d, learning_rate=%d)' %(subtype.title(), cnc, event, pp, lr))
                    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                    print "Writing %s" %tsne_plot_out
                    plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                    save_as_pdf(tsne_plot_out)
                    plt.close()
                    if DEBUG: sys.exit()

def run_smoking_procs():
    cancers = load_cancers()
    metadata_df = None
    for event in EVENT_LIST:
        for cnc in cancers:
            cnc_event_embed_dir = os.path.join(EMBED_DIR, event, 'cancer_only', cnc)
            assert os.path.exists(cnc_event_embed_dir)
            pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(cnc_event_embed_dir)
            if metadata_df is None or not metadata_df.index.equals(pca_embeds.index):
                metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
                metadata_df, subtype_names = utils.append_smoking_and_age(config.smoking_age_path, metadata_df)
            metadata_df = metadata_df.loc[metadata_df['is_tumor']]
            assert 'tobacco_smoking_history' in metadata_df.columns
            assert 'age_began_smoking_in_years' in metadata_df.columns
            subtype = 'tobacco_smoking_history'
            plotting_df, legend_kwargs, default_kwargs = \
                eplt.load_subtype_color_tumor_marker_kwargs(metadata_df, subtype, cnc)
            if plotting_df is None: continue
            default_kwargs.update({'alpha': .6, 's':20})

            pca_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_smoking_and_age_cnc_subset', cnc, '%s_pca_embeds.png' %(subtype))
            if SKIP_EXISTING_PLOTS and os.path.exists(pca_plot_out):
                pass
            else:
                index_with_meta = _assert_correct_index(metadata_df, pca_embeds, cnc)
                embeds = pca_embeds.loc[index_with_meta]
                plt.close()

                fig, axes = eplt.plot_pca_embeds(pca_embeds, plotting_df, **default_kwargs)
                eplt.add_legend(axes[-1][1], legend_kwargs, loc='bottom', ncol=5)
                eplt.add_pca_variance_expl(axes, pca_model)
                fig.suptitle('AltSplice %s PCA'%event.title())
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])

                plt.suptitle('PCA AltSplice %s %s-Only %s' %(subtype.title(), cnc, event))

                if not os.path.exists(os.path.dirname(pca_plot_out)): os.makedirs(os.path.dirname(pca_plot_out))
                print "Writing %s" %pca_plot_out
                plt.savefig(pca_plot_out, bbox_inches='tight', dpi=300)
                save_as_pdf(pca_plot_out)


            for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                if not TSNE_PP_PLOT_SET == 'all' and not pp in TSNE_PP_PLOT_SET: continue
                if not TSNE_LR_PLOT_SET == 'all' and not lr in TSNE_LR_PLOT_SET: continue
                index_with_meta = _assert_correct_index(metadata_df, tsne_embeds, cnc)

                tsne_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_smoking_and_age_cnc_subset', cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
                if SKIP_EXISTING_PLOTS and os.path.exists(tsne_plot_out): continue
                if not os.path.exists(os.path.dirname(tsne_plot_out)): os.makedirs(os.path.dirname(tsne_plot_out))

                embeds = tsne_embeds.loc[index_with_meta]
                fig, ax = eplt.plot_tsne_embeds(embeds, plotting_df, **default_kwargs)
                eplt.add_legend(ax, legend_kwargs, loc='bottom', ncol=5)
                fig.suptitle('AltSplice %s %s %s\ntSNE(perplexity=%d, learning_rate=%d)' %(subtype.title(), cnc, event, pp, lr))
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                print "Writing %s" %tsne_plot_out
                plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                save_as_pdf(tsne_plot_out)
                plt.close()
    return

def plot_age_embeddings(embeds, age_series):
    assert embeds.index.equals(age_series.index.intersection(embeds.index))
    null_index = age_series.index[age_series.isnull()]
    age_index = age_series.index.difference(null_index)
    age_values = age_series.loc[age_index].astype(int)
    null_index = null_index.append(age_index[age_values > 100])
    age_index = age_index[age_values <= 100]
    age_values = age_series.loc[age_index]

    if age_index.size == 0: return None
    if null_index.size > 0:
        # handle missing age
        pass

#    cmap = sns.dark_palette("#5d7ad4", as_cmap=True)
    cmap = sns.cubehelix_palette(as_cmap=True)
    age_values = age_series.loc[age_index].values.astype(int)
    x,y = embeds.columns[:2]
    fig, ax = plt.subplots()
    cax = ax.scatter(embeds.loc[age_index,x], embeds.loc[age_index,y], c=age_values, cmap=cmap)
    ax.set_xticks([])
    ax.set_yticks([])
    cbar = fig.colorbar(cax)
    return True

def run_age_procs():
    cancers = load_cancers()
    metadata_df = None
    for event in EVENT_LIST:
        for cnc in cancers:
            cnc_event_embed_dir = os.path.join(EMBED_DIR, event, 'cancer_only', cnc)
            assert os.path.exists(cnc_event_embed_dir)
            pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(cnc_event_embed_dir)
            metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
            metadata_df, subtype_names = utils.append_smoking_and_age(config.smoking_age_path, metadata_df)
            assert 'tobacco_smoking_history' in metadata_df.columns
            assert 'age_began_smoking_in_years' in metadata_df.columns
            subtype = 'age_began_smoking_in_years'

            pca_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_smoking_and_age_cnc_subset', cnc, '%s_pca_embeds.png' %(subtype))
            if SKIP_EXISTING_PLOTS and os.path.exists(pca_plot_out):
                pass
            else:
                index_with_meta = _assert_correct_index(metadata_df, pca_embeds, cnc)
                embeds = pca_embeds.loc[index_with_meta]
                plt.close()
                success = plot_age_embeddings(embeds, metadata_df[subtype])
                if success is None: continue
                plt.suptitle('PCA AltSplice %s %s-Only %s' %(subtype.title(), cnc, event))

                if not os.path.exists(os.path.dirname(pca_plot_out)): os.makedirs(os.path.dirname(pca_plot_out))
                print "Writing %s" %pca_plot_out
                plt.savefig(pca_plot_out, bbox_inches='tight', dpi=300)
                save_as_pdf(pca_plot_out)

            for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                if not TSNE_PP_PLOT_SET == 'all' and not pp in TSNE_PP_PLOT_SET: continue
                if not TSNE_LR_PLOT_SET == 'all' and not lr in TSNE_LR_PLOT_SET: continue

                tsne_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_smoking_and_age_cnc_subset', cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
                if SKIP_EXISTING_PLOTS and os.path.exists(tsne_plot_out):
                    pass
                else:
                    index_with_meta = _assert_correct_index(metadata_df, tsne_embeds, cnc)
                    embeds = tsne_embeds.loc[index_with_meta]
                    plt.close()
                    success = plot_age_embeddings(embeds, metadata_df[subtype])
                    if success is None: continue
                    plt.suptitle('AltSplice %s %s %s\ntSNE(perplexity=%d, learning_rate=%d)' %(subtype.title(), cnc, event, pp, lr))

                    tsne_plot_out = os.path.join(PLOT_DIR, event, 'embeds', 'highlights_smoking_and_age_cnc_subset', cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
                    if not os.path.exists(os.path.dirname(tsne_plot_out)): os.makedirs(os.path.dirname(tsne_plot_out))
                    print "Writing %s" %tsne_plot_out
                    plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                    save_as_pdf(tsne_plot_out)
                    plt.close()
    return

DEBUG = False
SKIP_EXISTING_PLOTS = False

if __name__ == '__main__':
    do_subtypes = True
    do_staging  = True
    do_smoking = False
    do_age = False

    if do_subtypes: run_subtype_embedding_procs()
    if do_staging: run_staging_embedding_procs()
    if do_smoking: run_smoking_procs()
    if do_age: run_age_procs()

