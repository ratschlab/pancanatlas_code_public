import os
import sys
import re
import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import embedding_plotter as eplt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

def save_as_pdf(png_path):
    assert png_path.endswith('.png')
    pdf_out = re.sub('.png$', '.pdf', png_path)
    print "Writing %s" %pdf_out
    plt.savefig(pdf_out, bbox_inches='tight')
    return

def plot_pca_embeddings(pca_embeds, pca_model, plotting_df, legend_kwargs, default_kwargs):
    fig, axes = eplt.plot_pca_embeds(pca_embeds, plotting_df, **default_kwargs)
    eplt.add_legend(axes[-1][1], legend_kwargs, loc='bottom', ncol=5)
    eplt.add_pca_variance_expl(axes, pca_model)
    fig.suptitle('AltSplice %s PCA'%event.title())
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axes

def plot_tsne_embeddings(tsne_embeds, plotting_df, legend_kwargs, default_kwargs):
    fig, ax = eplt.plot_tsne_embeds(tsne_embeds, plotting_df, **default_kwargs)
    eplt.add_legend(ax, legend_kwargs, loc='bottom', ncol=5)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax

def run_tumor_normal_coloring(tn_plot_dir, event, metadata_df, pca_model, pca_embeds, tsne_embeds_dict):
#    cmap = dict(zip([True, False], eplt.sns.color_palette("Set2", 2)))
    cmap = dict(zip([True, False], ['lightblue', 'orange']))
    default_kwargs = {'alpha': .25, 's':10}
    is_tumor_to_marker = {True:'o', False:'^'}

    plotting_df = pd.DataFrame(index=metadata_df.index)
    plotting_df['facecolors'] = metadata_df['is_tumor'].map(cmap)
    plotting_df['marker'] = metadata_df['is_tumor'].map(is_tumor_to_marker)

    # List of kwargs used to create a legend
    legend_kwargs_list = list()
    legend_defaults = {'lw':0, 's': 30, 'alpha':.75}
    legend_kwargs_list.append(
        eplt._copy_and_update(
            legend_defaults,
            {'c':cmap[True],
             'label': 'Tumor',
             'marker': is_tumor_to_marker[True]
            }))

    legend_kwargs_list.append(
        eplt._copy_and_update(
            legend_defaults,
            {'c':cmap[False],
             'label': 'Normal',
             'marker': is_tumor_to_marker[False]
            }))

    if not os.path.exists(tn_plot_dir): os.makedirs(tn_plot_dir)
    if do_pca:
        pca_plot_out = os.path.join(tn_plot_dir, 'pca_embeds.png')
        tumor_idx = metadata_df.index[metadata_df['is_tumor']]
        normal_idx = metadata_df.index[~metadata_df['is_tumor']]
        fig, axes = eplt.plot_pca_embeds(pca_embeds.loc[tumor_idx], plotting_df.loc[tumor_idx], **default_kwargs)
        fig, axes = eplt.plot_pca_embeds(pca_embeds.loc[normal_idx], plotting_df.loc[normal_idx], axes=axes, **default_kwargs)
        eplt.add_legend(axes[-1][1], legend_kwargs_list, loc='bottom', ncol=5)
        eplt.add_pca_variance_expl(axes, pca_model)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.suptitle('AltSplice %s PCA'%event.title())
        print "Writing %s" %pca_plot_out
        plt.savefig(pca_plot_out, bbox_inches='tight', dpi=300)
        if do_pdf: save_as_pdf(pca_plot_out)
        plt.close()

    if do_tsne:
        for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
            if not tsne_pp_plot_set == 'all' and not pp in tsne_pp_plot_set: continue
            if not tsne_lr_plot_set == 'all' and not lr in tsne_lr_plot_set: continue
            tsne_plot_out = os.path.join(tn_plot_dir, 'tsne_embeds_pp_%d_lr_%d.png' %(pp, lr))
            fig, ax = eplt.plot_tsne_embeds(tsne_embeds.loc[tumor_idx], plotting_df.loc[tumor_idx], **default_kwargs)
            fig, ax = eplt.plot_tsne_embeds(tsne_embeds.loc[normal_idx], plotting_df.loc[normal_idx], ax=ax, **default_kwargs)
            eplt.add_legend(ax, legend_kwargs_list, loc='bottom', ncol=5)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            fig.suptitle('AltSplice %s tSNE(perplexity=%d, learning_rate=%d)' %(event.title(), pp, lr))
            print "Writing %s" %tsne_plot_out
            plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
            if do_pdf: save_as_pdf(tsne_plot_out)
            plt.close()
    return

DEBUG = True
do_pdf = True
do_pca = True
do_tsne = True
tsne_pp_plot_set = 'all' #[50]
tsne_lr_plot_set = 'all' #[500]

do_regular_embeddings = False
do_hl_embeddings = False
do_tumor_normal = True


# INPUT:
# directory containing embeddings -- pca & tsne

if __name__ == '__main__':
    plot_base_dir = os.path.join(config.plot_dir, 'altsplice_2018')
    embed_base_dir = os.path.join(config.embed_dir, 'altsplice_2018')
    event_list = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'concatenated']
#    event_list = ['concatenated']
    for event in event_list:
        embed_dir = os.path.join(embed_base_dir, event)
        plot_dir = os.path.join(plot_base_dir, event, 'embeds')
        if not os.path.exists(plot_dir): os.makedirs(plot_dir)

        # Load embeds, metadata
        pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir)
        metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
        metadata_df, subtype_names = utils.append_subtype(config.subtype_path, metadata_df)
        if (DEBUG): sys.exit()

        if do_tumor_normal:
            tn_plot_dir = os.path.join(plot_base_dir, event, 'embeds', 'complete_tumor_normal')
            run_tumor_normal_coloring(tn_plot_dir, event, metadata_df, pca_model, pca_embeds, tsne_embeds_dict)

        if do_regular_embeddings:
            plotting_df, legend_kwargs, default_kwargs = \
                eplt.load_cnc_color_tumor_marker_kwargs(metadata_df)
            reg_plot_dir = os.path.join(plot_base_dir, event, 'embeds', 'complete')
            if not os.path.exists(reg_plot_dir): os.makedirs(reg_plot_dir)
            if do_pca:
                pca_plot_out = os.path.join(reg_plot_dir, 'pca_embeds.png')
                fig, axes = plot_pca_embeddings(pca_embeds, pca_model, plotting_df, legend_kwargs, default_kwargs)
                fig.suptitle('AltSplice %s PCA'%event.title())
                print "Writing %s" %pca_plot_out
                plt.savefig(pca_plot_out, bbox_inches='tight', dpi=300)
                if do_pdf: save_as_pdf(pca_plot_out)
                plt.close()

            if do_tsne:
                for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                    if not tsne_pp_plot_set == 'all' and not pp in tsne_pp_plot_set: continue
                    if not tsne_lr_plot_set == 'all' and not lr in tsne_lr_plot_set: continue
                    tsne_plot_out = os.path.join(reg_plot_dir, 'tsne_embeds_pp_%d_lr_%d.png' %(pp, lr))
                    fig, ax = plot_tsne_embeddings(tsne_embeds, plotting_df, legend_kwargs, default_kwargs)
                    fig.suptitle('AltSplice %s tSNE(perplexity=%d, learning_rate=%d)' %(event.title(), pp, lr))
                    print "Writing %s" %tsne_plot_out
                    plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                    if do_pdf: save_as_pdf(tsne_plot_out)
                    plt.close()

        if do_hl_embeddings:
            cnc_groups = metadata_df.groupby('cnc')
            for cnc in np.unique(metadata_df['cnc'].values):
                cnc_index = cnc_groups.get_group(cnc).index
                other_index = np.setdiff1d(metadata_df.index, cnc_index)
                for subtype in subtype_names:
                    plotting_df, legend_kwargs, default_kwargs = \
                        eplt.load_subtype_color_tumor_marker_kwargs(metadata_df, subtype, cnc)
                    if plotting_df is None: continue
                    default_kwargs.update({'alpha': .6, 's':20})
                    # call plot
                    for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                        if not tsne_pp_plot_set == 'all' and not pp in tsne_pp_plot_set: continue
                        if not tsne_lr_plot_set == 'all' and not lr in tsne_lr_plot_set: continue
                        tsne_plot_out = os.path.join(plot_dir, 'highlights_subtype', cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
                        if not os.path.exists(os.path.dirname(tsne_plot_out)): os.makedirs(os.path.dirname(tsne_plot_out))
                        embeds = tsne_embeds_dict[(pp, lr)]
                        bg_defaults = default_kwargs.copy(); bg_defaults['alpha'] = .4
                        fig, ax = eplt.plot_tsne_embeds(embeds.loc[other_index], plotting_df.loc[other_index], ax=None, **bg_defaults)
                        fig, ax = eplt.plot_tsne_embeds(embeds.loc[cnc_index], plotting_df.loc[cnc_index], ax=ax, **default_kwargs)
                        eplt.add_legend(ax, legend_kwargs, loc='bottom', ncol=5)
                        fig.suptitle('AltSplice %s %s %s\ntSNE(perplexity=%d, learning_rate=%d)' %(subtype.title(), cnc, event, pp, lr))
                        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                        print "Writing %s" %tsne_plot_out
                        plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                        if do_pdf: save_as_pdf(tsne_plot_out)
                        plt.close()
