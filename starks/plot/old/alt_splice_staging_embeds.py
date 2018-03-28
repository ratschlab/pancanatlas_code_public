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


if __name__ == '__main__':
    do_pca = True
    do_tsne = True
    tsne_pp_plot_set = [50]
    tsne_lr_plot_set = [500]

    plot_base_dir = os.path.join(config.plot_dir, 'altsplice')
    embed_base_dir = os.path.join(config.embed_dir, 'altsplice')
    event_list = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'concatenated']
    for event in event_list:
        embed_dir = os.path.join(embed_base_dir, event)
        plot_dir = os.path.join(plot_base_dir, event, 'embeds', 'highlight_staging')
        if not os.path.exists(plot_dir): os.makedirs(plot_dir)

        # Load embeds, metadata
        pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir)
        metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
        staging_df = utils.load_staging(config.staging_info, metadata_df.index)
        subtype_names = staging_df.columns
        metadata_df = metadata_df.loc[staging_df.index]
        metadata_df = pd.concat((metadata_df, staging_df), axis=1)

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
                    tsne_plot_out = os.path.join(plot_dir, cnc, '%s_tsne_embeds_pp_%d_lr_%d.png' %(subtype, pp, lr))
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
                    tsne_plot_out_pdf = re.sub('.png$', '.pdf', tsne_plot_out)
                    print "Writing %s" %tsne_plot_out_pdf
                    plt.savefig(tsne_plot_out_pdf, bbox_inches='tight', dpi=300)
                    plt.close()



