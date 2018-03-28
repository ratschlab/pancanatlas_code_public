import os
import sys
import re
import argparse

import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plot_args

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_utils

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import config
import utils

DEBUG = False
# input:
#    embed_dir: contains pca and tsne embeddings
#    base_plot_dir: where to write plots
#    description (e.g. altsplice-exon_skip, expression)

# tasks, fxn of embedding
#    color_cnc: color all cancer types
#    highlight_cnc: color only one cancer type (one plot per cancer)
#    tumor_normal: color tumor and normals
#    library_size: color samples by library size 
#    degscore: color samples by degragradtion score

# pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir)

# task: color_cnc, highlight_cnc

def run_embedding_proc(outdir, plot_configs, desc, tsne_embeds_dict, prefix=None):
    '''Create a scatter plot for each embedding in tsne_embeds_dict
    '''
    plotting_df, legend_kwargs_list, default_kwargs = plot_configs
    if plotting_df is None: return
    for (pp, lr), embed_df in tsne_embeds_dict.items():
        if prefix is None:
            basename = 'tsne_embeds_pp_%d_lr_%d.png' %(pp, lr)
        else:
            basename = '%s_tsne_embeds_pp_%d_lr_%d.png' %(prefix, pp, lr)
        outpath = os.path.join(outdir, basename)
        title = '%s tSNE(perplexity=%d, learning_rate=%d)' %(desc.title(), pp, lr)
        fig, ax = plt.subplots()
        plot_utils.plot_embeddings(embed_df, plotting_df,
                             default_kwargs=default_kwargs,
                             legend_kwargs_list=legend_kwargs_list,
                             ax=ax, title=title)
        plot_utils.save(outpath, do_pdf=DO_PDF)
        plt.close()
    return

def main(embed_dir, plot_dir, desc):
    '''Runs all tasks on a single embedding.

    embed_dir: location of pca & tsne embeddings
    plot_dir: directory to write plots
    desc: identifies embedding (used for e.g. plot titles)
    '''
    assert os.path.exists(embed_dir), embed_dir
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)

    rnadeg_df = utils.load_rnadeg(config.rnadeg_path) if DEGSCORE else None
    libsize_df = utils.load_libsize(config.libsize_path) if LIBSIZE else None

    pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir,
                                                                whitelist=WHITELIST,
                                                                pp_set=TSNE_PP_PLOT_SET,
                                                                lr_set=TSNE_LR_PLOT_SET)

    metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
    metadata_df, subtype_names = utils.append_subtype(config.subtype_path, metadata_df)

    if ALL_TASKS or COLOR_CNC:
        outdir = os.path.join(plot_dir, 'complete')
        if not os.path.exists(outdir): os.makedirs(outdir)
        plot_configs = plot_args.cnc_plotting_args(metadata_df)
        run_embedding_proc(outdir, plot_configs, desc, tsne_embeds_dict)

    if ALL_TASKS or HL_CNC is not None:
        if ALL_TASKS or HL_CNC == 'all':
            cnc_list = np.unique(metadata_df['cnc'].values)
        else:
            cnc_list = HL_CNC
        cnc_groups = metadata_df.groupby('cnc')
        for cnc in cnc_list:
            cnc_index = cnc_groups.get_group(cnc).index
            other_index = np.setdiff1d(metadata_df.index, cnc_index)
            for subtype in subtype_names:
                if DEBUG: print cnc, subtype
                subtype_configs = plot_args.load_subtype_color_tumor_marker_kwargs(metadata_df, subtype, cnc)
                outdir = os.path.join(plot_dir, 'highlights_subtype', cnc)
                if not os.path.exists(outdir): os.makedirs(outdir)
                subtype_desc = ' '.join([subtype, desc])
                run_embedding_proc(outdir, subtype_configs, subtype_desc, tsne_embeds_dict, prefix=subtype)
        pass

    if ALL_TASKS or TUMOR_NORMAL:
        if DEBUG: print "start tumor normal"
        tn_configs = plot_args.tumor_normal_plotting_args(metadata_df)
        outdir = os.path.join(plot_dir, 'complete_tumor_normal')
        if not os.path.exists(outdir): os.makedirs(outdir)
        tn_desc = '%s Tumor/Normal'%desc
        run_embedding_proc(outdir, tn_configs, tn_desc, tsne_embeds_dict)

    if ALL_TASKS or LIBSIZE:
        if DEBUG: print " start libsize"
        cbar_title = libsize_df.columns[0]
        outdir = os.path.join(plot_dir, 'qc', 'libsize')
        if not os.path.exists(outdir): os.makedirs(outdir)
        for (pp, lr), embed_df in tsne_embeds_dict.items():
            outpath = os.path.join(outdir, 'tsne_embeds_pp_%d_lr_%d.png' %(pp, lr))
            axis_title = '%s Library Size Effects\ntSNE(perplexity=%d, learning_rate=%d)' %(desc.title(), pp, lr)
            fig, ax = plt.subplots()
            plot_utils.plot_continuous_color_embeddings(embed_df, libsize_df.iloc[:, 0],
                                                        ax=ax, axis_title=axis_title,
                                                        cbar_title=cbar_title)
            plot_utils.save(outpath, do_pdf=DO_PDF)

    if ALL_TASKS or DEGSCORE:
        if DEBUG: print " start degscore"
        outdir = os.path.join(plot_dir, 'qc', 'degscore')
        if not os.path.exists(outdir): os.makedirs(outdir)
        cbar_title = rnadeg_df.columns[0]
        for (pp, lr), embed_df in tsne_embeds_dict.items():
            outpath = os.path.join(outdir, 'tsne_embeds_pp_%d_lr_%d.png' %(pp, lr))
            fig, ax = plt.subplots()
            axis_title = '%s RNADeg Effects\ntSNE(perplexity=%d, learning_rate=%d)' %(desc.title(), pp, lr)
            plot_utils.plot_continuous_color_embeddings(embed_df, rnadeg_df.iloc[:, 0],
                                                        ax=ax, axis_title=axis_title,
                                                        cbar_title=cbar_title)
            plot_utils.save(outpath, do_pdf=DO_PDF)

DEBUG = False

if __name__ == '__main__':

    WHITELIST = utils.load_whitelist(config.whitelist_path)
    TSNE_PP_PLOT_SET = [10, 50, 200, 500]
    TSNE_LR_PLOT_SET = [500]
    DO_PDF = True

    # Tasks
    ALL_TASKS = True

    COLOR_CNC = True
    HL_CNC = 'all'
    TUMOR_NORMAL = True
    LIBSIZE = True
    DEGSCORE = True

    #Inputs
    embed_dir_list = list()
    plot_dir_list = list()
    desc_list = list()

    # Add Expression
    if True:
        embed_dir_list.append(os.path.join(config.embed_dir, 'expression'))
        plot_dir_list.append(os.path.join(config.plot_dir, 'expression', 'embeds'))
        desc_list.append('Expression')

    # Add AltSplice
    if False:
        altsplice_event_list= ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']
        for event in altsplice_event_list:
            plot_dir_list.append(os.path.join(config.plot_dir, 'altsplice', event, 'embeds'))
            embed_dir_list.append(os.path.join(config.embed_dir, 'altsplice', event))
            desc_list.append('AltSplice %s'%event.title())

    for i, (embed_dir, plot_dir, desc) in enumerate(zip(embed_dir_list, plot_dir_list, desc_list)):
        assert os.path.exists(embed_dir), embed_dir
        print("[%d] %s Embeds: %s"%(i, desc, embed_dir))
        print("[%d] %s Plots:  %s"%(i, desc, plot_dir))

    for i, (embed_dir, plot_dir, desc) in enumerate(zip(embed_dir_list, plot_dir_list, desc_list)):
        main(embed_dir, plot_dir, desc)

