import os
import sys
import pandas as pd
import config
import alt_splice

EMBED_DIR = config.embed_dir

if __name__ == '__main__':

    map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                         'intron_retention': config.alt_splice_intron_retention_path,
                         'alt_3prime': config.alt_splce_alt_3prime_path,
                         'alt_5prime': config.alt_splce_alt_5prime_path}
    if False:
        for event, path in map_event_to_file.items():
            event_embed_dir = os.path.join(EMBED_DIR, 'altsplice_2018', event)
            if not os.path.exists(event_embed_dir): os.makedirs(event_embed_dir)
            psi, strains, gene_idx = alt_splice.load_data(path)
            psi_df = pd.DataFrame(psi, index=strains, columns=gene_idx.astype(int))
            if not os.path.exists(event_embed_dir): os.makedirs(event_embed_dir)
            outpath = os.path.join(event_embed_dir, '%s_filtered_psi.tsv'%event)
            print('Writing to %s' %outpath)
            psi_df.to_csv(outpath, sep='\t')

    combined_psi, strains, gene_idx = alt_splice.load_combined_events(map_event_to_file)
    comb_psi_df = pd.DataFrame(combined_psi, index=strains, columns=range(len(gene_idx)))
    event_embed_dir = os.path.join(EMBED_DIR, 'altsplice_2018', 'concatenated')
    if not os.path.exists(event_embed_dir): os.makedirs(event_embed_dir)
    outpath = os.path.join(event_embed_dir, 'combined_filtered_psi.tsv')
    print('Writing to %s' %outpath)
    psi_df.to_csv(outpath, sep='\t')

