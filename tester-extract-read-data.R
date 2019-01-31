rm(list=ls())
library(tailfinder)
guppy.single.flipflop.dna <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/data_for_readtype_exploration/1_guppy-single-flipflop-dna.fast5'
guppy.multi.flipflop.dna <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/data_for_readtype_exploration/2_guppy-multi-flipflop-dna.fast5'
albacore.single.standard.rna <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/data_for_readtype_exploration/3_albacore-single-standard-rna.fast5'
albacore.multi.standard.rna <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/data_for_readtype_exploration/6_albacore-multi-standard-rna.fast5'
guppy.single.standard.rna <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/data_for_readtype_exploration/4_guppy-single-standard-rna.fast5'
guppy.multi.standard.rna <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/data_for_readtype_exploration/5_guppy-multi-standard-rna.fast5'



extract_read_data(read_path = albacore.single.standard.rna,
                  read_id_fast5_file = NA,
                  plot_debug = F,
                  basecalled_with = 'albacore',
                  multifast5 = F,
                  model = 'standard')


extract_read_data(read_path = NA,
                  read_id_fast5_file = list(read_id='read_001a68d8-f83e-4a76-a6ea-a1ad3d001b75', fast5_file=albacore.multi.standard.rna),
                  plot_debug = F,
                  basecalled_with = 'albacore',
                  multifast5 = T,
                  model = 'standard')

extract_read_data(read_path = guppy.single.flipflop.dna,
                  read_id_fast5_file = NA,
                  plot_debug = F,
                  basecalled_with = 'guppy',
                  multifast5 = F,
                  model = 'flipflop')

extract_read_data(read_path = guppy.multi.flipflop.dna,
                  read_id_fast5_file = list(read_id='read_0015e3a5-f3a2-4d2f-b288-c5b712f2c0b1', fast5_file=guppy.multi.flipflop.dna),
                  plot_debug = F,
                  basecalled_with = 'guppy',
                  multifast5 = T,
                  model = 'flipflop')

extract_read_data(read_path = guppy.single.standard.rna,
                  read_id_fast5_file = NA,
                  plot_debug = F,
                  basecalled_with = 'guppy',
                  multifast5 = F,
                  model = 'standard')

extract_read_data(read_path = guppy.multi.standard.rna,
                  read_id_fast5_file = list(read_id='read_000ebe45-cfa9-43d5-9d32-51971c50c4c7', fast5_file=guppy.multi.standard.rna),
                  plot_debug = F,
                  basecalled_with = 'guppy',
                  multifast5 = T,
                  model = 'standard')



