rm(list=ls())
library(tailfinder)
fast5_dir <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/few_4000_per_fast5_reads'
data='pcr-dna'
save_dir <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop'
num_cores <- 120
save_plots=F
show_plots=F
plot_debug=F
multifast5=T
basecalled_with_flipflop=T
plotting_library='rbokeh'

a <- find_dna_tails_parallel(fast5_dir,
                        data=data,
                        save_dir=save_dir,
                        csv_filename='dna_tails.csv',
                        num_cores=num_cores,
                        save_plots=save_plots,
                        show_plots=show_plots,
                        plot_debug=plot_debug,
                        multifast5=multifast5,
                        basecalled_with_flipflop=basecalled_with_flipflop,
                        plotting_library=plotting_library)


