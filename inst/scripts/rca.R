#load the csv file
rm(list = ls())
library(dplyr)
library(tailfindr)
library(Biostrings)
df <- read.csv("/Users/adnaniazi/Documents/phd/code/tailfindr/inst/scripts/wholeGFP-hit.csv",
               header = TRUE, stringsAsFactors = FALSE)
df$filepath <-  gsub(pattern = "/Users/max/R/RCA/",
                     replacement = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/",
                     x = df$filepath)
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/4/b86a1494-ae62-4092-a104-5df8f64c22a1.fast5"
filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/13/9e10a4a2-f884-4861-ab6e-b0d73667d93a.fast5"
# Extract read data
read_data <- extract_read_data(file_path = filepath,
                               read_id_fast5_file = NA,
                               plot_debug = FALSE,
                               basecalled_with = 'Guppy',
                               basecall_group = 'Basecall_1D_000',
                               multifast5 = FALSE,
                               model = 'flipflop',
                               plotting_library = 'rbokeh')
fastq <- Biostrings::DNAString(read_data$fastq)

# define parameters for pairwise alignment
match <- 1
mismatch <- -1
type <- 'local'
gapOpening <- 0
gapExtension <- 1
submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                   mismatch = mismatch)


# poly(A) read

##############
### 1. GFP ###
##############
gfp <- Biostrings::DNAString('GGATCACTCTCGGCATGGACGAGCTGTACAAGTAG')
rc_gfp <- Biostrings::reverseComplement(gfp)

as_gfp <- Biostrings::pairwiseAlignment(pattern = gfp,
                                        subject = fastq,
                                        substitutionMatrix = submat,
                                        type = type,
                                        #scoreOnly = FALSE,
                                        gapOpening = gapOpening,
                                        gapExtension = gapExtension)
as_rc_gfp <- Biostrings::pairwiseAlignment(pattern = rc_gfp,
                                        subject = fastq,
                                        substitutionMatrix = submat,
                                        type = type,
                                        #scoreOnly = FALSE,
                                        gapOpening = gapOpening,
                                        gapExtension = gapExtension)
as_gfp_mp <- matchPattern(pattern = gfp,
                          subject = fastq,
                          max.mismatch=5,
                          min.mismatch=0,
                          with.indels=TRUE,
                          fixed=TRUE,
                          algorithm="indels")

as_rc_gfp_mp <- matchPattern(pattern = rc_gfp,
                             subject = fastq,
                             max.mismatch=5,
                             min.mismatch=0,
                             with.indels=TRUE,
                             fixed=TRUE,
                             algorithm="indels")

as_gfp
as_rc_gfp
as_gfp_mp
as_rc_gfp_mp

##################
### 2. Oligo T ###
##################
oligoT <- Biostrings::DNAString("GCCTGTCGCTCTATCTTC")
rc_oligoT <- Biostrings::DNAString("GAAGATAGAGCGACAGGC")
as_rc_oligoT <- Biostrings::pairwiseAlignment(pattern = rc_oligoT,
                                              subject = fastq,
                                              substitutionMatrix = submat,
                                              type = type,
                                              #scoreOnly = FALSE,
                                              gapOpening = gapOpening,
                                              gapExtension = gapExtension)

as_oligoT <- Biostrings::pairwiseAlignment(pattern = oligoT,
                                              subject = fastq,
                                              substitutionMatrix = submat,
                                              type = type,
                                              #scoreOnly = FALSE,
                                              gapOpening = gapOpening,
                                              gapExtension = gapExtension)

as_oligoT_mp <- matchPattern(pattern = oligoT,
                          subject = fastq,
                          max.mismatch=2,
                          min.mismatch=0,
                          with.indels=TRUE,
                          fixed=TRUE,
                          algorithm="indels")

as_rc_oligoT_mp <- matchPattern(pattern = rc_oligoT,
                      subject = fastq,
                      max.mismatch=2,
                      min.mismatch=0,
                      with.indels=TRUE,
                      fixed=TRUE,
                      algorithm="indels")
as_oligoT
as_rc_oligoT
as_oligoT_mp
as_rc_oligoT_mp


########################
### 2. Strand switch ###
########################
ss <- Biostrings::DNAString("AACGATGTACGGCGGGGAAGATAGAGCGACAGGC")
rc_ss <- Biostrings::DNAString("GCCTGTCGCTCTATCTTCCCCGCCGTACATCGTT")

as_ss <- Biostrings::pairwiseAlignment(pattern = ss,
                                              subject = fastq,
                                              substitutionMatrix = submat,
                                              type = type,
                                              #scoreOnly = FALSE,
                                              gapOpening = gapOpening,
                                              gapExtension = gapExtension)

as_rc_ss <- Biostrings::pairwiseAlignment(pattern = rc_ss,
                                           subject = fastq,
                                           substitutionMatrix = submat,
                                           type = type,
                                           #scoreOnly = FALSE,
                                           gapOpening = gapOpening,
                                           gapExtension = gapExtension)

as_ss_mp <- matchPattern(pattern = ss,
                             subject = fastq,
                             max.mismatch=8,
                             min.mismatch=0,
                             with.indels=TRUE,
                             fixed=TRUE,
                             algorithm="indels")

as_rc_ss_mp <- matchPattern(pattern = rc_ss,
                                subject = fastq,
                                max.mismatch=8,
                                min.mismatch=0,
                                with.indels=TRUE,
                                fixed=TRUE,
                                algorithm="indels")
as_ss
as_rc_ss
as_ss_mp
as_rc_ss_mp


## Get Event data
event_data <- read_data$event_data
event_data <- event_data %>%
    mutate(reduced_model_state = ifelse(move == 1, model_state, '')) %>%
    dplyr::slice(18000:27000) %>%
    mutate(a_base = ifelse(move == 1 & model_state == 'A', 1, 0)) %>%
    mutate(t_base = ifelse(move == 1 & model_state == 'T', 1, 0)) %>%
    mutate(c_base = ifelse(move == 1 & model_state == 'C', 1, 0)) %>%
    mutate(g_base = ifelse(move == 1 & model_state == 'G', 1, 0)) %>%
    mutate(atcg_base = 0) %>%
    mutate(atcg_base = ifelse(move == 1 & model_state == 'A', 1, atcg_base)) %>%
    mutate(atcg_base = ifelse(move == 1 & model_state == 'T', 2, atcg_base)) %>%
    mutate(atcg_base = ifelse(move == 1 & model_state == 'C', 3, atcg_base)) %>%
    mutate(atcg_base = ifelse(move == 1 & model_state == 'G', 4, atcg_base))

raw_data <- read_data$raw_data
raw_data <- raw_data[event_data$start[1]:event_data$start[nrow(event_data)]]
raw_data <- (raw_data - mean(raw_data))/sd(raw_data)
raw_data <- abs(raw_data)
index <- event_data$start[1]: event_data$start[nrow(event_data)]
data <- data.frame(index, raw_data)

library(rbokeh) # load echarts4r

p1 <- rbokeh::figure(data = data,
                     width = 1500,
                     height = 250) %>%
    rbokeh::ly_lines(x = index,
                     y = raw_data,
                     width = 1.5,
                     color = '#b2b2b2',
                     legend = "Raw data") %>%
    ly_text(x = start, y = 0, text = reduced_model_state, data = event_data,
            font_style = "bold", font_size = "8pt",
            align = "left", baseline = "middle") %>%
    rbokeh::tool_pan(dimensions = "width") %>%
    rbokeh::tool_wheel_zoom(dimensions = "width")
p1

#
# p2 <- rbokeh::figure(data = event_data, width=1500, height=100) %>%
#     rbokeh::ly_crect(x=start,
#                      y = (event_data$a_base)/2,
#                      width=rep(0.01, nrow(event_data)),
#                      height=event_data$a_base,
#                      color='#529c82',
#                      hover = list(start, model_state))
#
# p3 <- rbokeh::figure(data = event_data, width=1500, height=100) %>%
#     rbokeh::ly_crect(x=start,
#                      y = (event_data$t_base)/2,
#                      width=rep(0.01, nrow(event_data)),
#                      height=event_data$t_base,
#                      color='#529c82')
# p4 <- rbokeh::figure(data = event_data, width=1500, height=100) %>%
#     rbokeh::ly_crect(x=start,
#                      y = (event_data$c_base)/2,
#                      width=rep(0.01, nrow(event_data)),
#                      height=event_data$c_base,
#                      color='#529c82')
# p5 <- rbokeh::figure(data = event_data, width=1500, height=100) %>%
#     rbokeh::ly_crect(x=start,
#                      y = (event_data$g_base)/2,
#                      width=rep(0.01, nrow(event_data)),
#                      height=event_data$g_base,
#                      color='#529c82')
#
# p6 <- rbokeh::figure(width = 1500, height = 500) %>%
#     rbokeh::ly_crect(x = start,
#                      y = (event_data$atcg_base)/2,
#                      data = event_data,
#                      width = rep(0.2, nrow(event_data)),
#                      height = event_data$atcg_base,
#                      color = '#529c82',
#                      hover = list(event_data$start)) %>%
#     ly_text(x = start, y = 1, text = model_state, data = event_data,
#             font_style = "bold", font_size = "8pt",
#             align = "left", baseline = "middle")
# p6
#
# lst <- list(p1, p6, p2, p3, p4, p5)
# names(lst) <- c('Raw data', ' ATCG', 'A', 'T', 'C', 'G')
# nrow <- 6
#
#
#
# p <- rbokeh::grid_plot(lst,
#                        nrow = nrow,
#                        link_data = FALSE,
#                        same_axes=c(TRUE, FALSE))
# p
#


# sequences between polyA and polyT stretches
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 18925,
                                               end = 19190)
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 21650,
                                               end = 21925)
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 22875,
                                               end = 23120)

# sequence before polyA
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 20740,
                                               end = 21560)
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 22120,
                                               end = 22800)
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 24490,
                                               end = 25335)
tailfindr::extract_sequence_between_boundaries(event_data = event_data,
                                               start = 25775,
                                               end = 26640)
