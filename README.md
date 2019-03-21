
<!-- README.md is generated from README.Rmd. Please edit that file -->
tailfindr <a href=''><img src='man/figures/tailfindr-logo.png' align="right" height="250" /></a>
================================================================================================

<!-- badges: start -->
<!-- badges: end -->
What is tailfindr?
------------------

tailfindr is a R package for estimating poly(A)-tail lengths in Oxford Nanopore reads.

Features of tailfindr
---------------------

-   Works for both RNA and DNA reads. In the case of DNA reads, it estimates both poly(A)- and poly(T)-tail lengths.
-   Supports data that has been basecalled with Albacore or Guppy. It also support data that has been basecalled using the newer 'flipflop' model.
-   Can work on single or multi-fast5 file reads.

tailfindr has been developed at [Valen Lab](https://www.cbu.uib.no/valen/) in [Computational Biology Unit](https://www.cbu.uib.no/) at the [University of Bergen](https://www.uib.no/), Norway.

Installation
------------

#### Step 1. Installing HDF5 library

tailfindr depends on the HDF5 library for reading Fast5 files. For OS X and Linux, the HDF5 library needs to be installed via one of the (shell) commands specified below:

| System                                      | Command                            |
|:--------------------------------------------|:-----------------------------------|
| **OS X (using Homebrew)**                   | `brew install hdf5`                |
| **Debian-based systems (including Ubuntu)** | `sudo apt-get install libhdf5-dev` |
| **Systems supporting yum and RPMs**         | `sudo yum install hdf5-devel`      |

HDF5 1.8.14 has been pre-compiled for Windows and is available [here](https://github.com/mannau/h5-libwin) — thus no manual installation is required.

#### Step 2. Installing devtools

Currently, tailfindr is not listed on CRAN/Bioconductor, so you need to install it using `devtools`. To install `devtools` use the following command:

``` r
install.packages("devtools")
```

#### Step 3. Installing tailfindr

If you want to install tailfindr without building the vignette, then run the command below:

``` r
devtools::install_github("adnaniazi/tailfindr")
```

If you also want to build the vignette while installing tailfindr, then run the command below:

``` r
remotes::install_github('adnaniazi/tailfindr', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), force = TRUE)
```

Now you are ready to use tailfindr.

Usage
-----

#### 1. Minimal working example

`find_tails()` is the main function that you can use to find tail lengths in both RNA and DNA reads. It saves a CSV file containing all the tail-length data. Furthermore, it also returns the same data as a tibble.

Give below is a minimal use case in which we will run tailfindr on example RNA reads present in the tailfindr package.

``` r
library(tailfindr)
df <- find_tails(fast5_dir = system.file('extdata', 'rna', package = 'tailfindr'),
                 save_dir = '~/Downloads',
                 csv_filename = 'rna_tails.csv',
                 num_cores = 2)
```

In the above example, tailfindr returns a tibble containing the tail data which is then stored in the variable `df`. tailfindr also savs this dataframe as a csv file (`rna_tails.csv`) in the user-specified `save_dir`, which in this case is set to `~/Downloads`. A logfile is also saved in the `save_dir`. The parameter `num_cores` can be increased depending on the number of *physical* cores at your disposal.

#### 2. Plotting the tail

Additionally, tailfindr allows you to generate plots that show the tail location in the raw squiggle. You can save these plots as interactive `.html` files by using `'rbokeh'` as the `plotting_library`. You can zoom in on the tail region in the squiggle and see the exact location of the tail.

Give below is a minimal use case in which we will run tailfindr on example cDNA reads present in the tailfindr package, and also save the plots:

``` r
df <- find_tails(fast5_dir = fast5_dir <- system.file('extdata', 'cdna', package = 'tailfindr'),
                 save_dir = '~/Downloads',
                 csv_filename = 'cdna_tails.csv',
                 num_cores = 2,
                 save_plots = TRUE,
                 plotting_library = 'rbokeh')
```

![Poly(T) read squiggle plot](https://github.com/adnaniazi/tailfindr/raw/master/man/figures/poly_t_without_debug.gif)

However, note that generating plots can slow down the performace of tailfindr. We recommend that you generate these plots only for a small subset of your reads.

#### 3. Plotting the tail and debug traces

tailfindr can plot additional information that it used while deriving the tail boundaries. Please read our preprint to learn how tailfindr works. To plot this information, set the `plot_debug_traces` parameter to `TRUE`.

``` r
df <- find_tails(fast5_dir = fast5_dir <- system.file('extdata', 'cdna', package = 'tailfindr'),
                 save_dir = '~/Downloads',
                 csv_filename = 'cdna_tails.csv',
                 num_cores = 2,
                 save_plots = TRUE,
                 plot_debug_traces = TRUE,
                 plotting_library = 'rbokeh')
```

![Poly(A) read squiggle plot](https://github.com/adnaniazi/tailfindr/raw/master/man/figures/poly_a_with_debug.gif)

There are more options available in the find\_tails() function. Please see its [documentation](https://rdrr.io/github/adnaniazi/tailfindr/man/find_tails.html).

### Description of the CSV/Dataframe columns

tailfindr returns tail data in a dataframe and also saves this information in a user-specified CSV file. The columns generated depend on the whether tailfindr was run on RNA or DNA data. Below is a description of columns for both thses scenarios:

##### When input data is RNA

<table>
<colgroup>
<col width="12%" />
<col width="8%" />
<col width="78%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Column Names</th>
<th align="left">Datatype</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">read_id</td>
<td align="left">character</td>
<td align="left">Read ID as given in the Fast5 file</td>
</tr>
<tr class="even">
<td align="left">tail_start</td>
<td align="left">numeric</td>
<td align="left">Sample index of start site of the tail in raw data</td>
</tr>
<tr class="odd">
<td align="left">tail_end</td>
<td align="left">numeric</td>
<td align="left">Sample index of end site of the tail in raw data</td>
</tr>
<tr class="even">
<td align="left">samples_per_nt</td>
<td align="left">numeric</td>
<td align="left">Read rate in terms of samples per nucleotide</td>
</tr>
<tr class="odd">
<td align="left">tail_length</td>
<td align="left">numeric</td>
<td align="left">Tail length in nucleotides. It is the difference between <code>tail_end</code> and <code>tail_start</code> divided by <code>samples_per_nt</code></td>
</tr>
<tr class="even">
<td align="left">file_path</td>
<td align="left">character</td>
<td align="left">Absolute path of the Fast5 file</td>
</tr>
</tbody>
</table>

##### When input data is DNA

Here are the columns that you will get from tailfindr if you have run it on DNA data:

<table>
<colgroup>
<col width="7%" />
<col width="8%" />
<col width="83%" />
</colgroup>
<thead>
<tr class="header">
<th>Column Names</th>
<th>Datatype</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>read_id</td>
<td>character</td>
<td>Read ID as given in the Fast5 file</td>
</tr>
<tr class="even">
<td>read_type</td>
<td>character factor</td>
<td>Whether a read is <code>&quot;polyA&quot;</code>, <code>&quot;polyT&quot;</code>, or <code>&quot;invalid&quot;</code>. Invalid reads are those in which tailfindr wasn't able to find Nanopore primers with high confidence.</td>
</tr>
<tr class="odd">
<td>tail_is_valid</td>
<td>logical</td>
<td>Whether a poly(A) tail is a full-length read or not. This is important because a poly(A) tail is at the end of the read, and premature termination of reads is prevelant in cDNA.</td>
</tr>
<tr class="even">
<td>tail_start</td>
<td>numeric</td>
<td>Sample index of start site of the tail in raw data</td>
</tr>
<tr class="odd">
<td>tail_end</td>
<td>numeric</td>
<td>Sample index of end site of the tail in raw data</td>
</tr>
<tr class="even">
<td>samples_per_nt</td>
<td>numeric</td>
<td>Read rate in terms of samples per nucleotide</td>
</tr>
<tr class="odd">
<td>tail_length</td>
<td>numeric</td>
<td>Tail length in nucleotides. It is the difference between <code>tail_end</code> and <code>tail_start</code> divided by <code>samples_per_nt</code></td>
</tr>
<tr class="even">
<td>file_path</td>
<td>character</td>
<td>Absolute path of the Fast5 file</td>
</tr>
</tbody>
</table>

The devil👹 in the details
-------------------------

-   tailfindr currently works on data in the `/Analyses/Basecall_1D_000/BaseCalled_template/` path of the Fast5 file data hierarchy. It won't work on data present in, lets say, `/Analyses/Basecall_1D_001/BaseCalled_template/` path or `/Analyses/Basecall_1D_002/BaseCalled_template/` path; such paths are generated if you re-basecall already-basecalled data. To avoid this problem, use tailfindr on files that have been basecalled from the raw Fast5 files.
-   For DNA data, tailfindr decides whether a read is poly(A) or poly(T) based on finding Nanopore primers/adaptors. If you are using the flipflop model to basecall DNA data, please ensure that the nanopore adaptors are not trimmed off while basecalling. This can be done by turning off `enabling_trimming` option in the basecalling script. The script below shows you how we have basecalled our reads using the flipflop model

``` bash
#!/bin/sh
INPUT=/raw/fast5/files/path/
OUTPUT=/output/folder/path/
guppy_basecaller \
    --config dna_r9.4.1_450bps_flipflop.cfg \
    --input $INPUT \
    --save_path $OUTPUT \
    --recursive \
    --fast5_out \
    --hp_correct 1 \
    --disable_pings 1 \
    --enable_trimming 0 
```

Getting help
------------

If you encounter a clear bug, please file a minimal reproducible example on [github](https://github.com/adnaniazi/tailfindr/issues). For questions and other discussion, email me at <adnan.niazi@uib.no>.

License
-------

And of course:

GPL-3: <https://www.gnu.org/licenses/gpl-3.0.en.html>
