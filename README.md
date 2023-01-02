# PRAWNS: Pan-genome representation of whole genomes tool

PRAWNS is a fast and scalable tool that generates an efficient representation of closely related whole genomes to provide a concise list of genomic features or sequence entities shared by a user-specified fraction of the genomes. It is designed specifically for large-scale genomic studies like population genomics and comparative genomics with an emphasis on getting insights into the microbial genome biology. [(Link to publication)](https://doi.org/10.1093/bioinformatics/btac844)

<ins>PRAWNS relies on two main algorithmic innovations:</ins>
1. Locating the <b>*conserved regions*</b> shared across multiple genomes.\
*Metablocks* may contain inexact matches, and reduce the number of features to be considered by an order of magnitude over individual exact-matching *blocks*.
2. We introduce a new type of genomic feature called <b>*paired regions*</b>---these are pairs of *conserved regions* collocated in multiple genomes but that may vary in distance between each other within different genomes.\
Using such paired features, scientists can assess the influence of the collocation of genomic segments on phenotypes.

PRAWNS can be parallelized over multiple threads and uses disk-based storage, enabling it to scale to thousands of genomes.\
The input genomes could be draft asemblies or complete (circularized) genomes; these could represent closely related isolates, including bacterial species or species complex, fungal genomes, viruses, or even mobile elements like plasmid families.

PRAWNS fills up a gap in the current bioinformatics toolkit available to scientists for a holistic large-scale whole-genome population genomics and comparative genomics.\
PRAWNS provides a mechanism for large-scale association studies, including assessing the association of bacterial genomic features with phenotypes, such as antibiotic resistance.

<!-- Given a collection of whole genomes for the closely related isolates (e.g. from a bacterial species or species complex), the tool generates an efficient pan-genome representation for the corresponding isolates. The input genomes could be draft assemblies or complete genomes. The generated pan-genome provides a concise list of structural variants identified across these isolates, which can then be used in a downstream analysis. In addition to the detection of structural variants, the pan-genome also locates the variants that are collocated across multiple isolates—such paired occurrences and the separation between the corresponding structural variants are also provided which can aide the downstream analysis. -->

## Dependencies
1. Python 3.5 or later
2. C++11 or later
3. ```requirements.txt``` enlists the python package versions on which the script was tested (we expect the script to be supported by the later versions as well)

## Installation:
```bash
git clone https://github.com/KiranJavkar/PRAWNS.git
cd PRAWNS
make
```
Run the following command to ensure that all required Python libraries are installed:
```
pip install -r requirements.txt
```

## Running the command:
```
python run_prawns.py input.csv
```
Where ```input.csv``` comprises of upto 3 columns: (i) sample name (assembly/genome) (ii) fasta file path (iii) contig orientations file path (optional, needed if ```--use_oriented_links True```)

```
-bash-4.2$ python run_prawns.py -h
usage: PRAWNS [-h] -i INPUT [-n [NCORES]] [-K [KMER_LEN]] [-o [OUTDIR]]
              [-p [MIN_PERC]] [-l [USE_ORIENTED_LINKS]]
              [-b [MIN_GROUP_BLOCKS]] [-M [MAX_METABLOCK_MISMATCH]]
              [-N [MAX_NEIGHBOR_DISTANCE]] [-s [MIN_BLOCK_SIZE]]
              [-S [MAX_INTERVARIANT_SEPARATION]] [-m [MEM]] [-g [GENOME_LEN]]
              [-V]

PRAWNS: Pan-genome representation of whole genomes tool

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input csv file
  -n [NCORES], --ncores [NCORES]
                        Number of cores to be used (default: 8)
  -K [KMER_LEN], --kmer_len [KMER_LEN]
                        Length of kmers (default: 25)
  -o [OUTDIR], --outdir [OUTDIR]
                        Output directory
  -p [MIN_PERC], --min_perc [MIN_PERC]
                        Minimum % of genomes a variant would be present in
                        (default: 5.0)
  -l [USE_ORIENTED_LINKS], --use_oriented_links [USE_ORIENTED_LINKS]
                        Use MetaCarvel oriented links; if True, 3rd column in
                        input csv should be path to oriented links of
                        corresponding assembly (default: False)
  -b [MIN_GROUP_BLOCKS], --min_group_blocks [MIN_GROUP_BLOCKS]
                        Minimum number of exact matching regions (blocks) that
                        can be grouped into metablocks across the genomes
                        (default: 3)
  -M [MAX_METABLOCK_MISMATCH], --max_metablock_mismatch [MAX_METABLOCK_MISMATCH]
                        Maximum number of mismatches permitted to allow merger
                        and extension of metablocks across the genomes
                        (default: 25)
  -N [MAX_NEIGHBOR_DISTANCE], --max_neighbor_distance [MAX_NEIGHBOR_DISTANCE]
                        Maximum separation between neighboring blocks for
                        collocated blocks (components) identification
                        (default: 5)
  -s [MIN_BLOCK_SIZE], --min_block_size [MIN_BLOCK_SIZE]
                        Smallest size of a block retained as a conserved
                        region (default: 50)
  -S [MAX_INTERVARIANT_SEPARATION], --max_intervariant_separation [MAX_INTERVARIANT_SEPARATION]
                        Maximum number of bases between adjacent conserved
                        regions from a genome to get paired regions (default:
                        50)
  -m [MEM], --mem [MEM]
                        Upper limit for RAM memory usage. Can be in
                        mb/MB/gb/GB/tb/TB (case insensitive), default unit is
                        MB. (default: 36000MB)
  -g [GENOME_LEN], --genome_len [GENOME_LEN]
                        Average genome length. Can be in k/K/m/M/g/G (case
                        insensitive), default unit is M, i.e. 1x10^6 nt.
                        (default: 4M)
  -V, --version         show program's version number and exit
```

For the contig orientations files, we rely on the format used by [MetaCarvel](https://github.com/marbl/MetaCarvel) to represent the oriented links file.\
The oriented links can be generated by running MetaCarvel with `--keep True`. We recommend using `--bsize` argument to have higher confidence for the oriented links obtained for the sequenced isolate under consideration.

## Test dataset run
A test dataset of 10 *Acinetobacter baumannii* genomes has been provided in the `test/` directory. The `test/assemblies/` folder provides the corresponding 10 assemblies in FASTA format. The `test/oriented_links/` folder contains the oriented links obtained after contig scaffolding.\
The `test/input_prawns_test.csv` file forms the comma-separated input file to be given to PRAWNS. This file contains three columns, as described above; if the `--use_oriented_links` is not to be used, then the third column is not necessary.
#### Running PRAWNS with oriented links
```
python -u run_prawns.py --input test/input_prawns_test.csv --outdir test_results  --min_perc 30.0 --use_oriented_links True --ncores 2 -g 4M
```
The above command will identify the conserved and paired regions present in at least 30% (3 out of the 10 genomes) while using the contig scaffolding information. The results would be stored in `test_results/` output directory. The `--ncores` option can be set to the number of CPUs available to run PRAWNS. With 2 CPUs, it would take about 4 minutes for these 10 genomes of approximately 4 mega base pair (`-g 4M`) in length.

#### Running PRAWNS without oriented links
```
python -u run_prawns.py --input test/input_prawns_test.csv --outdir test_results  --min_perc 30.0 --ncores 2 -g 4M
```
This command would execute, similar to the previous one, but without assessing the contig orientations (irrespective of whether the third column is provided in the input file). Compared to the previous command, this command would not check for the presence of *composite metablocks* and paired regions across contigs (described below). The command can also be executed identically by setting `--use_oriented_links False`.

## Description of the output files:
### *Conserved regions*
The *conserved regions* are represented in the form of *metablocks* and *retained blocks*. *Metablocks* are the *conserved regions* obtained by aggregating the *collocated blocks* that have a consistent ordering over multiple genomes (refer PRAWNS manuscript for detailed description of *metablocks*). *Retained blocks* are those that are at least `min_block_size` long and were not aggregated into a metablock (in some genome).

The conserved regions are reported via two comma-separated files:
1. `*_coords.csv` contains the coordinates and strand (forward: 1, reverse: 0) on which the conserved region is present in the corresponding genome.\
If the conserved region is absent in a genome, it would be represented by the tuple 0,0,0.
3. `*_presence_absence.csv` denotes the binary presence (1) or absence (0) for the respective conserved region for the corresponding genome.

Additionally, the `metablocks.fasta` and `retained_blocks.fasta` files provide the conserved region sequences (as in the genome where the respective conserved regions were first encounted) in FASTA format.

#### *Metablocks*
The *metablock* files---`metablock_coords.csv` and `metablock_presence_absence.csv`---have their first column representing the *metablock index*, and the second column denoting the length of the *metablock*. The subsequent columns denote the coordinates-strand tuples or presence/absence respectively.\
The *metablock index* is composed of three parts: *(component id)\_(new metablock id)\_(split number)*. The `metablock_coords.csv` file will denote the entire *metablock index*, whereas `metablock_presence_absence.csv` not use the *split number*.

The *split number* is 0 by default, and would be non-zero only in the case of *composite metablocks*---this refers to a situation where a metablock is spanned across multiple contigs, and requires contig orientation information obtained via scaffolding.\
For instance, say a metablock spanned across 2 contigs in a genome (draft assembly) and is on just one contig in other genomes (whichever other genome contains it). In such case, the two constituent *metablocks* (separated by contig boundaries) would be reported; to denote fragmentation of a longer composite *metbalock*, the *split number* for one constituent *metablock* would be 0 while the other one will have 1. E.g. `3_2_0` and `3_2_1`.\
Observe that the binary presence-absence for constituent *metablocks* is the same and, hence, the *split number* is not required to report the presence/absence of the conserved regions

The information about the *component id* and *new metablock id* is not important for downstream analysis, as they are implicit to PRAWNS' internal computations. However, it is worth noting that all *metablocks* with the same *component id* are arising from the same *component of collocated blocks* and would, hence, be located in the proximity one another in their corresponding genomes.

#### *Retained blocks*
Similar to *metablocks*, the *block* files---`retained_block_coords.csv` and `retained_block_presence_absence.csv`---have the *block index* in the first column and its length in the second column. The subsequent columns denote the coordinates-strand tuples or presence/absence respectively.\
Unlike the *metablocks*, the *block index* is a single number assigned by PRAWNS internally.

If a *block* has been merged into a *metablock* and the corresponding *metablock* has been identified in a genome, the corresponding *retained block* would be deemed absent in the respective genome.\
Therefore, the entries would only denote the genomes in which the *metablocks* were not detected in their entirity but some of the constituent longer *blocks* were present, and the longer *blocks* which were present isolated (didn't have multiple collocated *blocks* shared over several genomes) and could not be aggregated into a *metablock*.


### *Paired regions*
The *paired regions* denote the conserved regions collocated in the given genomes. They provide the means to assess the impact of the genomic context between conserved regions in the context of some phenotype of interest.

The paired regions are presented using two comma-separated files:
1. `intrapair_separation.csv` provides the separation (number of nucleotides) between the collocated conserved regions.\
A negative separation suggests an overlap between the adjoining conserved regions (number of overlapping nucleotides or the nucleotides shared between the adjoining conserved regions) in the corresponding genome/s. If the paired region is missing in a genomes, the associated value for the separation would be 0.
2. `pair_presence_absence.csv` contains the binary presence (1) or absence (0) of the paired regions in the respective genomes.

In both the paired regions files, the first column denotes the *pair index*, and the subsequent columns provide the separation or presence/absence.\
The *pair index* follows the format: *(conserved region 1 index)*|*(conserved region 2 index)*|*(conserved region 1 strand)*|*(conserved region 2 strand)*\
Note that either of the two conserved regions could be *metablocks* or *retained blocks*; the type of conserved region forming the pair can be inferred from the respective *conserved region index* format. Also observe that the same pair of conserved regions but with different relative orientations would result in distint paired regions with distint *paired region indices*.

<!-- + --TO BE UPDATED--
- metablock_coords
- metablock_presence_absence
- similarly for block
- pair_presence_absence
- intrapair_separation -->

## Additional utility output

The `sv_coverages.txt` file provides the coverage of each genome using the *conserved regions*.\
It is a comma-separated file with a header and three columns; the columns correspond to sample name, total genome coverage using the conserved regions, and the number of conserved regions detected within the genome, respectively.

**Extracting a specific conserved region sequence using the coordinates from a particular genome**\
PRAWNS appends all contigs from a multi-fasta file and creates a single coordinate system for each genome. The conserved regions are located using the same coordinate system. The coordinates marking the contig ends can be found in the `contig_lengths` output directory; the file names correspond to the sample number (0-indexed) based on the sample sequence provide in the `input.csv` file.

To get the exact fasta sequence of a variant from a particular genome, run the following command:
`./kmer_variant_coords_fasta_display.o <PRAWNS_results_dir>/all_assembly_filepaths <sample_number(0-indexed)> <start_coordinate(from the corresponding conserved region's coords csv file)> <end_coordinate>`

E.g.: `./kmer_variant_coords_fasta_display.o PRAWNS_results/all_assembly_filepaths.txt 0 2818695 2818729` \
`./kmer_variant_coords_fasta_display.o test_results/all_assembly_filepaths.txt 4 3892906 3893510`


### Distance matrix for population structure adjustment

The *conserved regions* can be used to estimate a distance matrix for the given genomes which can then be used for population structure adjustments in subsequent analyses.

An auxiliary script, `get_distance_matrix_and_newick_tree.py`, has been provided for the distance matrix computation. For *N* genomes, the distance matrix will be a *N*x*N* matrix, where the genome indices (ordering) is same as that provided in the input csv file used for running PRAWNS. This script can be run as: \
`python get_distance_matrix_and_newick_tree.py -i <prawns_input.csv> -d <PRAWNS_results_dir>`

Additionally, this script generates a phylogenetic tree plot, a newick tree for the same, and clusters the genomes into groups of genomes denoting the underlying population structure likely to be exhibited by the given genomes. These groups are represented by distinct group labels (which are numbers starting with 0) and all genomes within the same group will have the same group label.
These group labels can been seen in the phylogenetic tree plot and the newick tree.

## Citation
If you use PRAWNS for your work, please cite the manuscript published in *Bioinformatics*:
```
@article{Javkar2022,
    author = {Javkar, Kiran and Rand, Hugh and Strain, Errol and Pop, Mihai},
    title = "{PRAWNS: Compact pan-genomic features for whole-genome population genomics}",
    journal = {Bioinformatics},
    year = {2022},
    month = {12},
    abstract = "{Scientists seeking to understand the genomic basis of bacterial phenotypes, such as antibiotic resistance, today have access to an unprecedented number of complete and nearly-complete genomes. Making sense of these data requires computational tools able to perform multiple-genome comparisons efficiently, yet currently available tools cannot scale beyond several tens of genomes.We describe PRAWNS, an efficient and scalable tool for multiple-genome analysis. PRAWNS defines a concise set of genomic features (metablocks), as well as pairwise relationships between them, which can be used as a basis for large-scale genotype-phenotype association studies. We demonstrate the effectiveness of PRAWNS by identifying genomic regions associated with antibiotic resistance in Acinetobacter baumannii.PRAWNS is implemented in C ++ and Python3, licensed under the GPLv3 license, and freely downloadable from GitHub (https://github.com/KiranJavkar/PRAWNS.git)Supplementary data are available at Bioinformatics online.}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btac844},
    url = {https://doi.org/10.1093/bioinformatics/btac844},
    note = {btac844},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btac844/48443118/btac844.pdf},
}
```

NOTE: This tool is still under active development and may produce errors while running. Please report any error encountered as a github issue so that we can fix it during the development. For any questions, please email kjavkar[AT]umd[DOT]edu or kiran.javkar2707[AT]gmail[DOT]com.
