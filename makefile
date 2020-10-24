DEST_DIR = ~/bin

CFLAGS =  -O3 -Wall -Wextra -std=c++11

# ALL =   kmer_positional_binning.o kmer_filtering.o kmer_pair_generator.o kmer_pair_grouping.o kmer_block_constructor.o kmer_block_pair_constructor.o  kmer_block_pair_feature_grouping.o kmer_block_pair_feature_grouped_evaluation.o
# ALL =   kmer_positional_binning.o kmer_filtering.o kmer_pair_generator.o kmer_pair_grouping.o kmer_block_constructor.o kmer_ranged_block_pair_constructor.o  kmer_ranged_block_pair_aggregator_and_filtering.o kmer_ranged_block_pair_updating.o
ALL =   kmer_positional_binning.o kmer_filtering.o kmer_pair_generator.o kmer_pair_grouping.o kmer_block_constructor.o kmer_block_and_neighbours_locator.o kmer_k_neigbour_subgraph_generator.o kmer_metablock_constructor.o \
		kmer_feature_filtering_and_pairing.o kmer_feature_aggregation_and_filtering.o kmer_distant_feature_pair_updating.o kmer_feature_fasta_generator.o kmer_variant_coords_fasta_display.o remove_files.o

all: $(ALL)

kmer_positional_binning.o:
		g++ $(CFLAGS) -o kmer_positional_binning.o kmer_positional_binning.cpp

kmer_filtering.o:
		g++ $(CFLAGS) -o kmer_filtering.o kmer_filtering.cpp

kmer_pair_generator.o:
		g++ $(CFLAGS) -o kmer_pair_generator.o kmer_pair_generator.cpp

kmer_pair_grouping.o:
		g++ $(CFLAGS) -o kmer_pair_grouping.o kmer_pair_grouping.cpp

kmer_block_constructor.o:
		g++ $(CFLAGS) -o kmer_block_constructor.o kmer_block_constructor.cpp

kmer_block_and_neighbours_locator.o:
		g++ $(CFLAGS) -o kmer_block_and_neighbours_locator.o kmer_block_and_neighbours_locator.cpp

kmer_k_neigbour_subgraph_generator.o:
		g++ $(CFLAGS) -o kmer_k_neigbour_subgraph_generator.o kmer_k_neigbour_subgraph_generator.cpp

kmer_metablock_constructor.o:
		g++ $(CFLAGS) -o kmer_metablock_constructor.o kmer_metablock_constructor.cpp

kmer_feature_filtering_and_pairing.o:
		g++ $(CFLAGS) -o kmer_feature_filtering_and_pairing.o kmer_feature_filtering_and_pairing.cpp

kmer_distant_feature_pair_updating.o:
		g++ $(CFLAGS) -o kmer_distant_feature_pair_updating.o kmer_distant_feature_pair_updating.cpp

kmer_feature_aggregation_and_filtering.o:
		g++ $(CFLAGS) -o kmer_feature_aggregation_and_filtering.o kmer_feature_aggregation_and_filtering.cpp

kmer_feature_fasta_generator.o:
		g++ $(CFLAGS) -o kmer_feature_fasta_generator.o kmer_feature_fasta_generator.cpp

kmer_variant_coords_fasta_display.o:
		g++ $(CFLAGS) -o kmer_variant_coords_fasta_display.o kmer_variant_coords_fasta_display.cpp

remove_files.o:
		g++ $(CFLAGS) -o remove_files.o remove_files.cpp

# kmer_block_pair_constructor.o:
# 		g++ $(CFLAGS) -o kmer_block_pair_constructor.o kmer_block_pair_constructor.cpp

# kmer_block_pair_feature_grouping.o:
# 		g++ $(CFLAGS) -o kmer_block_pair_feature_grouping.o kmer_block_pair_feature_grouping.cpp

# kmer_block_pair_feature_grouped_evaluation.o:
# 		g++ $(CFLAGS) -o kmer_block_pair_feature_grouped_evaluation.o kmer_block_pair_feature_grouped_evaluation.cpp

# kmer_ranged_block_pair_constructor.o:
# 		g++ $(CFLAGS) -o kmer_ranged_block_pair_constructor.o kmer_ranged_block_pair_constructor.cpp

# kmer_ranged_block_pair_aggregator_and_filtering.o:
# 		g++ $(CFLAGS) -o kmer_ranged_block_pair_aggregator_and_filtering.o kmer_ranged_block_pair_aggregator_and_filtering.cpp

# kmer_ranged_block_pair_updating.o:
# 		g++ $(CFLAGS) -o kmer_ranged_block_pair_updating.o kmer_ranged_block_pair_updating.cpp

clean:
		rm -f $(ALL)

install:
		cp $(ALL) $(DEST_DIR)