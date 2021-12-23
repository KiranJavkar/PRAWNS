#include <iostream>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <cassert>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <set>
#include <utility>
#include <math.h>

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using namespace std;

typedef pair<string,short> pss;
typedef unsigned long long int ulli;
typedef pair<short,list<ulli> > pslistulli;
typedef pair<ulli,ulli> pulliulli;
typedef pair<short,pulliulli > pspuu;
typedef pair<pulliulli,unsigned long > ppuuul;
typedef pair<unsigned long,ulli> pululli;
typedef pair<ulli,bool> pullib;
typedef pair<unsigned int,bool> puib;
typedef pair<int,bool> pib;
typedef tuple<ulli,int,bool> tup_uib;
typedef tuple<ulli,ulli,bool> tup_uub;
typedef pair<short,ulli> psulli;
typedef tuple<ulli,psulli,bool> tup_upsub;
typedef tuple<string,ulli,ulli,ulli> tup_suuu;
typedef tuple<unsigned int,ulli,ulli,ulli> tup_iuuu;
typedef pair<string, unsigned int> psui;
typedef pair<string, ulli> pstrulli;
typedef tuple<unsigned int,bool,unsigned int,bool,float,float,short> tup_ibibffs;
typedef tuple<psulli,bool,unsigned int> tup_psubi;
typedef tuple<psulli,psulli,bool,bool> tup_psupsubb;
typedef tuple<ulli,tup_psupsubb> tup_ulli_tpsupsubb;
typedef tuple<tup_psubi,tup_psubi> tup_tpsubitpsubi;
typedef pair<psulli,vector<bool> > ppsuvecb;
typedef pair<tup_psupsubb,vector<bool> > ptuppsuvecb;
typedef pair<psulli,list<ulli> > ppsulistulli;
typedef pair<tup_psupsubb,list<ulli> > ptuppsulistulli;
typedef pair<psulli,ulli > ppsuulli;
typedef tuple<ulli,ulli,bool,bool> tup_uubb;
typedef tuple<unsigned long,unsigned long,bool,bool> tup_ululbb;
typedef tuple<ulli,ulli,bool,bool,long> tup_uubbl;
typedef tuple<ulli,ulli,int,bool> tup_uuib;
typedef tuple<ulli,ulli,bool> tup_uub;
typedef tuple<unsigned long,ulli,ulli,short> tup_uluus;
typedef tuple<unsigned long,short,bool> tup_ulsb;
typedef tuple<unsigned long,short,short> tup_ulss;
typedef tuple<unsigned long,short,short,vector<short> > tup_ulssvecs;
typedef tuple<ulli,bool,int> tup_ullibi;
typedef tuple<tup_ullibi,tup_ullibi> tup_tubitubi;
typedef tuple<unsigned long,bool,int> tup_ulbi;
typedef tuple<tup_ulbi,tup_ulbi> tup_tulbitulbi;
typedef tuple<ulli,bool,long> tup_ubl;
typedef tuple<tup_ubl,tup_ubl> tup_tubltubl;
typedef pair<tup_ululbb,vector<bool> > ptupulvecb;
typedef pair<short,pullib> ps_pullib;
typedef pair<ulli,pullib> pullipullib;
typedef pair<unsigned int,unsigned int> puiui;
typedef tuple<bool,bool,float,float,short> tup_bbffs;
typedef tuple<ulli,ulli,bool,ulli> tup_uubu;


ulli strtoulli(string number){
    ulli result = 0;
    for(int pos = 0; pos< number.length(); pos++){
        result *=10;
        result += number[pos]-'0';
    }
    return result;
}


bool is_valid_split_string(string item){
    if(item.length()>1)
        return true;
    if(item.length()==0)
        return false;
    char c = item[0];
    switch(c){
        case ' ': case '\t': case '\n': case '\b': case '\r': case ',': case '.':
        return false;
        default:
        return true;
    }
}


template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        if(is_valid_split_string(item))
            *(result++) = item;
    }
}


list<string> split(const string &s, char delim) {
    list<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}


list<string> get_lines(const char* FILENAME){
    ifstream input_file(FILENAME);
    list<string> lines_list;
    string text;
    while (getline( input_file, text )){
        lines_list.push_back(text);
    }
    return lines_list;
}


list<string> get_filepaths(const char* FILENAME){
    ifstream input_file(FILENAME);
    string filepath;
    list<string> filepaths;
    while (getline( input_file, filepath )){
        filepaths.push_back(filepath);
    }
    return filepaths;
}


// struct Assemblywise_Collinear_Block{
//     // bool coordinate_update_ignore;
//     bool orientation;
//     tup_psupsubb start_kmer_pair_feature_tuple, end_kmer_pair_feature_tuple;
//     ulli block_idx, start_pos, end_pos; // In worst case, number of blocks ==> number of kmer pairs: upper bounded by genome length
// };


struct Collinear_Block_ends{
    ulli block_idx;
    tup_psupsubb start_kmer_pair_feature_tuple, end_kmer_pair_feature_tuple;
    bool ref_kmer_pair_strand; // to be used only in the case of a block with a single kmer pair
};


struct Assemblyspecific_Block{
    // bool coordinate_update_ignore;
    bool orientation;
    ulli block_idx, start_pos, end_pos;
};


void load_all_vector_files( const char* all_vector_blocks_filename, const char* all_vector_start_block_idx_filename,
                            list<string>& all_filename_list, list<ulli>& all_start_block_idx_list){
    ifstream input_file_vector_blocks_files(all_vector_blocks_filename);
    ifstream input_file_start_block_indices(all_vector_start_block_idx_filename);

    string curr_col_bool, curr_block_filename, curr_start_block_idx_str;

    while (getline( input_file_vector_blocks_files, curr_block_filename )){
        getline( input_file_start_block_indices, curr_start_block_idx_str );
        all_filename_list.push_back(curr_block_filename);
        all_start_block_idx_list.push_back(strtoulli(curr_start_block_idx_str));
    }

    input_file_vector_blocks_files.close();
    input_file_start_block_indices.close();
}


/*void load_kmer_pair_map_groups( map< tup_psupsubb, ulli >& kmer_pair_map, string binned_kmer_pairs_dir,
                                string assembly_idx_str, short feature_partitions){
    cout<<"load_kmer_pair_map_groups started: "<<assembly_idx_str<<" "<<kmer_pair_map.size()<<" "<<feature_partitions<<"\n";
    ifstream inFile_kmer_pair;
    string binned_pair_file_prefix = binned_kmer_pairs_dir + assembly_idx_str + "_";
    tup_psupsubb kmer_pair_feature_tuple;
    ulli position, map_feature_count;

    // init_map_group_lists(kmer_pair_map_list, feature_partitions);

    for(short map_group_idx=0; map_group_idx<feature_partitions; map_group_idx++){
        // cout<<"\t"<<map_group_idx<<"\n";

        inFile_kmer_pair.open(binned_pair_file_prefix+to_string(map_group_idx),ios::binary);
        inFile_kmer_pair.read((char*) (&map_feature_count), sizeof(map_feature_count));

        // cout<<"\t"<<map_group_idx<<" "<<map_feature_count<<"\n";

        while(map_feature_count--){
            inFile_kmer_pair.read((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
            inFile_kmer_pair.read((char*) (&position), sizeof(position));
            kmer_pair_map[kmer_pair_feature_tuple] = position;
        }

        inFile_kmer_pair.close();
    }
    cout<<"load_kmer_pair_map_groups: "<<binned_pair_file_prefix<<" "<<kmer_pair_map.size()<<"\n";
}*/


void load_kmer_pair_map_groups( map< tup_psupsubb, pullib >& kmer_pair_map, string binned_kmer_pairs_dir,
                                string assembly_idx_str, short feature_partitions){
    ifstream inFile_kmer_pair;
    string binned_pair_file_prefix = binned_kmer_pairs_dir + assembly_idx_str + "_";
    tup_psupsubb kmer_pair_feature_tuple;
    pullib position_strand_pair;
    ulli map_feature_count;

    // init_map_group_lists(kmer_pair_map_list, feature_partitions);

    for(short map_group_idx=0; map_group_idx<feature_partitions; map_group_idx++){

        inFile_kmer_pair.open(binned_pair_file_prefix+to_string(map_group_idx),ios::binary);
        inFile_kmer_pair.read((char*) (&map_feature_count), sizeof(map_feature_count));

        while(map_feature_count--){
            inFile_kmer_pair.read((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
            inFile_kmer_pair.read((char*) (&position_strand_pair), sizeof(position_strand_pair));
            kmer_pair_map[kmer_pair_feature_tuple] = position_strand_pair;
        }

        inFile_kmer_pair.close();
    }
    cout<<"load_kmer_pair_map_groups: "<<binned_pair_file_prefix<<" "<<kmer_pair_map.size()<<"\n";
}


/*void locate_and_save_assemblyspecific_collinear_blocks(list<Assemblyspecific_Block>& block_list, list<string> all_filename_list,
                                                    list<ulli> all_start_block_idx_list, unsigned long assembly_idx,
                                                    const char* bool_col_filename, map< tup_psupsubb, ulli > kmer_pair_map,
                                                    short kmer_len, ulli total_block_count, string assembly_blocks_out_dir,
                                                    map< string, list<Collinear_Block_ends> >& loaded_blocks_map){

    cout<<"locate_and_save_assemblyspecific_collinear_blocks started: "<<bool_col_filename<<" "<<assembly_idx<<" "<<kmer_pair_map.size()<<"\n";

    string hash_string, curr_col_bool, curr_block_filename;
    ulli curr_start_block_idx, block_count, start_pos, end_pos, pos_1, pos_2, missing_idx;

    list<string>::iterator it_blocks_file = all_filename_list.begin();
    list<ulli>::iterator it_start_idx = all_start_block_idx_list.begin();
    map< string, list<Collinear_Block_ends> >::iterator it_loaded_map;
    list<Collinear_Block_ends> current_blocks_list;
    Collinear_Block_ends current_block;
    Assemblyspecific_Block current_assemblyspecific_block;
    list<Collinear_Block_ends>::iterator it_block_ends;
    list<string> lsplit;

    unsigned long block_assembly_idx;

    ifstream input_file_bool_col(bool_col_filename);
    ifstream inFile;

    string outstring = "";
    ofstream outFile;

    while (getline( input_file_bool_col, curr_col_bool )){
        // cout<<(*it_start_idx)<<" "<<curr_col_bool<<"\n";

        if(curr_col_bool == "1"){
            curr_block_filename = *it_blocks_file;
            lsplit = split(curr_block_filename, '_');
            hash_string = lsplit.back();
            curr_start_block_idx = *it_start_idx;

            // Check if the blocks from current filemap have already been loaded in the map
            // else (1st load of the filemap blocks) read the blocks and save them in the map
            //          if a block is not already loaded, it is possible that the first genome requiring the block to be loaded
            //          is actually the first genome containing the block ==> ends coordinates directly available while reading the block

            it_loaded_map = loaded_blocks_map.find(hash_string); 

            if(it_loaded_map != loaded_blocks_map.end()){
                current_blocks_list = it_loaded_map->second;

                for(it_block_ends = current_blocks_list.begin(); it_block_ends != current_blocks_list.end(); it_block_ends++){
                    current_assemblyspecific_block.block_idx = it_block_ends->block_idx;

                    pos_1 = kmer_pair_map[it_block_ends->start_kmer_pair_feature_tuple];
                    pos_2 = kmer_pair_map[it_block_ends->end_kmer_pair_feature_tuple];

                    current_assemblyspecific_block.orientation = (pos_1<pos_2);

                    if(current_assemblyspecific_block.orientation){
                        current_assemblyspecific_block.start_pos = pos_1;
                        current_assemblyspecific_block.end_pos = pos_2+kmer_len;
                    }
                    else{
                        current_assemblyspecific_block.start_pos = pos_2;
                        current_assemblyspecific_block.end_pos = pos_1+kmer_len;
                    }

                    // if(block_list.size()%2500==0){
                    //     cout<<"\t"<<block_list.size()<<" "<<current_assemblyspecific_block.start_pos<<" ";
                    //     cout<<current_assemblyspecific_block.end_pos<<" "<<current_assemblyspecific_block.orientation<<"\n";
                    // }

                    block_list.push_back(current_assemblyspecific_block);
                }
            }
            else{

                current_blocks_list.clear();
                inFile.open(curr_block_filename,ios::binary);
                inFile.read((char*) (&block_count), sizeof(block_count));
                inFile.read((char*) (&block_assembly_idx), sizeof(block_assembly_idx));

                for(ulli idx=0; idx<block_count; idx++){

                    current_block.block_idx = curr_start_block_idx + idx;
                    current_assemblyspecific_block.block_idx = current_block.block_idx;

                    inFile.read((char*) (&current_block.start_kmer_pair_feature_tuple), sizeof(current_block.start_kmer_pair_feature_tuple));
                    inFile.read((char*) (&current_block.end_kmer_pair_feature_tuple), sizeof(current_block.end_kmer_pair_feature_tuple));

                    current_blocks_list.push_back(current_block);

                    inFile.read((char*) (&start_pos), sizeof(start_pos));
                    inFile.read((char*) (&end_pos), sizeof(end_pos));

                    if(block_assembly_idx == assembly_idx){
                        current_assemblyspecific_block.start_pos = start_pos;
                        current_assemblyspecific_block.end_pos = end_pos;
                        current_assemblyspecific_block.orientation = (start_pos<end_pos);
                    }
                    else{
                        pos_1 = kmer_pair_map[current_block.start_kmer_pair_feature_tuple];
                        pos_2 = kmer_pair_map[current_block.end_kmer_pair_feature_tuple];


                        current_assemblyspecific_block.orientation = (pos_1<pos_2);

                        if(current_assemblyspecific_block.orientation){
                            current_assemblyspecific_block.start_pos = pos_1;
                            current_assemblyspecific_block.end_pos = pos_2+kmer_len;
                        }
                        else{
                            current_assemblyspecific_block.start_pos = pos_2;
                            current_assemblyspecific_block.end_pos = pos_1+kmer_len;
                        }
                    }

                    // if(block_list.size()%2500==0){
                    //     cout<<"\t"<<block_list.size()<<" "<<current_assemblyspecific_block.start_pos<<" ";
                    //     cout<<current_assemblyspecific_block.end_pos<<" "<<current_assemblyspecific_block.orientation<<"\n";
                    // }

                    block_list.push_back(current_assemblyspecific_block);
                    // outstring += to_string(current_assemblyspecific_block.start_pos) + ",";
                    // outstring += to_string(current_assemblyspecific_block.end_pos) + ",";
                    // outstring += (current_assemblyspecific_block.orientation)?"1,":"0,";
                }

                loaded_blocks_map[hash_string] = current_blocks_list;
                inFile.close();
            }
            
        }
        // else{
        //     // cout<<"\t"<<outstring<<"\n\n";
        //     missing_idx = *it_start_idx;
        //     if(next(it_start_idx) != all_start_block_idx_list.end()){
        //         while(missing_idx < *(next(it_start_idx))){
        //             outstring += "0,0,0,";
        //             missing_idx++;
        //         }
        //     }
        //     else{
        //         while(missing_idx<total_block_count){
        //             outstring += "0,0,0,";
        //             missing_idx++;
        //         }
        //     }
        //     // cout<<"\t\t"<<outstring<<"\n\n\n";
        // }

        it_blocks_file++;
        it_start_idx++;
    }

    cout<<"individual blocks located: "<<assembly_idx<<": "<<block_list.size()<<"\n";

    remove(bool_col_filename);

    ulli prev_block_idx = 0, diff;
    list<Assemblyspecific_Block>::iterator it_blocks = block_list.begin();

    while(prev_block_idx < it_blocks->block_idx){
        outstring += "0,0,0,";
        prev_block_idx++;
    }

    while(it_blocks != block_list.end()){
        diff = it_blocks->block_idx - prev_block_idx;
        while(diff>1){
            outstring += "0,0,0,";
            diff--;
        }
        outstring += to_string(it_blocks->start_pos) + ",";
        outstring += to_string(it_blocks->end_pos) + ",";
        outstring += (it_blocks->orientation)?"1,":"0,";

        prev_block_idx = it_blocks->block_idx;
        it_blocks++;
    }

    while(prev_block_idx++ < total_block_count)
        outstring += "0,0,0,";

    outstring.pop_back();
    outstring += "\n";
    outFile.open(assembly_blocks_out_dir + to_string(assembly_idx), ios::out);
    outFile << outstring;
    outFile.close();
    cout<<"locate_and_save_assemblyspecific_collinear_blocks ended: "<<assembly_idx<<": "<<block_list.size()<<"\n";
}*/


void locate_and_save_assemblyspecific_collinear_blocks(list<Assemblyspecific_Block>& block_list, list<string> all_filename_list,
                                                    list<ulli> all_start_block_idx_list, unsigned long assembly_idx,
                                                    const char* bool_col_filename, map< tup_psupsubb, pullib > kmer_pair_map,
                                                    short kmer_len, ulli total_block_count, string assembly_blocks_out_dir,
                                                    map< string, list<Collinear_Block_ends> >& loaded_blocks_map, unsigned long min_sv_size){

    cout<<"locate_and_save_assemblyspecific_collinear_blocks started: "<<bool_col_filename<<" "<<assembly_idx<<" "<<kmer_pair_map.size()<<"\n";

    string hash_string, curr_col_bool, curr_block_filename;
    ulli curr_start_block_idx, block_count, start_pos, end_pos, pos_1, pos_2, missing_idx;
    pullib pos_pair_1, pos_pair_2;
    bool ref_kmer_pair_strand;

    list<string>::iterator it_blocks_file = all_filename_list.begin();
    list<ulli>::iterator it_start_idx = all_start_block_idx_list.begin();
    map< string, list<Collinear_Block_ends> >::iterator it_loaded_map;
    list<Collinear_Block_ends> current_blocks_list;
    Collinear_Block_ends current_block;
    Assemblyspecific_Block current_assemblyspecific_block;
    list<Collinear_Block_ends>::iterator it_block_ends;
    list<string> lsplit;

    list<tup_uubu> first_occurrence_blocks_list;
    list<tup_uubu>::iterator it_first_occurrence_blocks_list;
    tup_uubu block_coords;
    tup_uub coords;

    unsigned long block_assembly_idx;

    ifstream input_file_bool_col(bool_col_filename);
    ifstream inFile;

    string outstring = "";
    ofstream outFile;

    while (getline( input_file_bool_col, curr_col_bool )){
        // cout<<(*it_start_idx)<<" "<<curr_col_bool<<"\n";

        if(curr_col_bool == "1"){
            curr_block_filename = *it_blocks_file;
            lsplit = split(curr_block_filename, '_');
            hash_string = lsplit.back();
            curr_start_block_idx = *it_start_idx;

            // Check if the blocks from current filemap have already been loaded in the map
            // else (1st load of the filemap blocks) read the blocks and save them in the map
            //          if a block is not already loaded, it is possible that the first genome requiring the block to be loaded
            //          is actually the first genome containing the block ==> ends coordinates directly available while reading the block

            it_loaded_map = loaded_blocks_map.find(hash_string); 

            if(it_loaded_map != loaded_blocks_map.end()){
                current_blocks_list = it_loaded_map->second;

                for(it_block_ends = current_blocks_list.begin(); it_block_ends != current_blocks_list.end(); it_block_ends++){
                    current_assemblyspecific_block.block_idx = it_block_ends->block_idx;

                    pos_pair_1 = kmer_pair_map[it_block_ends->start_kmer_pair_feature_tuple];
                    pos_pair_2 = kmer_pair_map[it_block_ends->end_kmer_pair_feature_tuple];
                    pos_1 = pos_pair_1.first;
                    pos_2 = pos_pair_2.first;

                    if(pos_1 == pos_2){
                        // Block with single kmer pair
                        // By convention, the first occurence of the block is considered to be in the forward strand
                        // Check if the strand (pos_pair_1.second) matches that of the pair that was encountered the first time
                        // Note that although the block is to be considered to be in the forward strand, the orientation designated to the pair
                        // at the time of pair construction may not be the same ==> check if the current pair orientation matches that orientation
                        // If it does, then the block is on the forward strand for the current assembly, else it is in the reverse strand.
                        
                        current_assemblyspecific_block.orientation = (it_block_ends->ref_kmer_pair_strand == pos_pair_1.second);
                    }
                    else
                        current_assemblyspecific_block.orientation = (pos_1<pos_2);

                    if(current_assemblyspecific_block.orientation){
                        current_assemblyspecific_block.start_pos = pos_1;
                        current_assemblyspecific_block.end_pos = pos_2+kmer_len;
                    }
                    else{
                        current_assemblyspecific_block.start_pos = pos_2;
                        current_assemblyspecific_block.end_pos = pos_1+kmer_len;
                    }

                    // if(block_list.size()%2500==0){
                    //     cout<<"\t"<<block_list.size()<<" "<<current_assemblyspecific_block.start_pos<<" ";
                    //     cout<<current_assemblyspecific_block.end_pos<<" "<<current_assemblyspecific_block.orientation<<"\n";
                    // }

                    block_list.push_back(current_assemblyspecific_block);
                }
            }
            else{

                current_blocks_list.clear();
                inFile.open(curr_block_filename,ios::binary);
                inFile.read((char*) (&block_count), sizeof(block_count));
                inFile.read((char*) (&block_assembly_idx), sizeof(block_assembly_idx));

                for(ulli idx=0; idx<block_count; idx++){
                    current_block.block_idx = curr_start_block_idx + idx;
                    current_assemblyspecific_block.block_idx = current_block.block_idx;

                    inFile.read((char*) (&current_block.start_kmer_pair_feature_tuple), sizeof(current_block.start_kmer_pair_feature_tuple));
                    inFile.read((char*) (&current_block.end_kmer_pair_feature_tuple), sizeof(current_block.end_kmer_pair_feature_tuple));

                    inFile.read((char*) (&pos_1), sizeof(pos_1));
                    inFile.read((char*) (&pos_2), sizeof(pos_2));

                    inFile.read((char*) (&ref_kmer_pair_strand), sizeof(ref_kmer_pair_strand));

                    current_block.ref_kmer_pair_strand = ref_kmer_pair_strand;

                    if(block_assembly_idx == assembly_idx){
                        current_assemblyspecific_block.start_pos = pos_1;
                        current_assemblyspecific_block.end_pos = pos_2;
                        // First occurrence of the block is considered to be on the forward strand
                        current_assemblyspecific_block.orientation = true;

                        if( pos_2-pos_1+1 >= min_sv_size){
                            // All blocks above the min_sv_size that occur for first time in the current assembly
                            block_coords = tup_uubu(pos_1, pos_2, true, current_block.block_idx);
                            first_occurrence_blocks_list.push_back(block_coords);
                        }
                    }
                    else{
                        // pos_1 = kmer_pair_map[current_block.start_kmer_pair_feature_tuple];
                        // pos_2 = kmer_pair_map[current_block.end_kmer_pair_feature_tuple];

                        pos_pair_1 = kmer_pair_map[current_block.start_kmer_pair_feature_tuple];
                        pos_pair_2 = kmer_pair_map[current_block.end_kmer_pair_feature_tuple];
                        pos_1 = pos_pair_1.first;
                        pos_2 = pos_pair_2.first;

                        if(pos_1 == pos_2) // Block with single kmer pair
                            current_assemblyspecific_block.orientation = (current_block.ref_kmer_pair_strand == pos_pair_1.second);
                        else
                            current_assemblyspecific_block.orientation = (pos_1<pos_2);

                        if(current_assemblyspecific_block.orientation){
                            current_assemblyspecific_block.start_pos = pos_1;
                            current_assemblyspecific_block.end_pos = pos_2+kmer_len;
                        }
                        else{
                            current_assemblyspecific_block.start_pos = pos_2;
                            current_assemblyspecific_block.end_pos = pos_1+kmer_len;
                        }
                    }

                    current_blocks_list.push_back(current_block);
                    block_list.push_back(current_assemblyspecific_block);
                }

                loaded_blocks_map[hash_string] = current_blocks_list;
                inFile.close();
            }
            
        }
        // else{
        //     // cout<<"\t"<<outstring<<"\n\n";
        //     missing_idx = *it_start_idx;
        //     if(next(it_start_idx) != all_start_block_idx_list.end()){
        //         while(missing_idx < *(next(it_start_idx))){
        //             outstring += "0,0,0,";
        //             missing_idx++;
        //         }
        //     }
        //     else{
        //         while(missing_idx<total_block_count){
        //             outstring += "0,0,0,";
        //             missing_idx++;
        //         }
        //     }
        //     // cout<<"\t\t"<<outstring<<"\n\n\n";
        // }

        it_blocks_file++;
        it_start_idx++;
    }

    cout<<"individual blocks located: "<<assembly_idx<<": "<<block_list.size()<<"\n";

    remove(bool_col_filename);

    ulli prev_block_idx = 0, diff;
    list<Assemblyspecific_Block>::iterator it_blocks = block_list.begin();

    if(block_list.size()>0){
        while(prev_block_idx < it_blocks->block_idx){
            outstring += "0,0,0,";
            prev_block_idx++;
        }

        while(it_blocks != block_list.end()){
            diff = it_blocks->block_idx - prev_block_idx;
            while(diff>1){
                outstring += "0,0,0,";
                diff--;
            }
            outstring += to_string(it_blocks->start_pos) + ",";
            outstring += to_string(it_blocks->end_pos) + ",";
            outstring += (it_blocks->orientation)?"1,":"0,";

            prev_block_idx = it_blocks->block_idx;
            it_blocks++;
        }
    }

    while(prev_block_idx++ < total_block_count)
        outstring += "0,0,0,";

    outstring.pop_back();
    outstring += "\n";
    outFile.open(assembly_blocks_out_dir + to_string(assembly_idx), ios::out);
    outFile << outstring;
    outFile.close();

    ulli block_idx, count=0;
    block_count = first_occurrence_blocks_list.size();
    outFile.open(assembly_blocks_out_dir + "fo_"+ to_string(assembly_idx), ios::binary);
    outFile.write((char*) (&block_count), sizeof(block_count));

    // first_occurrence_blocks_list is sorted on block_idx by design

    for(it_first_occurrence_blocks_list = first_occurrence_blocks_list.begin();
            it_first_occurrence_blocks_list != first_occurrence_blocks_list.end(); it_first_occurrence_blocks_list++){

        coords = tup_uub( get<0>( *it_first_occurrence_blocks_list), get<1>( *it_first_occurrence_blocks_list),
                        get<2>( *it_first_occurrence_blocks_list));
        block_idx = get<3>( *it_first_occurrence_blocks_list);

        outFile.write((char*)(&block_idx), sizeof(block_idx));
        outFile.write((char*) (&coords), sizeof(coords));

        if(count++%100==0){
            cout<< block_idx << " (" << get<0>( *it_first_occurrence_blocks_list) << ",";
            cout<< get<1>( *it_first_occurrence_blocks_list) <<","<< get<2>( *it_first_occurrence_blocks_list) << ") ";
        }
    }
    outFile.close();

    cout<<"locate_and_save_assemblyspecific_collinear_blocks ended: "<<assembly_idx<<": "<<block_list.size()<<"\n";
}


/*vector< ulli > get_contig_ends_vector(string contig_len_filename){
    list<string> lines = get_lines(contig_len_filename.c_str());
    list<string> lsplit;
    vector< ulli > length_ends_vector;
    unsigned int contig_count = 0;
    ulli contig_end;

    for(list<string> :: iterator line_it = lines.begin(); line_it != lines.end(); line_it++){
        lsplit = split(*line_it, '\t');
        contig_end = strtoulli(lsplit.back());
        length_ends_vector.push_back(contig_end);
        contig_count++;
    }
    return length_ends_vector;
}*/


vector< ulli > get_contig_ends_vector(string contig_len_filename, map<string, unsigned int>& contig_name_map){
    list<string> lines = get_lines(contig_len_filename.c_str());
    list<string> lsplit;
    vector< ulli > length_ends_list;
    unsigned int contig_count = 0;
    string contig_name;
    ulli contig_end;

    for(list<string> :: iterator line_it = lines.begin(); line_it != lines.end(); line_it++){
        lsplit = split(*line_it, '\t');
        contig_name = lsplit.front().substr(1, lsplit.front().length());
        contig_end = strtoulli(lsplit.back());
        contig_name_map.insert(psui(contig_name, contig_count));
        length_ends_list.push_back(contig_end);
        // cout<<contig_count<<" "<<contig_name<<" "<<contig_end<<"\n";
        contig_count++;
    }
    // cout<<"get_contig_lengths_list: "<<contig_len_filename<<" "<<contig_name_map.size()<<"\n\n";

    // for(const auto &keyval_pair : contig_name_map){
    //     cout<< keyval_pair.first << " : " << keyval_pair.second << "\n";
    // }

    return length_ends_list;
}


bool block_list_comparator(Assemblyspecific_Block& block_1, Assemblyspecific_Block& block_2){
    return block_1.start_pos < block_2.start_pos;
}


bool block_list_index_comparator(Assemblyspecific_Block& block_1, Assemblyspecific_Block& block_2){
    return block_1.block_idx < block_2.block_idx;
}


bool get_contig_orientation(string orientation){
    return orientation=="B";
}


list<tup_ibibffs> get_metacarvel_oriented_links_tup_list(string oriented_links_filename, map<string, unsigned int>& contig_name_map,
                                                        bool filter_rows=true, float max_overlap=-1000.0){
    cout<<"get_metacarvel_oriented_links_tup_list started: "<<oriented_links_filename<<" "<<contig_name_map.size()<<"\n";
    list<tup_ibibffs> oriented_links_list;
    if(oriented_links_filename.length()<1)
        return oriented_links_list;

    list<string> lines = get_lines(oriented_links_filename.c_str());
    list<string> lsplit;
    list<string>::iterator col_it;
    unsigned int contig_1_idx, contig_2_idx;
    bool contig_1_rel_orientation, contig_2_rel_orientation;
    float separation_mean, separation_std;
    short mapped_links;
    map<string, unsigned int>::iterator contig_map_iter;
    tup_ibibffs current_row;

    for(list<string>::iterator row_it=lines.begin(); row_it!=lines.end(); row_it++){
        lsplit = split(*row_it, '\t');
        col_it = lsplit.begin();
        contig_map_iter = contig_name_map.find(*col_it++);
        if(contig_map_iter == contig_name_map.end()){
            cout<<"COL 1:"<< *(prev(col_it)) <<" NOT FOUND!!!!!!!!!!!!!!!!!\n";
        }
        contig_1_idx = contig_map_iter->second;
        contig_1_rel_orientation = get_contig_orientation(*col_it++);
        contig_map_iter = contig_name_map.find(*col_it++);
        if(contig_map_iter == contig_name_map.end()){
            cout<<"COL 2:"<< *(prev(col_it)) <<" NOT FOUND!!!!!!!!!!!!!!!!!\n";
        }
        contig_2_idx = contig_map_iter->second;
        contig_2_rel_orientation = get_contig_orientation(*col_it++);
        separation_mean = stof(*col_it++);
        separation_std = stof(*col_it++);
        mapped_links = stoi(*col_it++);
        current_row = make_tuple(contig_1_idx, contig_1_rel_orientation, contig_2_idx,
                                contig_2_rel_orientation, separation_mean, separation_std, mapped_links);

        if(filter_rows){
            if(separation_mean>=max_overlap){
                oriented_links_list.push_back(current_row);
            }
            else{
                cout<<"Discarded: ";
                cout<<contig_1_idx<<"  "<<contig_1_rel_orientation<<"  ";
                cout<<contig_2_idx<<"  "<<contig_2_rel_orientation<<"  ";
                cout<<separation_mean<<"  "<<separation_std<<"  "<<mapped_links<<"\n";
            }
        }
        else{
            oriented_links_list.push_back(current_row);
        }
    }
    cout<<"get_metacarvel_oriented_links_tup_list ended: "<<oriented_links_filename<<" "<<contig_name_map.size()<<" "<<oriented_links_list.size()<<"\n";
    return oriented_links_list;
}


void generate_and_save_partitioned_k_nearest_neighbour_pairs(unsigned long assembly_idx, list< Assemblyspecific_Block > sorted_block_list,
                                                        vector<ulli> contig_ends_vector, short k_neighbours, int max_separation, short block_pair_partitions,
                                                        ulli blocks_per_partition, string neighbour_pair_partition_outdir){

    cout<<"generate_and_save_partitioned_k_nearest_neighbour_pairs started: "<<assembly_idx<<"\t"<<sorted_block_list.size()<<"\n";

    ulli block_idx_offset=0, block_rel_idx, block_start, block_end, block_idx;
    list< Assemblyspecific_Block >::iterator it_blocks;
    // vector<ulli> contig_ends_vector;
    unsigned int current_contig_idx;
    list< pulliulli > previous_k_blocks;
    list< pulliulli >::iterator it_prev;
    short pair_partition_idx;
    list< pspuu > neighbour_pairs_list;
    list< pspuu >::iterator it_pairs;
    string outstring;
    ofstream outFile;


    // Compute k-neighbouring pairs preceeding each block -> conversely, k neighbouring pairs succeeding the blocks are also computed
    // contig_ends_vector = get_contig_ends_vector(contig_len_dir + to_string(assembly_idx));
    cout<<contig_ends_vector.size()<<"\t";
    current_contig_idx = 0;

    for(it_blocks = sorted_block_list.begin(); it_blocks != sorted_block_list.end(); it_blocks++){
        // block_start = it_blocks->first;
        block_start = it_blocks->start_pos;
        block_end = it_blocks->end_pos;
        block_idx = it_blocks->block_idx;

        while(block_start > contig_ends_vector[current_contig_idx]){
            previous_k_blocks.clear();
            current_contig_idx++;
        }

        pair_partition_idx = static_cast<short>(block_idx/blocks_per_partition);

        for(it_prev = previous_k_blocks.begin(); it_prev != previous_k_blocks.end(); /*it_prev++*/){
            if( /*(block_start > it_prev->second + 2) &&*/ (block_start > it_prev->second + 1 + max_separation) )
                it_prev = previous_k_blocks.erase(it_prev);
            else{
                neighbour_pairs_list.push_back( make_pair(pair_partition_idx, make_pair( block_idx , it_prev->first ) ) );
                neighbour_pairs_list.push_back( make_pair(static_cast<short>( (it_prev->first)/blocks_per_partition ),
                                                    make_pair( it_prev->first, block_idx ) ) );
                it_prev++;
            }
        }

        if(previous_k_blocks.size() >= k_neighbours)
            previous_k_blocks.pop_front();

        previous_k_blocks.push_back( make_pair(block_idx, block_end) );
    }

    cout<<neighbour_pairs_list.size()<<"\t";

    neighbour_pairs_list.sort();
    it_pairs = neighbour_pairs_list.begin();
    for(short count=0; count<5 && it_pairs!=neighbour_pairs_list.end(); it_pairs++, count++){
        cout<<"("<< it_pairs->first << ",(" << (it_pairs->second).first << "," << (it_pairs->second).second <<"))\t";
    }
    cout<<"\n";
    
    pair_partition_idx = 0;
    outstring = "";
    for(it_pairs = neighbour_pairs_list.begin(); it_pairs != neighbour_pairs_list.end(); it_pairs++){
        while(pair_partition_idx != it_pairs->first){
            if(outstring.length()>0){
                outstring.pop_back();
                outstring += "\n";
            }
            outFile.open(neighbour_pair_partition_outdir + to_string(pair_partition_idx) + "_" + to_string(assembly_idx), ios::out);
            outFile << outstring;
            outFile.close();
            pair_partition_idx++;
            outstring = "";
            // outFile.open(neighbour_pair_partition_outdir + to_string(pair_partition_idx) + "_" + to_string(assembly_idx), ios::out);
        }
        outstring += to_string((it_pairs->second).first) + "," + to_string((it_pairs->second).second) + ",";
    }
    for( ;pair_partition_idx<block_pair_partitions; pair_partition_idx++){
        if(outstring.length()>0){
            outstring.pop_back();
            outstring += "\n";
        }
        outFile.open(neighbour_pair_partition_outdir + to_string(pair_partition_idx) + "_" + to_string(assembly_idx), ios::out);
        outFile << outstring;
        outFile.close();
    }

    cout<<"generate_and_save_partitioned_k_nearest_neighbour_pairs ended: "<<assembly_idx<<"\t"<<sorted_block_list.size()<<"\n";
}


void save_oriented_links_list(list< tup_ibibffs > oriented_links_list, string assembly_idx_str, string oriented_links_outdir){
    cout<<"save_oriented_links_list started: "<<assembly_idx_str<<" "<<oriented_links_list.size()<<"\n";

    string oll_filename = oriented_links_outdir + assembly_idx_str;

    ofstream outFile;
    outFile.open(oll_filename, ios::binary);

    ulli oll_count = oriented_links_list.size();
    outFile.write((char*) (&oll_count), sizeof(oll_count));

    for(list< tup_ibibffs >::iterator list_it = oriented_links_list.begin(); list_it != oriented_links_list.end(); list_it++){
        tup_ibibffs &current_row = (*list_it);
        outFile.write((char*) (&current_row), sizeof(current_row));
    }
    outFile.close();
    cout<<"save_oriented_links_list ended: "<<oll_filename<<" "<<oll_count<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_block_and_neighbours_locator.cpp "<<argc<<"\n";
    string bool_col_file_prefix = argv[1];
    string start_assembly_idx_str = argv[2];
    unsigned long start_assembly_idx = stoul(start_assembly_idx_str);
    string end_assembly_idx_str = argv[3];
    unsigned long end_assembly_idx = stoul(end_assembly_idx_str);
    string all_vector_blocks_filename = argv[4];
    string all_vector_start_block_idx_filename = argv[5];
    short feature_partitions = static_cast<short>(stoi(argv[6]));
    short block_pair_partitions = static_cast<short>(stoi(argv[7]));
    ulli total_block_count = strtoulli(argv[8]);
    string binned_kmer_pairs_dir = argv[9];
    string assembly_blocks_out_dir = argv[10];
    short kmer_len = static_cast<short>(stoi(argv[11]));
    string contig_len_dir = argv[12];
    short k_neighbours = static_cast<short>(stoi(argv[13]));
    short max_separation = static_cast<short>(stoi(argv[14]));
    string neighbour_pair_partition_outdir = argv[15];
    unsigned long min_sv_size = stoul(argv[16]);
    bool use_oriented_links = strncmp(argv[17],"1",1)==0;
    string oriented_links_outdir, contig_len_filename;

    list<string> oriented_links_filename_list;
    if(use_oriented_links){
        oriented_links_filename_list = get_filepaths(argv[18]);
        remove(argv[18]);
        oriented_links_outdir = argv[19];
    }

    // cout<<start_assembly_idx_str<<" "<<end_assembly_idx_str<<"\n";


    list<string> all_filename_list;
    list<ulli> all_start_block_idx_list;
    load_all_vector_files(all_vector_blocks_filename.c_str(), all_vector_start_block_idx_filename.c_str(),
                        all_filename_list, all_start_block_idx_list);


    vector<ulli> bin_arange;
    ulli step = (total_block_count+block_pair_partitions-1)/block_pair_partitions;
    for(ulli i = step; i < total_block_count; i += step)
        bin_arange.push_back(i);
    bin_arange.push_back(total_block_count);

    for(auto check = bin_arange.begin(); check != bin_arange.end(); check++)
        cout<< *check <<" ";
    cout<<"\n";


    // list<string> required_vector_file_list;
    // list<ulli> required_start_block_idx_list;
    // map< tup_psupsubb, ulli > kmer_pair_map;
    map< tup_psupsubb, pullib > kmer_pair_map;
    map< string, list<Collinear_Block_ends> > loaded_blocks_map;

    string bool_col_filename;

    list<string>::iterator it_file_list;
    list<ulli>::iterator it_idx_list;
    list<string>::iterator it_ol_names_list = oriented_links_filename_list.begin();
    // list< Assemblywise_Collinear_Block > assemblywise_collinear_block_list;
    list<Assemblyspecific_Block> assemblyspecific_collinear_block_list;
    vector<ulli> contig_ends_vector;
    map<string, unsigned int> contig_name_map;
    list< tup_ibibffs > oriented_links_list;
    string assembly_idx_str;

    ulli blocks_per_partition = ceil(float(total_block_count)/block_pair_partitions);


    // init_map_group_lists(block_pair_map_list, block_pair_partitions);


    for(unsigned long assembly_idx = start_assembly_idx; assembly_idx <= end_assembly_idx; assembly_idx++){
        // kmer_pair_map.clear();
        contig_name_map.clear();
        assemblyspecific_collinear_block_list.clear();

        assembly_idx_str = to_string(assembly_idx);

        bool_col_filename = bool_col_file_prefix + assembly_idx_str + ".txt";

        load_kmer_pair_map_groups(kmer_pair_map, binned_kmer_pairs_dir, assembly_idx_str, feature_partitions);


        locate_and_save_assemblyspecific_collinear_blocks(assemblyspecific_collinear_block_list, all_filename_list, all_start_block_idx_list,
                                                        assembly_idx, bool_col_filename.c_str(), kmer_pair_map, kmer_len, total_block_count,
                                                        assembly_blocks_out_dir, loaded_blocks_map, min_sv_size);
        kmer_pair_map.clear();

        assemblyspecific_collinear_block_list.sort(block_list_comparator);
        cout<<"Sorted collinear block list "<<assembly_idx<<" "<<assemblyspecific_collinear_block_list.size()<<"\n";

        contig_ends_vector = get_contig_ends_vector(contig_len_dir + to_string(assembly_idx), contig_name_map);

        generate_and_save_partitioned_k_nearest_neighbour_pairs(assembly_idx, assemblyspecific_collinear_block_list, contig_ends_vector, k_neighbours,
                                                                max_separation, block_pair_partitions, blocks_per_partition, neighbour_pair_partition_outdir);

        if(use_oriented_links){
            // contig_len_filename = contig_len_dir + assembly_idx_str;
            // contig_length_ends_list = get_contig_lengths_list(contig_len_filename, contig_name_map);
            // cout<<"Contig ends list: "<< contig_len_filename<<" "<<contig_name_map.size()<<"\n";

            try{
                oriented_links_list = get_metacarvel_oriented_links_tup_list(*it_ol_names_list, contig_name_map);
            }
            catch(int e){
                cerr<<"Exception caught. Error in parsing oriented links for "<<assembly_idx_str<<"\n";
            }
            it_ol_names_list++;

            save_oriented_links_list(oriented_links_list, assembly_idx_str, oriented_links_outdir);
        }
        cout<<assembly_idx_str<<" completed\n";
    }

    return 0;
}