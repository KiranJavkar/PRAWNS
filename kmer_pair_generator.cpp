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
typedef tuple<ulli,ulli,int,bool> tup_uuib;
typedef tuple<ulli,bool,int> tup_ullibi;
typedef tuple<tup_ullibi,tup_ullibi> tup_tubitubi;
typedef tuple<unsigned long,bool,int> tup_ulbi;
typedef tuple<tup_ulbi,tup_ulbi> tup_tulbitulbi;
typedef pair<tup_ululbb,vector<bool> > ptupulvecb;
typedef pair<short,pullib> ps_pullib;
typedef pair<ulli,pullib> pullipullib;


// list< ulli > get_suffix_list(string filename){
//     ifstream inFile;
//     inFile.open(filename,ios::binary);//|ios::in);
//     ulli retained_count;
//     inFile.read((char*) (&retained_count), sizeof(retained_count));
//     list< ulli > suffix_list;
//     for (ulli idx=0; idx<retained_count; idx++){
//         ulli suffix_idx;
//         inFile.read((char*) (&suffix_idx), sizeof(suffix_idx));
//         suffix_list.push_back(suffix_idx);
//     }
//     inFile.close();
//     return suffix_list;
// }


string get_nucleotide_string(ulli value, short length){
    string result = string(length, 'A');
    int pos = length-1;
    while(value>0){
        char nuc;
        switch(value%4){
            case 0: nuc='A'; break;
            case 1: nuc='C'; break;
            case 2: nuc='G'; break;
            case 3: nuc='T'; break;
        }
        result[pos--] = nuc;
        value = value>>2;
    }
    return result;
}


string get_kmer_string(int prefix_idx, ulli suffix_idx, short kmer_len=25, short prefix_len=5){
    return get_nucleotide_string(prefix_idx, prefix_len) + get_nucleotide_string(suffix_idx, kmer_len-prefix_len);
}


void get_kmer_position_tuple_list(list< tup_upsub >& kmer_position_tuple_list, string binned_kmer_dir, string filtered_kmer_dir,
                                    string assembly_idx_str, short prefix_count){
    cout<<"get_kmer_position_tuple_list started: "<<assembly_idx_str<<"\n";
    // list< tup_upsub > kmer_position_tuple_list;
    ulli filtered_suffix_idx, assembly_suffix_idx, filtered_count, binned_count, check, current_filtered_loc, current_binned_loc;
    pullib pos_strand_pair;
    ifstream inFile_filtered, inFile_binned_mer, inFile_binned_pos;
    string binned_kmer_file_prefix = binned_kmer_dir + assembly_idx_str + "_";

    for(short prefix_idx=0; prefix_idx<prefix_count; prefix_idx++){
        // cout<<binned_kmer_file_prefix<<prefix_idx<<"\n";
        inFile_filtered.open(filtered_kmer_dir + to_string(prefix_idx), ios::binary);
        inFile_binned_mer.open(binned_kmer_file_prefix + to_string(prefix_idx) + "_mer", ios::binary);
        inFile_binned_pos.open(binned_kmer_file_prefix + to_string(prefix_idx) + "_pos", ios::binary);

        inFile_filtered.read((char*) (&filtered_count), sizeof(filtered_count));
        inFile_binned_mer.read((char*) (&binned_count), sizeof(binned_count));
        inFile_binned_pos.read((char*) (&check), sizeof(check));

        assert(binned_count==check);

        current_filtered_loc = 0;
        current_binned_loc = 0;
        inFile_filtered.read((char*) (&filtered_suffix_idx), sizeof(filtered_suffix_idx));
        inFile_binned_mer.read((char*) (&assembly_suffix_idx), sizeof(assembly_suffix_idx));
        inFile_binned_pos.read((char*) (&pos_strand_pair), sizeof(pos_strand_pair));

        while(current_filtered_loc<filtered_count && current_binned_loc<binned_count){
            if(filtered_suffix_idx==assembly_suffix_idx){
                kmer_position_tuple_list.push_back( make_tuple( pos_strand_pair.first, // kmer position
                                                                psulli(prefix_idx, assembly_suffix_idx), // kmer prefix-suffix pair
                                                                pos_strand_pair.second)); // kmer orientation

                inFile_filtered.read((char*) (&filtered_suffix_idx), sizeof(filtered_suffix_idx));
                inFile_binned_mer.read((char*) (&assembly_suffix_idx), sizeof(assembly_suffix_idx));
                inFile_binned_pos.read((char*) (&pos_strand_pair), sizeof(pos_strand_pair));
                current_filtered_loc++;
                current_binned_loc++;
            }
            else if(filtered_suffix_idx<assembly_suffix_idx){
                inFile_filtered.read((char*) (&filtered_suffix_idx), sizeof(filtered_suffix_idx));
                current_filtered_loc++;
            }
            else{
                inFile_binned_mer.read((char*) (&assembly_suffix_idx), sizeof(assembly_suffix_idx));
                inFile_binned_pos.read((char*) (&pos_strand_pair), sizeof(pos_strand_pair));
                current_binned_loc++;
            }
        }

        inFile_filtered.close();
        inFile_binned_mer.close();
        inFile_binned_pos.close();

        // // Remove binned kmer files
        // string remove_file_cmd = "rm -rf " + binned_kmer_file_prefix + "*";
        // system(remove_file_cmd.c_str());
    }

    kmer_position_tuple_list.sort();

    cout<<"get_kmer_position_tuple_list ended: "<<assembly_idx_str<<": "<<kmer_position_tuple_list.size()<<"\n";
    // return kmer_position_tuple_list;
}


// Only overlapping kmers are used to form kmer pairs
// kmer pairs can be uniquely located by start position of the pair ==> smaller start position between the constituent kmers
// By using mathematical induction, the kmer pair with consecutive start locations can be merged to form contiguous blocks
/*void init_map_group_lists(list< map< tup_psupsubb, ulli > >& kmer_pair_map_list, short partitions){
    cout<<"init_map_group_lists "<<partitions<<"\n";
    for(short idx=0; idx<partitions; idx++){
        map< tup_psupsubb, ulli > new_kmer_pair_map;
        kmer_pair_map_list.push_back(new_kmer_pair_map);
    }
    // cout<<"init_map_group_lists completed\n";
}*/
void init_map_group_lists(list< map< tup_psupsubb, pullib > >& kmer_pair_map_list, short partitions){
    cout<<"init_map_group_lists "<<partitions<<"\n";
    for(short idx=0; idx<partitions; idx++){
        map< tup_psupsubb, pullib > new_kmer_pair_map;
        kmer_pair_map_list.push_back(new_kmer_pair_map);
    }
    // cout<<"init_map_group_lists completed\n";
}


/*void add_pair_to_map_group(psulli kmer_1_idx, psulli kmer_2_idx, bool kmer_1_strand, bool kmer_2_strand, ulli kmer_pair_start_pos,
                            list< map< tup_psupsubb, ulli > >& kmer_pair_map_list, double prefix_group_size){
    psulli new_kmer_1_idx, new_kmer_2_idx;
    bool kmer_1_relative_orientation_bool, kmer_2_relative_orientation_bool;

    if(kmer_1_idx.first<=kmer_2_idx.first){
        new_kmer_1_idx = kmer_1_idx;
        new_kmer_2_idx = kmer_2_idx;
        kmer_1_relative_orientation_bool = kmer_1_strand;
        kmer_2_relative_orientation_bool = kmer_2_strand;
    }
    else{
        new_kmer_1_idx = kmer_2_idx;
        new_kmer_2_idx = kmer_1_idx;
        kmer_1_relative_orientation_bool = !kmer_2_strand;
        kmer_2_relative_orientation_bool = !kmer_1_strand;
    }

    short map_group_idx = static_cast<short> (log2(new_kmer_1_idx.first)/prefix_group_size);
    // cout<<new_kmer_1_idx.first<<" "<<new_kmer_2_idx.first<<" "<<map_group_idx<<"\n";
    auto map_iter = kmer_pair_map_list.begin();
    advance(map_iter, map_group_idx);
    // Add to map
    (*map_iter)[make_tuple(new_kmer_1_idx, new_kmer_2_idx, kmer_1_relative_orientation_bool, kmer_2_relative_orientation_bool)] = kmer_pair_start_pos;
}*/


void add_pair_to_map_group(psulli kmer_1_idx, psulli kmer_2_idx, bool kmer_1_strand, bool kmer_2_strand, ulli kmer_pair_start_pos,
                            list< map< tup_psupsubb, pullib > >& kmer_pair_map_list, double prefix_group_size){
    psulli new_kmer_1_idx, new_kmer_2_idx;
    bool kmer_1_relative_orientation_bool, kmer_2_relative_orientation_bool, pair_strand;

    if(kmer_1_idx.first < kmer_2_idx.first || (kmer_1_idx.first == kmer_2_idx.first && kmer_1_idx.second <= kmer_2_idx.second)){
        new_kmer_1_idx = kmer_1_idx;
        new_kmer_2_idx = kmer_2_idx;
        kmer_1_relative_orientation_bool = kmer_1_strand;
        kmer_2_relative_orientation_bool = kmer_2_strand;
        pair_strand = true;
    }
    else{
        new_kmer_1_idx = kmer_2_idx;
        new_kmer_2_idx = kmer_1_idx;
        kmer_1_relative_orientation_bool = !kmer_2_strand;
        kmer_2_relative_orientation_bool = !kmer_1_strand;
        pair_strand = false;
    }

    short map_group_idx = static_cast<short> (log2(new_kmer_1_idx.first)/prefix_group_size);
    // cout<<new_kmer_1_idx.first<<" "<<new_kmer_2_idx.first<<" "<<map_group_idx<<"\n";
    auto map_iter = kmer_pair_map_list.begin();
    advance(map_iter, map_group_idx);
    // Add to map
    (*map_iter)[make_tuple(new_kmer_1_idx, new_kmer_2_idx, kmer_1_relative_orientation_bool, kmer_2_relative_orientation_bool)] = pullib(kmer_pair_start_pos, pair_strand);
}


/*void save_retained_kmers_and_extracted_pairs(list< tup_upsub >& kmer_position_tuple_list, short prefix_count, short partitions,
                                            string assembly_idx_str, string retained_kmer_dir, string binned_kmer_pairs_dir){
    list< map< tup_psupsubb, ulli > > kmer_pair_map_list;
    init_map_group_lists(kmer_pair_map_list, partitions);

    psulli kmer_idx, previous_kmer_idx;
    bool strand, previous_strand;
    double prefix_group_size = log2 (prefix_count)/partitions;

    string retained_kmer_filename = retained_kmer_dir + assembly_idx_str, binned_pair_file_prefix = binned_kmer_pairs_dir + assembly_idx_str + "_";

    ofstream outFile, outFile_pair;
    ulli retained_count = kmer_position_tuple_list.size(), position, previous_position=-1;//by design, the smallest possible start position is 1
    cout<<"save_retained_kmers_and_extracted_pairs: "<<retained_kmer_filename<<" "<<retained_count<<"\n";
    outFile.open(retained_kmer_filename,ios::binary);//|ios::out);
    outFile.write((char*) (&retained_count), sizeof(retained_count));

    // short counter = 25;

    for (list<tup_upsub>::iterator it = kmer_position_tuple_list.begin(); it != kmer_position_tuple_list.end(); it++){
        position = get<0>(*it);
        kmer_idx = get<1>(*it);
        strand = get<2>(*it);
        // Save retained kmers
        outFile.write((char*) (&position), sizeof(position));
        outFile.write((char*) (&kmer_idx), sizeof(kmer_idx));

        // Only overlapping kmers are used to form kmer pairs
        // kmer pairs can be uniquely located by start position of the pair ==> smaller start position between the constituent kmers
        // By using mathematical induction, the kmer pair with consecutive start locations can be merged to form contiguous blocks
        if(previous_position+1==position){
            add_pair_to_map_group(previous_kmer_idx, kmer_idx, previous_strand, strand, previous_position, kmer_pair_map_list, prefix_group_size);
        }

        // if(counter>0){
        //     cout<<position<<" ("<<kmer_idx.first<<","<<kmer_idx.second<<") "<<get_kmer_string(kmer_idx.first, kmer_idx.second)<<" "<<strand<<"\n";
        //     counter--;
        // }

        previous_kmer_idx = kmer_idx;
        previous_position = position;
        previous_strand = strand;
    }
    outFile.close();


    // Save binned pairs
    short map_group_idx = 0;
    string outfilename;
    ulli map_feature_count;
    tup_psupsubb kmer_pair_feature_tuple;

    for(auto map_iter = kmer_pair_map_list.begin(); map_iter != kmer_pair_map_list.end(); map_iter++){
        // Do this for each map in the list of maps
        outfilename = binned_pair_file_prefix + to_string(map_group_idx);
        map_feature_count = (*map_iter).size();

        cout<<outfilename<<"\t"<<map_feature_count<<"\n";
        outFile_pair.open(outfilename,ios::binary);//|ios::out);
        outFile_pair.write((char*) (&map_feature_count), sizeof(map_feature_count));

        for(const auto &feature_keyval_pair : (*map_iter)){
            kmer_pair_feature_tuple = feature_keyval_pair.first;
            position = feature_keyval_pair.second;
            outFile_pair.write((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
            outFile_pair.write((char*) (&position), sizeof(position));
        }
        outFile_pair.close();
        map_group_idx++;
    }
}*/


void save_retained_kmers_and_extracted_pairs(list< tup_upsub >& kmer_position_tuple_list, short prefix_count, short partitions,
                                            string assembly_idx_str, string retained_kmer_dir, string binned_kmer_pairs_dir){
    list< map< tup_psupsubb, pullib > > kmer_pair_map_list;
    init_map_group_lists(kmer_pair_map_list, partitions);

    psulli kmer_idx, previous_kmer_idx;
    bool strand, previous_strand;
    double prefix_group_size = log2 (prefix_count)/partitions;

    string retained_kmer_filename = retained_kmer_dir + assembly_idx_str, binned_pair_file_prefix = binned_kmer_pairs_dir + assembly_idx_str + "_";

    ofstream outFile, outFile_pair;
    ulli retained_count = kmer_position_tuple_list.size(), position, previous_position=-1;//by design, the smallest possible start position is 1
    cout<<"save_retained_kmers_and_extracted_pairs: "<<retained_kmer_filename<<" "<<retained_count<<"\n";
    outFile.open(retained_kmer_filename,ios::binary);//|ios::out);
    outFile.write((char*) (&retained_count), sizeof(retained_count));

    // short counter = 25;

    for (list<tup_upsub>::iterator it = kmer_position_tuple_list.begin(); it != kmer_position_tuple_list.end(); it++){
        position = get<0>(*it);
        kmer_idx = get<1>(*it);
        strand = get<2>(*it);
        // Save retained kmers
        outFile.write((char*) (&position), sizeof(position));
        outFile.write((char*) (&kmer_idx), sizeof(kmer_idx));

        // Only overlapping kmers are used to form kmer pairs
        // kmer pairs can be uniquely located by start position of the pair ==> smaller start position between the constituent kmers
        // By using mathematical induction, the kmer pair with consecutive start locations can be merged to form contiguous blocks
        if(previous_position+1==position){
            add_pair_to_map_group(previous_kmer_idx, kmer_idx, previous_strand, strand, previous_position, kmer_pair_map_list, prefix_group_size);
        }

        // if(counter>0){
        //     cout<<position<<" ("<<kmer_idx.first<<","<<kmer_idx.second<<") "<<get_kmer_string(kmer_idx.first, kmer_idx.second)<<" "<<strand<<"\n";
        //     counter--;
        // }

        previous_kmer_idx = kmer_idx;
        previous_position = position;
        previous_strand = strand;
    }
    outFile.close();


    // Save binned pairs
    short map_group_idx = 0;
    string outfilename;
    ulli map_feature_count;
    tup_psupsubb kmer_pair_feature_tuple;
    pullib position_strand_pair;

    for(auto map_iter = kmer_pair_map_list.begin(); map_iter != kmer_pair_map_list.end(); map_iter++){
        // Do this for each map in the list of maps
        outfilename = binned_pair_file_prefix + to_string(map_group_idx);
        map_feature_count = (*map_iter).size();

        cout<<outfilename<<"\t"<<map_feature_count<<"\n";
        outFile_pair.open(outfilename,ios::binary);//|ios::out);
        outFile_pair.write((char*) (&map_feature_count), sizeof(map_feature_count));

        for(const auto &feature_keyval_pair : (*map_iter)){
            kmer_pair_feature_tuple = feature_keyval_pair.first;
            position_strand_pair = feature_keyval_pair.second;
            outFile_pair.write((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
            outFile_pair.write((char*) (&position_strand_pair), sizeof(position_strand_pair));
        }
        outFile_pair.close();
        map_group_idx++;
    }
}


void remove_file(string filename){
    string remove_file_cmd = "rm -rf " + filename;
    system(remove_file_cmd.c_str());
}


void remove_file(string filename_prefix, short prefix_count){
    string prefix_idx_str;
    for (short prefix_idx = 0; prefix_idx<prefix_count; prefix_idx++){
        prefix_idx_str = to_string(prefix_idx);
        remove((filename_prefix+prefix_idx_str+"_mer").c_str());
        remove((filename_prefix+prefix_idx_str+"_pos").c_str());
    }
}


// void save_retained_kmers(list< tup_upsub >& kmer_position_tuple_list, string retained_kmer_filename){
//     ofstream outFile;
//     ulli retained_count = kmer_position_tuple_list.size(), position;
//     psulli kmer_idx;
//     cout<<"save_retained_kmers: "<<retained_kmer_filename<<" "<<retained_count<<"\n";
//     outFile.open(retained_kmer_filename,ios::binary);//|ios::out);
//     outFile.write((char*) (&retained_count), sizeof(retained_count));

//     short counter = 25;

//     for (list<tup_upsub>::iterator it = kmer_position_tuple_list.begin(); it != kmer_position_tuple_list.end(); it++){
//         position = get<0>(*it);
//         kmer_idx = get<1>(*it);
//         outFile.write((char*) (&position), sizeof(position));
//         outFile.write((char*) (&kmer_idx), sizeof(kmer_idx));

//         if(counter>0){
//             cout<<position<<" ("<<kmer_idx.first<<","<<kmer_idx.second<<") "<<get_kmer_string(kmer_idx.first, kmer_idx.second)<<" "<<get<2>(*it)<<"\n";
//             counter--;
//         }
//     }
//     outFile.close();
// }


void kmer_pair_generator_wrapper(string binned_kmer_dir, string filtered_kmer_dir, string assembly_idx_str,
                                short prefix_count, string retained_kmer_dir, string binned_kmer_pairs_dir, short partitions){
    // cout<<"kmer_pair_generator_wrapper started: "<<assembly_idx_str<<"\n";
    list< tup_upsub > kmer_position_tuple_list;
    get_kmer_position_tuple_list(kmer_position_tuple_list, binned_kmer_dir, filtered_kmer_dir, assembly_idx_str, prefix_count);
    // save_retained_kmers(kmer_position_tuple_list, retained_kmer_dir+assembly_idx_str);
    save_retained_kmers_and_extracted_pairs(kmer_position_tuple_list, prefix_count, partitions, assembly_idx_str,
                                            retained_kmer_dir, binned_kmer_pairs_dir);
    // remove_file(binned_kmer_dir + assembly_idx_str + "_*");
    remove_file(binned_kmer_dir + assembly_idx_str + "_", prefix_count);
    // cout<<"kmer_pair_generator_wrapper ended\n";
}


int main(int argc, char** argv){
    cout<<"kmer_pair_generator.cpp "<<argv<<"\n";
    string assembly_idx_str = argv[1];
    short prefix_count = static_cast<short>(stoi(argv[2]));
    string binned_kmer_dir = argv[3];
    string filtered_kmer_dir = argv[4];
    string retained_kmer_dir = argv[5];
    string binned_kmer_pairs_dir = argv[6];
    short partitions = static_cast<short>(stoi(argv[7]));

    // kmer_pair_generator_wrapper(binned_kmer_dir, filtered_kmer_dir, assembly_idx_str, prefix_count, retained_kmer_dir);
    kmer_pair_generator_wrapper(binned_kmer_dir, filtered_kmer_dir, assembly_idx_str, prefix_count, retained_kmer_dir,
                                binned_kmer_pairs_dir, partitions);

    return 0;
}