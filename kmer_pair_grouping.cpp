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
typedef pair<unsigned long,pullib> pulpullib;
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


/*void set_kmer_pair_presence_vector_multimap(multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map,
                                            map< tup_psupsubb, pululli >& kmer_pair_feature_locator_map,
                                            tup_psupsubb kmer_pair_feature_tuple, ulli pair_start_pos,
                                            unsigned long assembly_idx, unsigned long assembly_count){

    multimap< tup_psupsubb, vector<bool> >::iterator it = kmer_pair_presence_map.find(kmer_pair_feature_tuple);

    if(it!=kmer_pair_presence_map.end()){
        (it->second)[assembly_idx] = true;
    }
    else{
        vector<bool> presence_vector(assembly_count);
        presence_vector[assembly_idx] = true;
        kmer_pair_presence_map.insert(ptuppsuvecb(kmer_pair_feature_tuple, presence_vector));

        // Map the kmer pair to the first assembly containing the pair and the location of occurence in the same
        kmer_pair_feature_locator_map[kmer_pair_feature_tuple] = pululli(assembly_idx, pair_start_pos);
    }
}*/


void set_kmer_pair_presence_vector_multimap(multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map,
                                            map< tup_psupsubb, pulpullib >& kmer_pair_feature_locator_map,
                                            tup_psupsubb kmer_pair_feature_tuple, pullib position_strand_pair,
                                            unsigned long assembly_idx, unsigned long assembly_count){

    multimap< tup_psupsubb, vector<bool> >::iterator it = kmer_pair_presence_map.find(kmer_pair_feature_tuple);

    if(it!=kmer_pair_presence_map.end()){
        (it->second)[assembly_idx] = true;
    }
    else{
        vector<bool> presence_vector(assembly_count);
        presence_vector[assembly_idx] = true;
        kmer_pair_presence_map.insert(ptuppsuvecb(kmer_pair_feature_tuple, presence_vector));

        // Map the kmer pair to the first assembly containing the pair and the location of occurence in the same
        kmer_pair_feature_locator_map[kmer_pair_feature_tuple] = pulpullib(assembly_idx, position_strand_pair);
    }
}


/*void load_binned_kmer_pairs_and_update_presence_vector( multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map,
                                                        map< tup_psupsubb, pululli >& kmer_pair_feature_locator_map,
                                                        string binned_kmer_pairs_dir, unsigned long assembly_count, string map_group_idx_str){
    // cout<<"load_binned_kmer_pairs_and_update_presence_vector started: "<<map_group_idx_str<<"\n";
    tup_psupsubb kmer_pair_feature_tuple;
    ulli pair_start_pos, map_feature_count;
    ifstream inFile;

    for(unsigned long assembly_idx=0; assembly_idx<assembly_count; assembly_idx++){
        inFile.open(binned_kmer_pairs_dir + to_string(assembly_idx) + "_" + map_group_idx_str, ios::binary);
        inFile.read((char*) (&map_feature_count), sizeof(map_feature_count));

        while(map_feature_count--){
            inFile.read((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
            inFile.read((char*) (&pair_start_pos), sizeof(pair_start_pos));

            set_kmer_pair_presence_vector_multimap( kmer_pair_presence_map, kmer_pair_feature_locator_map,
                                                    kmer_pair_feature_tuple, pair_start_pos, assembly_idx, assembly_count);
        }
        inFile.close();
    }
    cout<<"load_binned_kmer_pairs_and_update_presence_vector ended: "<<map_group_idx_str<<":"<<kmer_pair_feature_locator_map.size()<<"\n";
}*/


void load_binned_kmer_pairs_and_update_presence_vector( multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map,
                                                        map< tup_psupsubb, pulpullib >& kmer_pair_feature_locator_map,
                                                        string binned_kmer_pairs_dir, unsigned long assembly_count, string map_group_idx_str){
    // cout<<"load_binned_kmer_pairs_and_update_presence_vector started: "<<map_group_idx_str<<"\n";
    tup_psupsubb kmer_pair_feature_tuple;
    pullib position_strand_pair;
    ulli map_feature_count;
    ifstream inFile;

    for(unsigned long assembly_idx=0; assembly_idx<assembly_count; assembly_idx++){
        inFile.open(binned_kmer_pairs_dir + to_string(assembly_idx) + "_" + map_group_idx_str, ios::binary);
        inFile.read((char*) (&map_feature_count), sizeof(map_feature_count));

        while(map_feature_count--){
            inFile.read((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
            inFile.read((char*) (&position_strand_pair), sizeof(position_strand_pair));

            set_kmer_pair_presence_vector_multimap( kmer_pair_presence_map, kmer_pair_feature_locator_map, kmer_pair_feature_tuple,
                                                    position_strand_pair, assembly_idx, assembly_count);
        }
        inFile.close();
    }
    cout<<"load_binned_kmer_pairs_and_update_presence_vector ended: "<<map_group_idx_str<<":"<<kmer_pair_feature_locator_map.size()<<"\n";
}


/*void filter_min_presence_kmer_pairs(multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map,
                                    map< tup_psupsubb, pululli >& kmer_pair_feature_locator_map,
                                    unsigned long min_presence_count){
    // cout<<"filter_min_presence_kmer_pairs started: "<<min_presence_count<<"\n";

    multimap< tup_psupsubb, vector<bool> >::iterator map_it = kmer_pair_presence_map.begin();
    map< tup_psupsubb, pululli >::iterator feat_loc_map_it = kmer_pair_feature_locator_map.begin();

    tup_psupsubb kmer_pair_feature_tuple;
    vector<bool> current_presence_vector;
    unsigned long presence_count;

    while(map_it != kmer_pair_presence_map.end()){
        kmer_pair_feature_tuple = map_it->first;
        current_presence_vector = map_it->second;
        presence_count = count(current_presence_vector.begin(), current_presence_vector.end(), true);

        if(presence_count<min_presence_count){
            map_it = kmer_pair_presence_map.erase(map_it);
            feat_loc_map_it = kmer_pair_feature_locator_map.erase(feat_loc_map_it);
        }
        else{
            map_it++;
            feat_loc_map_it++;
        }
    }
    cout<<"filter_min_presence_kmer_pairs ended: "<<min_presence_count<<" "<<kmer_pair_feature_locator_map.size()<<"\n";
}*/


void filter_min_presence_kmer_pairs(multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map,
                                    map< tup_psupsubb, pulpullib >& kmer_pair_feature_locator_map,
                                    unsigned long min_presence_count){
    // cout<<"filter_min_presence_kmer_pairs started: "<<min_presence_count<<"\n";

    multimap< tup_psupsubb, vector<bool> >::iterator map_it = kmer_pair_presence_map.begin();
    map< tup_psupsubb, pulpullib >::iterator feat_loc_map_it = kmer_pair_feature_locator_map.begin();

    tup_psupsubb kmer_pair_feature_tuple;
    vector<bool> current_presence_vector;
    unsigned long presence_count;

    while(map_it != kmer_pair_presence_map.end()){
        kmer_pair_feature_tuple = map_it->first;
        current_presence_vector = map_it->second;
        presence_count = count(current_presence_vector.begin(), current_presence_vector.end(), true);

        if(presence_count<min_presence_count){
            map_it = kmer_pair_presence_map.erase(map_it);
            feat_loc_map_it = kmer_pair_feature_locator_map.erase(feat_loc_map_it);
        }
        else{
            map_it++;
            feat_loc_map_it++;
        }
    }
    cout<<"filter_min_presence_kmer_pairs ended: "<<min_presence_count<<" "<<kmer_pair_feature_locator_map.size()<<"\n";
}


void group_kmer_pairs_by_presence_vector(multimap< vector<bool>, list<tup_psupsubb> >& presence_block_kmer_pair_map,
                                        multimap< tup_psupsubb, vector<bool> >& kmer_pair_presence_map){

    // cout<<"group_kmer_pairs_by_presence_vector started:\n";

    multimap< vector<bool>, list<tup_psupsubb> >::iterator block_map_it;
    tup_psupsubb kmer_pair_feature_tuple;
    vector<bool> current_presence_vector;

    for(multimap< tup_psupsubb, vector<bool> >::iterator map_it = kmer_pair_presence_map.begin(); map_it != kmer_pair_presence_map.end(); map_it++){
        kmer_pair_feature_tuple = map_it->first;
        current_presence_vector = map_it->second;
        block_map_it = presence_block_kmer_pair_map.find(current_presence_vector);

        if(block_map_it!=presence_block_kmer_pair_map.end()){
            (block_map_it->second).push_back(kmer_pair_feature_tuple);
        }
        else{
            list<tup_psupsubb> block_kmer_pair_list;
            block_kmer_pair_list.push_back(kmer_pair_feature_tuple);
            presence_block_kmer_pair_map.insert( pair< vector<bool>, list<tup_psupsubb> >(current_presence_vector, block_kmer_pair_list));
        }
    }
    cout<<"group_kmer_pairs_by_presence_vector ended: "<<presence_block_kmer_pair_map.size()<<"\n";
}


/*void save_grouped_kmer_pairs(vector<bool>& current_presence_vector, list<tup_psupsubb>& kmer_pair_tuple_list, unsigned long block_assembly_idx,
                            map< tup_psupsubb, pululli >& kmer_pair_feature_locator_map, string outfilename){

    ofstream outFile;
    // cout<<outfilename<<"\t"<<block_assembly_idx<<"\n";
    outFile.open(outfilename,ios::binary);//|ios::out);

    bool presence;

    for(vector<bool>::iterator it = current_presence_vector.begin(); it != current_presence_vector.end(); it++){
        presence = (*it);
        outFile.write((char*) (&presence), sizeof(presence));
    }

    outFile.write((char*) (&block_assembly_idx), sizeof(block_assembly_idx));

    ulli list_size = kmer_pair_tuple_list.size();
    outFile.write((char*) (&list_size), sizeof(list_size));

    map< tup_psupsubb, pululli >::iterator feat_loc_map_it;
    for(list<tup_psupsubb>::iterator it_kmer_pair = kmer_pair_tuple_list.begin(); it_kmer_pair != kmer_pair_tuple_list.end(); it_kmer_pair++){
        tup_psupsubb &kmer_pair_feature_tuple = (*it_kmer_pair);
        feat_loc_map_it = kmer_pair_feature_locator_map.find(kmer_pair_feature_tuple);
        ulli &kmer_pair_pos = (feat_loc_map_it->second).second;
        outFile.write((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
        outFile.write((char*) (&kmer_pair_pos), sizeof(kmer_pair_pos));
    }
}*/


void save_grouped_kmer_pairs(vector<bool>& current_presence_vector, list<tup_psupsubb>& kmer_pair_tuple_list, unsigned long block_assembly_idx,
                            map< tup_psupsubb, pulpullib >& kmer_pair_feature_locator_map, string outfilename){

    ofstream outFile;
    // cout<<outfilename<<"\t"<<block_assembly_idx<<"\n";
    outFile.open(outfilename,ios::binary);//|ios::out);

    bool presence;

    for(vector<bool>::iterator it = current_presence_vector.begin(); it != current_presence_vector.end(); it++){
        presence = (*it);
        outFile.write((char*) (&presence), sizeof(presence));
    }

    outFile.write((char*) (&block_assembly_idx), sizeof(block_assembly_idx));

    ulli list_size = kmer_pair_tuple_list.size();
    outFile.write((char*) (&list_size), sizeof(list_size));

    map< tup_psupsubb, pulpullib >::iterator feat_loc_map_it;
    for(list<tup_psupsubb>::iterator it_kmer_pair = kmer_pair_tuple_list.begin(); it_kmer_pair != kmer_pair_tuple_list.end(); it_kmer_pair++){
        tup_psupsubb &kmer_pair_feature_tuple = (*it_kmer_pair);
        feat_loc_map_it = kmer_pair_feature_locator_map.find(kmer_pair_feature_tuple);
        pullib &position_strand_pair = (feat_loc_map_it->second).second;
        outFile.write((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
        outFile.write((char*) (&position_strand_pair), sizeof(position_strand_pair));
    }
}


/*void save_grouped_kmer_pairs_wrapper(multimap< vector<bool>, list<tup_psupsubb> >& presence_block_kmer_pair_map,
                                    map< tup_psupsubb, pululli >& kmer_pair_feature_locator_map,
                                    string grouped_pairs_out_dir, string map_group_idx_str){
    // cout<<"save_grouped_kmer_pairs_wrapper started: "<<map_group_idx_str<<"\n";

    map< tup_psupsubb, pululli >::iterator feat_loc_map_it;
    unsigned long block_assembly_idx; // The smallest assembly index that has all the features constituting this block
    hash< vector<bool> > hasher;
    vector<bool> current_presence_vector;
    list<tup_psupsubb> kmer_pair_tuple_list;
    string hash_string, outfilename, group_files_outfilename;
    ofstream outFile;
    multimap< vector<bool>, list<tup_psupsubb> >::iterator block_map_it;

    for(block_map_it = presence_block_kmer_pair_map.begin(); block_map_it != presence_block_kmer_pair_map.end(); block_map_it++){
        current_presence_vector = block_map_it->first;
        kmer_pair_tuple_list = block_map_it->second;

        feat_loc_map_it = kmer_pair_feature_locator_map.find(kmer_pair_tuple_list.front());
        block_assembly_idx = (feat_loc_map_it->second).first;

        hash_string = to_string(hasher(current_presence_vector));
        outfilename = grouped_pairs_out_dir + hash_string + "_" + map_group_idx_str;

        save_grouped_kmer_pairs(current_presence_vector, kmer_pair_tuple_list, block_assembly_idx, kmer_pair_feature_locator_map, outfilename);
        
        group_files_outfilename = grouped_pairs_out_dir + "group_files_" + hash_string;
        outFile.open(group_files_outfilename, ios_base::app);
        outFile<<outfilename<<"\n";
        outFile.close();
    }
    cout<<"save_grouped_kmer_pairs_wrapper ended: "<<map_group_idx_str<<"\n";
}*/


void save_grouped_kmer_pairs_wrapper(multimap< vector<bool>, list<tup_psupsubb> >& presence_block_kmer_pair_map,
                                    map< tup_psupsubb, pulpullib >& kmer_pair_feature_locator_map,
                                    string grouped_pairs_out_dir, string map_group_idx_str){
    // cout<<"save_grouped_kmer_pairs_wrapper started: "<<map_group_idx_str<<"\n";

    map< tup_psupsubb, pulpullib >::iterator feat_loc_map_it;
    unsigned long block_assembly_idx; // The smallest assembly index that has all the features constituting this block
    hash< vector<bool> > hasher;
    vector<bool> current_presence_vector;
    list<tup_psupsubb> kmer_pair_tuple_list;
    string hash_string, outfilename, group_files_outfilename;
    ofstream outFile;
    multimap< vector<bool>, list<tup_psupsubb> >::iterator block_map_it;

    for(block_map_it = presence_block_kmer_pair_map.begin(); block_map_it != presence_block_kmer_pair_map.end(); block_map_it++){
        current_presence_vector = block_map_it->first;
        kmer_pair_tuple_list = block_map_it->second;

        feat_loc_map_it = kmer_pair_feature_locator_map.find(kmer_pair_tuple_list.front());
        block_assembly_idx = (feat_loc_map_it->second).first;

        hash_string = to_string(hasher(current_presence_vector));
        outfilename = grouped_pairs_out_dir + hash_string + "_" + map_group_idx_str;

        save_grouped_kmer_pairs(current_presence_vector, kmer_pair_tuple_list, block_assembly_idx, kmer_pair_feature_locator_map, outfilename);
        
        group_files_outfilename = grouped_pairs_out_dir + "group_files_" + hash_string;
        outFile.open(group_files_outfilename, ios_base::app);
        outFile<<outfilename<<"\n";
        outFile.close();
    }
    cout<<"save_grouped_kmer_pairs_wrapper ended: "<<map_group_idx_str<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_pair_grouping.cpp "<<argc<<"\n";
    string map_group_idx_str = argv[1];
    string binned_kmer_pairs_dir = argv[2];
    unsigned long assembly_count = stoul(argv[3]);
    string grouped_pairs_out_dir = argv[4];
    // Only filter out kmer pairs below a certain value
    // kmer pairs present in (nearly) all genomes can still have different phenotypic associations based on blocks collocation
    unsigned long min_presence_count = stoul(argv[5]);
    short kmer_len = static_cast<short>(stoi(argv[6]));

    multimap< tup_psupsubb, vector<bool> > kmer_pair_presence_map;
    // map< tup_psupsubb, pululli > kmer_pair_feature_locator_map;
    map< tup_psupsubb, pulpullib > kmer_pair_feature_locator_map;
    multimap< vector<bool>, list<tup_psupsubb> > presence_block_kmer_pair_map;


    load_binned_kmer_pairs_and_update_presence_vector(  kmer_pair_presence_map, kmer_pair_feature_locator_map, binned_kmer_pairs_dir,
                                                        assembly_count, map_group_idx_str);

    filter_min_presence_kmer_pairs(kmer_pair_presence_map, kmer_pair_feature_locator_map, min_presence_count);

    group_kmer_pairs_by_presence_vector(presence_block_kmer_pair_map, kmer_pair_presence_map);

    save_grouped_kmer_pairs_wrapper(presence_block_kmer_pair_map, kmer_pair_feature_locator_map, grouped_pairs_out_dir, map_group_idx_str);

    return 0;
}