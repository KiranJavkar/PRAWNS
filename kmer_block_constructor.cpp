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
typedef tuple<pullib,tup_psupsubb> tup_pullib_tpsupsubb;
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


vector<string> get_kmer_pair_grouped_filenames(const char* FILENAME){
    ifstream input_file(FILENAME);
    string curr_filename;
    vector<string> filenames;
    while (getline( input_file, curr_filename )){
        filenames.push_back(curr_filename);
    }
    return filenames;
}


/*bool compare_kmer_pair_pos_tuple(const tup_ulli_tpsupsubb& a, const tup_ulli_tpsupsubb& b) 
{ 
    return (get<0>(a) <= get<0>(b)); 
}*/ 


bool compare_kmer_pair_pos_tuple(const tup_pullib_tpsupsubb& a, const tup_pullib_tpsupsubb& b) 
{ 
    return ((get<0>(a)).first <= (get<0>(b)).first); 
} 


/*void load_kmer_pair_tuple_list(string input_filename, list<tup_ulli_tpsupsubb>& kmer_pair_location_tuple_list, bool is_initial_load,
                                unsigned long assembly_count, vector<bool>& presence_vector, unsigned long& block_assembly_idx){
    // cout<<"load_kmer_pair_tuple_list started: "<<input_filename<<" "<<block_assembly_idx<<"\n";

    ifstream inFile;
    inFile.open(input_filename,ios::binary);//|ios::in);
    unsigned long current_block_assembly_idx;
    vector<bool> current_presence_vector;
    bool presence;
    tup_psupsubb kmer_pair_feature_tuple;
    ulli list_size, kmer_pair_pos;
    tup_ulli_tpsupsubb kmer_pair_pos_tuple;

    for(unsigned long idx=0; idx<assembly_count; idx++){
        inFile.read((char*) (&presence), sizeof(presence));
        current_presence_vector.push_back(presence);
    }

    inFile.read((char*) (&current_block_assembly_idx), sizeof(current_block_assembly_idx));
    inFile.read((char*) (&list_size), sizeof(list_size));
    // cout<<input_filename<<"\t"<<current_block_assembly_idx<<"\t"<<list_size<<"\n";

    if(is_initial_load){
        presence_vector = current_presence_vector;
        block_assembly_idx = current_block_assembly_idx;
    }
    else{
        assert(block_assembly_idx==current_block_assembly_idx);
        assert(presence_vector==current_presence_vector);
    }

    for (ulli idx=0; idx<list_size; idx++){
        inFile.read((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
        inFile.read((char*) (&kmer_pair_pos), sizeof(kmer_pair_pos));

        kmer_pair_pos_tuple = make_tuple(kmer_pair_pos, kmer_pair_feature_tuple);
        kmer_pair_location_tuple_list.push_back(kmer_pair_pos_tuple);
    }
    inFile.close();
    cout<<"load_kmer_pair_tuple_list ended: "<<input_filename<<" "<<kmer_pair_location_tuple_list.size()<<"\n";
}*/


void load_kmer_pair_tuple_list(string input_filename, list<tup_pullib_tpsupsubb>& kmer_pair_location_tuple_list, bool is_initial_load,
                                unsigned long assembly_count, vector<bool>& presence_vector, unsigned long& block_assembly_idx){
    // cout<<"load_kmer_pair_tuple_list started: "<<input_filename<<" "<<block_assembly_idx<<"\n";

    ifstream inFile;
    inFile.open(input_filename,ios::binary);//|ios::in);
    unsigned long current_block_assembly_idx;
    vector<bool> current_presence_vector;
    bool presence;
    tup_psupsubb kmer_pair_feature_tuple;
    pullib position_strand_pair;
    ulli list_size;
    tup_pullib_tpsupsubb kmer_pair_pos_tuple;

    for(unsigned long idx=0; idx<assembly_count; idx++){
        inFile.read((char*) (&presence), sizeof(presence));
        current_presence_vector.push_back(presence);
    }

    inFile.read((char*) (&current_block_assembly_idx), sizeof(current_block_assembly_idx));
    inFile.read((char*) (&list_size), sizeof(list_size));
    // cout<<input_filename<<"\t"<<current_block_assembly_idx<<"\t"<<list_size<<"\n";

    if(is_initial_load){
        presence_vector = current_presence_vector;
        block_assembly_idx = current_block_assembly_idx;
    }
    else{
        assert(block_assembly_idx==current_block_assembly_idx);
        assert(presence_vector==current_presence_vector);
    }

    for (ulli idx=0; idx<list_size; idx++){
        inFile.read((char*) (&kmer_pair_feature_tuple), sizeof(kmer_pair_feature_tuple));
        inFile.read((char*) (&position_strand_pair), sizeof(position_strand_pair));

        kmer_pair_pos_tuple = make_tuple(position_strand_pair, kmer_pair_feature_tuple);
        kmer_pair_location_tuple_list.push_back(kmer_pair_pos_tuple);
    }
    inFile.close();
    cout<<"load_kmer_pair_tuple_list ended: "<<input_filename<<" "<<kmer_pair_location_tuple_list.size()<<"\n";
}


struct Collinear_Block{
    vector<bool> presence_vector;
    // psulli start_kmer_1, start_kmer_2, end_kmer_1, end_kmer_2;
    unsigned long block_assembly_idx;
    ulli start_pos, end_pos;
    tup_psupsubb start_kmer_pair_feature_tuple, end_kmer_pair_feature_tuple;
    bool ref_kmer_pair_strand; // Required for block with single kmer pair
};


/*Collinear_Block init_collinear_block(vector<bool>& presence_vector, unsigned long block_assembly_idx, tup_ulli_tpsupsubb kmer_pair_pos_tuple){
    Collinear_Block new_block;
    new_block.presence_vector = presence_vector;
    new_block.block_assembly_idx = block_assembly_idx;
    new_block.start_pos = get<0>(kmer_pair_pos_tuple);
    new_block.start_kmer_pair_feature_tuple = get<1>(kmer_pair_pos_tuple);
    // tup_psupsubb kmer_pair_feature_tuple = get<1>(kmer_pair_pos_tuple);
    // new_block.start_kmer_1 = get<0>(kmer_pair_feature_tuple);
    // new_block.start_kmer_2 = get<1>(kmer_pair_feature_tuple);
    return new_block;
}*/


Collinear_Block init_collinear_block(vector<bool>& presence_vector, unsigned long block_assembly_idx, tup_pullib_tpsupsubb kmer_pair_pos_tuple){
    Collinear_Block new_block;
    new_block.presence_vector = presence_vector;
    new_block.block_assembly_idx = block_assembly_idx;
    new_block.start_pos = (get<0>(kmer_pair_pos_tuple)).first;
    new_block.ref_kmer_pair_strand = (get<0>(kmer_pair_pos_tuple)).second;
    new_block.start_kmer_pair_feature_tuple = get<1>(kmer_pair_pos_tuple);
    // tup_psupsubb kmer_pair_feature_tuple = get<1>(kmer_pair_pos_tuple);
    // new_block.start_kmer_1 = get<0>(kmer_pair_feature_tuple);
    // new_block.start_kmer_2 = get<1>(kmer_pair_feature_tuple);
    return new_block;
}


/*list< Collinear_Block > extract_collinear_blocks(list<tup_ulli_tpsupsubb>& kmer_pair_location_tuple_list, vector<bool>& presence_vector,
                                                unsigned long block_assembly_idx, short kmer_len){
    
    assert(kmer_pair_location_tuple_list.size()>0);

    list< Collinear_Block > block_list;
    list<tup_ulli_tpsupsubb>::iterator it_loc_tup = kmer_pair_location_tuple_list.begin();
    Collinear_Block current_block = init_collinear_block(presence_vector, block_assembly_idx, (*it_loc_tup));

    ulli previous_pos = current_block.start_pos;
    ulli current_pos;

    tup_psupsubb kmer_pair_feature_tuple;

    it_loc_tup++;

    while(it_loc_tup != kmer_pair_location_tuple_list.end()){
        current_pos = get<0>(*it_loc_tup);

        if(previous_pos+1 != current_pos){
            current_block.end_pos = (previous_pos+static_cast<ulli>(kmer_len));

            // kmer_pair_feature_tuple = get<1>(*prev(it_loc_tup));
            // current_block.end_kmer_1 = get<0>(kmer_pair_feature_tuple);
            // current_block.end_kmer_2 = get<1>(kmer_pair_feature_tuple);
            current_block.end_kmer_pair_feature_tuple = get<1>(*prev(it_loc_tup));

            // cout<< block_list.size() << ": " << current_block.start_pos << " " << current_block.end_pos <<"\n";

            block_list.push_back(current_block);

            current_block = init_collinear_block(presence_vector, block_assembly_idx, *it_loc_tup);
        }
        previous_pos = current_pos;
        it_loc_tup++;
    }

    current_block.end_pos = (previous_pos+static_cast<ulli>(kmer_len));
    // kmer_pair_feature_tuple = get<1>(*prev(kmer_pair_location_tuple_list.end())); //get<1>(*prev(it_loc_tup));
    // current_block.end_kmer_1 = get<0>(kmer_pair_feature_tuple);
    // current_block.end_kmer_2 = get<1>(kmer_pair_feature_tuple);
    current_block.end_kmer_pair_feature_tuple = get<1>(*prev(it_loc_tup));
    block_list.push_back(current_block);

    cout<<"extract_collinear_blocks: "<<block_assembly_idx<<" "<<block_list.size()<<"\n";

    return block_list;
}*/


list< Collinear_Block > extract_collinear_blocks(list<tup_pullib_tpsupsubb>& kmer_pair_location_tuple_list, vector<bool>& presence_vector,
                                                unsigned long block_assembly_idx, short kmer_len){
    
    assert(kmer_pair_location_tuple_list.size()>0);

    list< Collinear_Block > block_list;
    list<tup_pullib_tpsupsubb>::iterator it_loc_tup = kmer_pair_location_tuple_list.begin();
    Collinear_Block current_block = init_collinear_block(presence_vector, block_assembly_idx, (*it_loc_tup));

    ulli previous_pos = current_block.start_pos;
    ulli current_pos;

    tup_psupsubb kmer_pair_feature_tuple;

    it_loc_tup++;

    while(it_loc_tup != kmer_pair_location_tuple_list.end()){
        current_pos = (get<0>(*it_loc_tup)).first;

        if(previous_pos+1 != current_pos){
            current_block.end_pos = (previous_pos+static_cast<ulli>(kmer_len));

            // kmer_pair_feature_tuple = get<1>(*prev(it_loc_tup));
            // current_block.end_kmer_1 = get<0>(kmer_pair_feature_tuple);
            // current_block.end_kmer_2 = get<1>(kmer_pair_feature_tuple);
            current_block.end_kmer_pair_feature_tuple = get<1>(*prev(it_loc_tup));

            // cout<< block_list.size() << ": " << current_block.start_pos << " " << current_block.end_pos <<"\n";

            block_list.push_back(current_block);

            current_block = init_collinear_block(presence_vector, block_assembly_idx, *it_loc_tup);
        }
        previous_pos = current_pos;
        it_loc_tup++;
    }

    current_block.end_pos = (previous_pos+static_cast<ulli>(kmer_len));
    // kmer_pair_feature_tuple = get<1>(*prev(kmer_pair_location_tuple_list.end())); //get<1>(*prev(it_loc_tup));
    // current_block.end_kmer_1 = get<0>(kmer_pair_feature_tuple);
    // current_block.end_kmer_2 = get<1>(kmer_pair_feature_tuple);
    current_block.end_kmer_pair_feature_tuple = get<1>(*prev(it_loc_tup));
    block_list.push_back(current_block);

    cout<<"extract_collinear_blocks: "<<block_assembly_idx<<" "<<block_list.size()<<"\n";

    return block_list;
}


/*void save_collinear_blocks( list< Collinear_Block >&block_list, vector<bool>& presence_vector, unsigned long block_assembly_idx,
                            string outdir, string hash_string){
    cout<<"save_collinear_blocks started: "<<outdir<<" "<<hash_string<<"\n";

    string block_filename = outdir + "blocks_" + hash_string;
    string filemap_filename = outdir + "filemap_" + hash_string;

    ofstream outFile;
    outFile.open(block_filename, ios::binary);

    ulli block_count = block_list.size();
    ulli total_blocks_length = 0;

    outFile.write((char*) (&block_count), sizeof(block_count));
    outFile.write((char*) (&block_assembly_idx), sizeof(block_assembly_idx));

    for(list< Collinear_Block >::iterator list_it = block_list.begin(); list_it != block_list.end(); list_it++){
        Collinear_Block &current_block = (*list_it);
        // outFile.write((char*) (&current_block), sizeof(current_block));
        
        // outFile.write((char*) (&current_block.start_kmer_1), sizeof(current_block.start_kmer_1));
        // outFile.write((char*) (&current_block.start_kmer_2), sizeof(current_block.start_kmer_2));
        // outFile.write((char*) (&current_block.end_kmer_1), sizeof(current_block.end_kmer_1));
        // outFile.write((char*) (&current_block.end_kmer_2), sizeof(current_block.end_kmer_2));
        outFile.write((char*) (&current_block.start_kmer_pair_feature_tuple), sizeof(current_block.start_kmer_pair_feature_tuple));
        outFile.write((char*) (&current_block.end_kmer_pair_feature_tuple), sizeof(current_block.end_kmer_pair_feature_tuple));
        outFile.write((char*) (&current_block.start_pos), sizeof(current_block.start_pos));
        outFile.write((char*) (&current_block.end_pos), sizeof(current_block.end_pos));
        total_blocks_length += (current_block.end_pos - current_block.start_pos + 1);
    }

    cout<<"save_collinear_blocks: "<<block_filename<<" "<<block_count<<" "<<total_blocks_length<<"\n";

    ofstream outFile_2;
    outFile_2.open(filemap_filename, ios::out);

    bool presence;

    for(unsigned long idx=0; idx<presence_vector.size(); idx++){
        presence = presence_vector[idx];
        outFile_2 << (presence?"1":"0") << ",";
        outFile.write((char*) (&presence), sizeof(presence));
    }
    outFile_2 << block_filename << "," << block_count << "," << total_blocks_length << "\n";

    cout<<"save_collinear_blocks ended: "<<filemap_filename<<"\n";
}*/


void save_collinear_blocks( list< Collinear_Block >&block_list, vector<bool>& presence_vector, unsigned long block_assembly_idx,
                            string outdir, string hash_string){
    cout<<"save_collinear_blocks started: "<<outdir<<" "<<hash_string<<"\n";

    string block_filename = outdir + "blocks_" + hash_string;
    string filemap_filename = outdir + "filemap_" + hash_string;

    ofstream outFile;
    outFile.open(block_filename, ios::binary);

    ulli block_count = block_list.size();
    ulli total_blocks_length = 0;

    outFile.write((char*) (&block_count), sizeof(block_count));
    outFile.write((char*) (&block_assembly_idx), sizeof(block_assembly_idx));

    for(list< Collinear_Block >::iterator list_it = block_list.begin(); list_it != block_list.end(); list_it++){
        Collinear_Block &current_block = (*list_it);
        outFile.write((char*) (&current_block.start_kmer_pair_feature_tuple), sizeof(current_block.start_kmer_pair_feature_tuple));
        outFile.write((char*) (&current_block.end_kmer_pair_feature_tuple), sizeof(current_block.end_kmer_pair_feature_tuple));
        outFile.write((char*) (&current_block.start_pos), sizeof(current_block.start_pos));
        outFile.write((char*) (&current_block.end_pos), sizeof(current_block.end_pos));
        outFile.write((char*) (&current_block.ref_kmer_pair_strand), sizeof(current_block.ref_kmer_pair_strand));
        total_blocks_length += (current_block.end_pos - current_block.start_pos + 1);
    }

    cout<<"save_collinear_blocks: "<<block_filename<<" "<<block_count<<" "<<total_blocks_length<<"\n";

    ofstream outFile_2;
    outFile_2.open(filemap_filename, ios::out);

    bool presence;

    for(unsigned long idx=0; idx<presence_vector.size(); idx++){
        presence = presence_vector[idx];
        outFile_2 << (presence?"1":"0") << ",";
        outFile.write((char*) (&presence), sizeof(presence));
    }
    outFile_2 << block_filename << "," << block_count << "," << total_blocks_length << "\n";

    cout<<"save_collinear_blocks ended: "<<filemap_filename<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_block_constructor.cpp "<<argc<<"\n";
    string hash_string = argv[1];
    string hash_paths_file = argv[2];
    unsigned long assembly_count = stoul(argv[3]);
    short kmer_len = static_cast<short>(stoi(argv[4]));
    string blocks_output_dir = argv[5];

    // list<tup_ulli_tpsupsubb> kmer_pair_location_tuple_list;
    list<tup_pullib_tpsupsubb> kmer_pair_location_tuple_list;
    vector<bool> presence_vector;
    unsigned long block_assembly_idx;

    vector<string> grouped_kmer_pair_files = get_kmer_pair_grouped_filenames(hash_paths_file.c_str()); //glob(hash_string);
    for(int i=0; i<grouped_kmer_pair_files.size(); i++){
        cout<<grouped_kmer_pair_files[i]<<"\n";
        load_kmer_pair_tuple_list(grouped_kmer_pair_files[i], kmer_pair_location_tuple_list, (i==0), assembly_count, presence_vector,
                                block_assembly_idx);

        // Remove grouped kmer pair file
        remove(grouped_kmer_pair_files[i].c_str());
    }
    remove(hash_paths_file.c_str());

    // cout<<block_assembly_idx<<"\t"<<kmer_pair_location_tuple_list.size()<<"\n";
    // for(unsigned long idx=0; idx<assembly_count; idx++)
    //     cout<<presence_vector[idx]<<" ";
    // cout<<'\n';

    // ulli counter=0;
    // list<tup_ulli_tpsupsubb>::iterator it;
    // for(it = kmer_pair_location_tuple_list.begin(); it != kmer_pair_location_tuple_list.end() && counter<10; it++, counter++){
    //     cout<<"("<<get<0>(*it)<<", ";
    //     tup_psupsubb kmer_pair_feature_tuple = get<1>(*it);
    //     cout<<"( ("<<get<0>(kmer_pair_feature_tuple).first<<","<<get<0>(kmer_pair_feature_tuple).second<<"),";
    //     cout<<" ("<<get<1>(kmer_pair_feature_tuple).first<<","<<get<1>(kmer_pair_feature_tuple).second<<"),";
    //     cout<<" "<<get<2>(kmer_pair_feature_tuple)<<", "<<get<3>(kmer_pair_feature_tuple)<<" )"<<"\n";
    // }

    kmer_pair_location_tuple_list.sort(compare_kmer_pair_pos_tuple);
    // cout<<"\nSORTED\n\n";

    // for(it = kmer_pair_location_tuple_list.begin(); it != kmer_pair_location_tuple_list.end() /*&& counter<10*/; it++/*, counter++*/){
    //     cout<<"("<<get<0>(*it)<<", ";
    //     tup_psupsubb kmer_pair_feature_tuple = get<1>(*it);
    //     cout<<"( ("<<get<0>(kmer_pair_feature_tuple).first<<","<<get<0>(kmer_pair_feature_tuple).second<<"),";
    //     cout<<" ("<<get<1>(kmer_pair_feature_tuple).first<<","<<get<1>(kmer_pair_feature_tuple).second<<"),";
    //     cout<<" "<<get<2>(kmer_pair_feature_tuple)<<", "<<get<3>(kmer_pair_feature_tuple)<<" )"<<"\n";
    // }

    list< Collinear_Block > collinear_block_list = extract_collinear_blocks(kmer_pair_location_tuple_list, presence_vector, block_assembly_idx,
                                                                            kmer_len);
    // list< Collinear_Block >::iterator it_cb;

    // for(it_cb = collinear_block_list.begin(); it_cb != collinear_block_list.end(); it_cb++){
    //     cout<< "Collinear Block: " << it_cb->start_pos << "\t" << it_cb->end_pos << "\t" << (it_cb->end_pos - it_cb->start_pos + 1) << "\n";
    //     cout<<"\tStart: ("<< (it_cb->start_kmer_1).first << "," << (it_cb->start_kmer_1).second << "),";
    //     cout<<" ("<< (it_cb->start_kmer_2).first << "," << (it_cb->start_kmer_2).second << ")\n";
    //     cout<<"\tEnd: ("<< (it_cb->end_kmer_1).first << "," << (it_cb->end_kmer_1).second << "),";
    //     cout<<" ("<< (it_cb->end_kmer_2).first << "," << (it_cb->end_kmer_2).second << ")\n\n";
    // }

    save_collinear_blocks(collinear_block_list, presence_vector, block_assembly_idx, blocks_output_dir, hash_string);

    return 0;
}