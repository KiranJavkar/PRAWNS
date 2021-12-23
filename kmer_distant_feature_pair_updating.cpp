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
typedef pair<bool,bool> pbb;
typedef pair<unsigned int,bool> puib;
typedef pair<int,bool> pib;
typedef tuple<ulli,int,bool> tup_uib;
typedef tuple<ulli,ulli,bool> tup_uub;
typedef tuple<ulli,ulli,ulli> tup_uuu;
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
typedef tuple<ulli,ulli,bool,ulli,bool> tup_uubub;
typedef pair<tup_uubb,tup_uubb> ptupuubbtupuubb;
typedef tuple<ulli,bool,ulli,bool,bool,bool,long> tup_ububbbl;
typedef tuple<bool,bool,bool,bool,long> tup_bbbbl;
typedef tuple<string,bool,string,bool,bool,bool> tup_sbsbbb;


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


vector< ulli > get_contig_ends_vector(string contig_len_filename){
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
}


map< puiui, tup_bbffs > load_oriented_links_map(string oll_filename){
    // cout<<"load_oriented_links_map started: "<<oll_filename<<"\n";

    // list< tup_ibibffs > oriented_links_list;
    map< puiui, tup_bbffs> oriented_links_map;
    ulli oll_count;
    tup_ibibffs current_row;

    ifstream inFile;
    inFile.open(oll_filename, ios::binary);

    inFile.read((char*) (&oll_count), sizeof(oll_count));

    while(oll_count--){
        inFile.read((char*) (&current_row), sizeof(current_row));
        // oriented_links_list.push_back(current_row);
        oriented_links_map[ puiui( get<0>(current_row), get<2>(current_row) ) ] = make_tuple(get<1>(current_row), get<3>(current_row),
                                                                                        get<4>(current_row), get<5>(current_row), get<6>(current_row));
    }
    inFile.close();
    // cout<<"load_oriented_links_map ended: "<<oriented_links_map.size()<<"\n";
    return oriented_links_map;
}


unsigned int get_contig_index_by_coordinate(vector<ulli> contig_ends_vector, ulli start_coordinate, unsigned int start, unsigned int end){
    if(start==end)
        return start;
    unsigned int mid_contig = (start+end)/2;
    ulli mid_value = contig_ends_vector[mid_contig];

    if(start_coordinate == mid_value)
        return mid_contig;

    if(start_coordinate > mid_value)
        return get_contig_index_by_coordinate(contig_ends_vector, start_coordinate, mid_contig+1, end);

    return get_contig_index_by_coordinate(contig_ends_vector, start_coordinate, start, mid_contig);
}


// \/ load the feature coords (separated by whether it is a metablock or block) for each assembly chosen within the current partition (delete loaded files)
// \/ load the contig ends vectors and oriented links (if provided) for the corresponding assemblies
// XXXXX load the maps that denote the pairs to check from the respective assemblies--#partitions files
// \/ load all files containing the pair tuples detected
// \/ for the assemblies to be checked, load the pair presence files -- the pairs not marked as present are to be checked
// \/ check if a valid pair can be formed: features should exists and the subsequent pair to be formed should match the required relative orientations
// \/      if pair exists, set the observed separation in the corresponding update map
// \/ once update maps are generated, load the existing pair separation files
// \/ for each entry where the presence is 0 (false), check if an entry exists in the corresponding update map
// \/      if it does, set the presence string to 2 and separation string to updated separation distance
// \/      else, convert the existing values to string and add to the respective strings
// \/ Save the presence and separation files


// get the metablock and block maps ==> probably merge the function into above function to avoid stack/heap copies of map lists
// load the pair tuples within a partition :: maintain a list of lists so that the pair tuples are separated by the partition
//      note that the file containing the pair tuples is NOT to be deleted here 
//      the output files will only contain partitioned values (1,0,2::distant presence) or separation values
//      another c++ program will take these files and the pair tuple partition file to create an intermediate text output :: delete the files then
// Load the pair presence file (partition)
// check if the current pair would exists in the current assembly (if its reported presence is 0)
//      check   if both features exists,
//              if so check their contig nos <== requires the contig ends vectors to be loaded
//              if contig nos are different, check if contig orientations are available <== required oriented links to be loaded
//              if relative orientations are consistent with the expected orientations from the pair
//      if these conditions are satisfied, set the update map entry for the pair: <..., current_assembly_distant_separation, .... > (else 0)
// save the loaded and updated presence combination in a map pointing to list of pair<bool,bool>(is_near, is_far)
//      (1,0): existing pair, (0,1): updated pair ==> output would show 2, (0,0): missing pair


void load_features_and_update(  string filtered_feature_dir, string paired_feature_dir, string start_assembly_idx_str, string end_assembly_idx_str,
                                short partition_count, long max_inter_feature_separation, string contig_len_dir, bool use_oriented_links,
                                string oriented_links_dir){
    cout<<"load_features_and_update started: "<<start_assembly_idx_str<<" "<<end_assembly_idx_str<<" "<<partition_count<<"\n";

    ulli partitioned_assembly_feature_count, block_idx, pair_count, idx, distant_pairs_count = 0, new_update_count = 0;//, feature_no;
    short feature_tuple_string_length, partition_idx;

    tup_uubb current_coords;
    tup_uub current_feature_coords;
    
    string feature_tuple_string, filename, feature_1_tuple_string, feature_2_tuple_string, opstr, opstr_presence;
    string filename_prefix = filtered_feature_dir + start_assembly_idx_str + "_" + end_assembly_idx_str + "_";

    unsigned long assembly_idx, start_assembly_idx = stoul(start_assembly_idx_str), end_assembly_idx = stoul(end_assembly_idx_str);

    list< map<string, tup_uub > > metablock_coords_map_list;
    list< map<ulli, tup_uub > > block_coords_map_list;

    list< map<string, tup_uub > >::iterator it_metablock_coords_map_list;
    list< map<ulli, tup_uub > >::iterator it_block_coords_map_list;

    map<string, tup_uub > metablock_coords_map;
    map<ulli, tup_uub > block_coords_map;

    map<string, tup_uub >::iterator it_metablock_coords_map;
    map<ulli, tup_uub >::iterator it_block_coords_map;

    list< vector<ulli> > contig_ends_vector_list;
    list< vector<ulli> >::iterator it_contig_ends_vector_list;

    vector<ulli> contig_ends_vector;

    list< map< puiui, tup_bbffs > > oriented_links_map_list;
    map< puiui, tup_bbffs > oriented_links_map;

    list< map< puiui, tup_bbffs > >::iterator it_ol_map_list;
    map< puiui, tup_bbffs >::iterator it_ol_map;

    list< list< tup_sbsbbb > > pair_tuples_outerlist;
    list< list< tup_sbsbbb > >::iterator it_pair_tuples_outerlist;
    list< tup_sbsbbb > pair_tuples_innerlist;
    list< tup_sbsbbb >::iterator it_pair_tuples_innerlist;

    list< list< list< pbb > > > pair_presence_pair_outerlist;
    list< list< list< pbb > > >::iterator it_pair_presence_pair_outerlist;
    list< list< pbb > > current_pair_presence_pair_list;
    list< list< pbb > >::iterator it_current_pair_presence_pair_list;
    list< pbb > pair_presence_pair_innerlist;
    list< pbb >::iterator it_pair_presence_pair_innerlist;

    list< list< list< long > > > pair_updated_separation_outerlist;
    list< list< list< long > > >::iterator it_pair_updated_separation_outerlist;
    list< list< long > > current_pair_updated_separation_list;
    list< list< long > >::iterator it_current_pair_updated_separation_list;
    list< long > pair_updated_separation_innerlist;
    list< long >::iterator it_pair_updated_separation_innerlist;

    list<string> outstrings_presence, outstrings_separation;
    list<string>::iterator it_outstrings_presence, it_outstrings_separation;

    tup_sbsbbb pair_tuple;
    bool is_1_metablock, is_2_metablock, relative_orientation_1, relative_orientation_2, presence, valid;
    pbb pair_presence;
    pbb update_presence_pair;
    unsigned int contig_idx_1, contig_idx_2;
    tup_uub coords_1, coords_2;
    long separation;
    puiui contig_pair;


    for(assembly_idx=start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++){
        metablock_coords_map_list.push_back(metablock_coords_map);
        block_coords_map_list.push_back(block_coords_map);

        // Load contig ends vectors and oriented links (if available)

        contig_ends_vector = get_contig_ends_vector(contig_len_dir + to_string(assembly_idx));
        contig_ends_vector_list.push_back(contig_ends_vector);

        if(use_oriented_links){
            oriented_links_map = load_oriented_links_map(oriented_links_dir + to_string(assembly_idx));
            oriented_links_map_list.push_back(oriented_links_map);
        }

        // outstrings_presence.push_back("");
        // outstrings_separation.push_back("");
    }

    ifstream inFile;
    ofstream outFile;


    // Load Individual features (metablocks and retained blocks)

    for(partition_idx = 0; partition_idx < partition_count; partition_idx++){

        filename =  filename_prefix + to_string(partition_idx);
        cout<<"\t"<<filename<<"\n";
        inFile.open(filename.c_str(), ios::binary);
        if(!inFile){
            cout<<"ALERT!!! ERROR IN LOADING "<<filename<<"!!!";
            return;
        }

        it_metablock_coords_map_list = metablock_coords_map_list.begin();
        it_block_coords_map_list = block_coords_map_list.begin();

        for(assembly_idx=start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++, it_metablock_coords_map_list++, it_block_coords_map_list++){
            // cout<<" "<<assembly_idx<<" ";
            inFile.read((char*) (&partitioned_assembly_feature_count), sizeof(partitioned_assembly_feature_count));
            cout<<"\t\t"<<partition_idx<<"\t"<<partitioned_assembly_feature_count;
            
            while(partitioned_assembly_feature_count--){
                inFile.read((char*) (&feature_tuple_string_length), sizeof(feature_tuple_string_length));
                feature_tuple_string.resize(feature_tuple_string_length);
                inFile.read( &feature_tuple_string[0], feature_tuple_string_length);
                inFile.read((char*) (&current_coords), sizeof(current_coords));

                current_feature_coords = make_tuple( get<0>(current_coords), get<1>(current_coords), get<2>(current_coords) );

                /*if(partitioned_assembly_feature_count%500==0){
                    cout<<" "<<feature_tuple_string<<" ("<< get<0>(current_feature_coords) <<","<< get<1>(current_feature_coords) <<",";
                    cout<< get<2>(current_feature_coords) <<") ";
                }*/

                if( get<3>(current_coords) ){
                    // metablock
                    metablock_coords_map[feature_tuple_string] = current_feature_coords;

                }
                else{
                    // block
                    block_idx = strtoulli(feature_tuple_string);
                    (*it_block_coords_map_list)[block_idx] = current_feature_coords;
                }
            }
            cout<<"  "<< (*it_metablock_coords_map_list).size()<<"  "<<(*it_block_coords_map_list).size();
        }
        cout<<"\n";

        inFile.close();
        remove(filename.c_str()); // To be removed here
    }

    cout<<"Individual features loaded for "<<start_assembly_idx<<" to "<<end_assembly_idx<<"\n";


    filename_prefix = paired_feature_dir + start_assembly_idx_str + "_" + end_assembly_idx_str + "_";

    for(partition_idx = 0; partition_idx < partition_count; partition_idx++){
        
        // Load the pair tuples within each partition

        pair_tuples_innerlist.clear();
        inFile.open((paired_feature_dir + "pair_tuples_" + to_string(partition_idx)).c_str(), ios::binary);
        if(!inFile){
            cout<<"ALERT!!! ERROR IN LOADING PAIR TUPLES FILE "<<filename<<"!!! "<<" : "<<pair_count<<'\n';
            continue;
            // return;
        }
        inFile.read((char*) (&pair_count), sizeof(pair_count));

        for(idx=0; idx<pair_count; idx++){
            // Feature 1
            inFile.read((char*)(&feature_tuple_string_length), sizeof(feature_tuple_string_length));
            feature_1_tuple_string.resize(feature_tuple_string_length);
            inFile.read( &feature_1_tuple_string[0], feature_tuple_string_length);
            inFile.read((char*)(&is_1_metablock), sizeof(is_1_metablock));

            // Feature 2
            inFile.read((char*)(&feature_tuple_string_length), sizeof(feature_tuple_string_length));
            feature_2_tuple_string.resize(feature_tuple_string_length);
            inFile.read( &feature_2_tuple_string[0], feature_tuple_string_length);
            inFile.read((char*)(&is_2_metablock), sizeof(is_2_metablock));

            inFile.read((char*)(&relative_orientation_1), sizeof(relative_orientation_1));
            inFile.read((char*)(&relative_orientation_2), sizeof(relative_orientation_2));

            pair_tuples_innerlist.push_back(make_tuple( feature_1_tuple_string, is_1_metablock, feature_2_tuple_string, is_2_metablock,
                                                        relative_orientation_1, relative_orientation_2));
        }

        inFile.close();

        // pair_tuples_outerlist.append(pair_tuples_innerlist);

        cout<<"\tPair tuples loaded from "<<paired_feature_dir<<"pair_tuples_"<<partition_idx<<" : "<<pair_tuples_innerlist.size()<<"\n";


        // Load pair presence file corresponding to last loaded partitioned pair tuples

        current_pair_presence_pair_list.clear();

        /*it_outstrings_presence = outstrings_presence.begin();
        it_outstrings_separation = outstrings_separation.begin();
        for(assembly_idx=start_assembly_idx; assembly_idx <= end_assembly_idx; assembly_idx++){
            (*it_outstrings_presence) = "";
            (*it_outstrings_separation) = "";

            it_outstrings_presence++;
            it_outstrings_separation++;
        }*/

        if(pair_count==0){
            ofstream { (filename_prefix + "_" + to_string(partition_idx) + "_updated_presence").c_str() };
            ofstream { (filename_prefix + "_" + to_string(partition_idx) + "_updated_separation").c_str() };
            continue;
        }

        filename = filename_prefix + to_string(partition_idx) + "_presence";
        cout<<"\t"<<filename<<"\n";
        inFile.open(filename.c_str(), ios::binary);
        if(!inFile){
            cout<<"ALERT!!! ERROR IN LOADING PARTITIONED PAIR PRESENCE FILE "<<filename<<"!!! "<<" : "<<pair_count<<'\n';
            ofstream { (filename_prefix + "_" + to_string(partition_idx) + "_updated_presence").c_str() };
            ofstream { (filename_prefix + "_" + to_string(partition_idx) + "_updated_separation").c_str() };
            continue;
            // return;
        }

        for(idx=0; idx<pair_count; idx++){
            pair_presence_pair_innerlist.clear();
            for(assembly_idx=start_assembly_idx; assembly_idx <= end_assembly_idx; assembly_idx++){
                inFile.read((char*)(&presence), sizeof(presence));
                pair_presence_pair_innerlist.push_back(pbb(presence, false));
                // inFile.read((char*)(&pair_presence), sizeof(pair_presence));
                // pair_presence_pair_innerlist.push_back(pair_presence);
            }
            current_pair_presence_pair_list.push_back(pair_presence_pair_innerlist);
        }
        inFile.close();

        // pair_presence_pair_outerlist.push_back(current_pair_presence_pair_list);

        remove(filename.c_str()); // To be removed here

        cout<<"\tExisting near pair presence loaded from "<<filename<<" : "<<current_pair_presence_pair_list.size()<<"\n";


        // Check if the pairs from the currently loaded partition exist as distant pairs in the current assembly partition

        it_pair_tuples_innerlist = pair_tuples_innerlist.begin();
        it_current_pair_presence_pair_list = current_pair_presence_pair_list.begin();
        current_pair_updated_separation_list.clear();
        opstr_presence = "";
        new_update_count = 0;

        while(it_pair_tuples_innerlist != pair_tuples_innerlist.end() && it_current_pair_presence_pair_list != current_pair_presence_pair_list.end()){

            it_pair_presence_pair_innerlist = (*it_current_pair_presence_pair_list).begin();
            it_metablock_coords_map_list = metablock_coords_map_list.begin();
            it_block_coords_map_list = block_coords_map_list.begin();
            it_contig_ends_vector_list = contig_ends_vector_list.begin();
            if(use_oriented_links)
                it_ol_map_list = oriented_links_map_list.begin();

            pair_updated_separation_innerlist.clear();

            // it_outstrings_presence = outstrings_presence.begin();
            if(start_assembly_idx == 0){
                opstr_presence += get<0>(*it_pair_tuples_innerlist) + "|";
                opstr_presence += get<2>(*it_pair_tuples_innerlist) + "|";
                opstr_presence += ( ( get<3>(*it_pair_tuples_innerlist) )?"1|":"0|" );
                opstr_presence += ( ( get<4>(*it_pair_tuples_innerlist) )?"1|,":"0|," );
            }

            for(assembly_idx = start_assembly_idx;
                    assembly_idx <= end_assembly_idx && it_pair_presence_pair_innerlist != (*it_current_pair_presence_pair_list).end();
                    assembly_idx++){

                update_presence_pair = *it_pair_presence_pair_innerlist;
                valid = false;

                if(!(update_presence_pair.first)){
                    // if near pair (existing) is absent, check for distant pair presence

                    // check if the features exist
                    // begin by checking if feature 1 exists:

                    feature_1_tuple_string = get<0>( *it_pair_tuples_innerlist );
                    is_1_metablock = get<1>( *it_pair_tuples_innerlist );

                    if(is_1_metablock){
                        it_metablock_coords_map = (*it_metablock_coords_map_list).find(feature_1_tuple_string);
                        valid = ( it_metablock_coords_map != (*it_metablock_coords_map_list).end() );
                    }
                    else{
                        block_idx = strtoulli(feature_1_tuple_string);
                        it_block_coords_map = (*it_block_coords_map_list).find(block_idx);
                        valid = ( it_block_coords_map != (*it_block_coords_map_list).end() );
                    }

                    if(valid){

                        coords_1 = ((is_1_metablock)?(it_metablock_coords_map->second):(it_block_coords_map->second));

                        // check if feature 2 exists:

                        feature_2_tuple_string = get<2>( *it_pair_tuples_innerlist );
                        is_2_metablock = get<3>( *it_pair_tuples_innerlist );

                        if(is_2_metablock){
                            it_metablock_coords_map = (*it_metablock_coords_map_list).find(feature_2_tuple_string);
                            valid = ( it_metablock_coords_map != (*it_metablock_coords_map_list).end() );
                        }
                        else{
                            block_idx = strtoulli(feature_2_tuple_string);
                            it_block_coords_map = (*it_block_coords_map_list).find(block_idx);
                            valid = ( it_block_coords_map != (*it_block_coords_map_list).end() );
                        }

                        if(valid){

                            coords_2 = ((is_2_metablock)?(it_metablock_coords_map->second):(it_block_coords_map->second));

                            // both constituent features exists -- now check if they match the expected relative orientations
                            // get required relative orientations

                            relative_orientation_1 = get<4>( *it_pair_tuples_innerlist );
                            relative_orientation_2 = get<5>( *it_pair_tuples_innerlist );

                            // get feature contigs

                            contig_idx_1 = get_contig_index_by_coordinate((*it_contig_ends_vector_list), get<0>(coords_1), 0,
                                                    (*it_contig_ends_vector_list).size() - 1 );
                            contig_idx_2 = get_contig_index_by_coordinate((*it_contig_ends_vector_list), get<0>(coords_2), 0,
                                                    (*it_contig_ends_vector_list).size() - 1 );

                            // check if both features are from the same contig:
                            if(contig_idx_1 == contig_idx_2){
                                // check if the feature relative orientations are satisfied with the expected relative orientations

                                if( get<0>(coords_1) < get<0>(coords_2) ){
                                    // features are positionally ordered; check for orientations
                                    valid = ( (relative_orientation_1 == get<2>(coords_1)) && (relative_orientation_2 == get<2>(coords_2)) );
                                    separation = get<0>(coords_2) - get<1>(coords_1) - 1; // start(2) - end(1) - 1
                                }
                                else{
                                    valid = ( (relative_orientation_1 != get<2>(coords_1)) && (relative_orientation_2 != get<2>(coords_2)) );
                                    separation = get<0>(coords_1) - get<1>(coords_2) - 1; // start(1) - end(2) - 1
                                }

                                // if yes:  *it_pair_presence_pair_innerlist = pbb(false, true);
                                if(valid && separation > max_inter_feature_separation){
                                    *it_pair_presence_pair_innerlist = pbb(false, true);
                                    // add to list only if the distant pair exists
                                    pair_updated_separation_innerlist.push_back(separation);
                                    distant_pairs_count++;

                                    cout<< "\t\t";
                                    cout<< get<0>(*it_pair_tuples_innerlist) << "|";
                                    cout<< get<2>(*it_pair_tuples_innerlist) << "|";
                                    cout<< ( ( get<3>(*it_pair_tuples_innerlist) )?"1|":"0|" );
                                    cout<< ( ( get<4>(*it_pair_tuples_innerlist) )?"1|,":"0|," );
                                    cout<<":"<<assembly_idx<<":(same contig)"<<contig_idx_1<<","<<contig_idx_2<<":"<<separation<<"\n";
                                }
                                else
                                    valid = false;
                            }
                            else if(use_oriented_links){
                                // if use_oriented_links
                                //      check if the contigs are oriented
                                //          check if the feature relative orientations are satisfied with the expected relative orientations

                                if(contig_idx_1 < contig_idx_2){
                                    contig_pair = puiui(contig_idx_1, contig_idx_2);
                                    it_ol_map = (*it_ol_map_list).find(contig_pair);
                                    if( it_ol_map != (*it_ol_map_list).end() ){
                                        // contig oriented link pair found
                                        // check if the feature relative orientations are satisfied:
                                        // if relative_orientation_1 == true:   B-------[--(1)-->]---->E(0) or E<----[<--(0)--]----B(1)
                                        // if relative_orientation_1 == false:  B-------[<-(0)---]---->E(0) or E<----[---(1)->]----B(1)
                                        // ==> ( relative_orientation_1 == (contig_orientation XOR feature orientation) ) 

                                        // valid = (relative_orientation_1 == (( get<2>(coords_1) && ~(get<0>( (*it_ol_map).second )) ) ||
                                        //                                     ( ~get<2>(coords_1) && (get<0>( (*it_ol_map).second )) ) ));
                                        valid = (relative_orientation_1 == (( get<2>(coords_1) && !(get<0>( (*it_ol_map).second )) ) ||
                                                                            ( !get<2>(coords_1) && (get<0>( (*it_ol_map).second )) ) ));

                                        // if relative_orientation_2 == true:   (1)B-------[--(1)-->]---->E or (0)E<----[<--(0)--]----B
                                        // if relative_orientation_2 == false:  (1)B-------[<-(0)---]---->E or (0)E<----[---(1)->]----B
                                        // ==> ( relative_orientation_2 == (contig_orientation XNOR feature orientation) )

                                        // valid = valid && (relative_orientation_2 == (( get<2>(coords_2) && (get<1>( (*it_ol_map).second )) ) ||
                                        //                                             ( ~get<2>(coords_2) && ~(get<1>( (*it_ol_map).second )) ) ));
                                        valid = valid && (relative_orientation_2 == (( get<2>(coords_2) && (get<1>( (*it_ol_map).second )) ) ||
                                                                                    ( !get<2>(coords_2) && !(get<1>( (*it_ol_map).second )) ) ));
                                    }
                                    else
                                        valid = false;
                                }
                                else{
                                    contig_pair = puiui(contig_idx_2, contig_idx_1);
                                    it_ol_map = (*it_ol_map_list).find(contig_pair);
                                    if( it_ol_map != (*it_ol_map_list).end() ){
                                        // contig oriented link pair found
                                        // check if the feature relative orientations are satisfied:
                                        // if relative_orientation_2 == true:   B-------[--(1)-->]---->E(0) or E<----[<--(0)--]----B(1)
                                        // if relative_orientation_2 == false:  B-------[<-(0)---]---->E(0) or E<----[---(1)->]----B(1)
                                        // ==> ( relative_orientation_2 == (contig_orientation XOR feature orientation) ) 

                                        // valid = (relative_orientation_2 == (( get<2>(coords_2) && ~(get<0>( (*it_ol_map).second )) ) ||
                                        //                                     ( ~get<2>(coords_2) && (get<0>( (*it_ol_map).second )) ) ));
                                        valid = (relative_orientation_2 == (( get<2>(coords_2) && !(get<0>( (*it_ol_map).second )) ) ||
                                                                            ( !get<2>(coords_2) && (get<0>( (*it_ol_map).second )) ) ));

                                        // if relative_orientation_1 == true:   (1)B-------[--(1)-->]---->E or (0)E<----[<--(0)--]----B
                                        // if relative_orientation_1 == false:  (1)B-------[<-(0)---]---->E or (0)E<----[---(1)->]----B
                                        // ==> ( relative_orientation_1 == (contig_orientation XNOR feature orientation) )

                                        // valid = valid && (relative_orientation_1 == (( get<2>(coords_1) && (get<1>( (*it_ol_map).second )) ) ||
                                        //                                             ( ~get<2>(coords_1) && ~(get<1>( (*it_ol_map).second )) ) ));
                                        valid = valid && (relative_orientation_1 == (( get<2>(coords_1) && (get<1>( (*it_ol_map).second )) ) ||
                                                                                    ( !get<2>(coords_1) && !(get<1>( (*it_ol_map).second )) ) ));
                                    }
                                    else
                                        valid = false;
                                }

                                // if yes:  *it_pair_presence_pair_innerlist = pbb(false, true);

                                if(valid){

                                    // distance from 1st feature contig
                                    if(get<0>( (*it_ol_map).second )){
                                        // distance to contig start
                                        if(contig_idx_1 == 0)
                                            separation = get<0>(coords_1) - 1;
                                        else
                                            separation = get<0>(coords_1) - 1 - (*it_contig_ends_vector_list)[contig_idx_1-1];
                                    }
                                    else{
                                        // distance to contig end
                                        separation = (*it_contig_ends_vector_list)[contig_idx_1] - get<1>(coords_1);
                                    }

                                    // distance from 2nd feature contig
                                    if(get<1>( (*it_ol_map).second )){
                                        // distance to contig start
                                        if(contig_idx_2 == 0)
                                            separation += get<0>(coords_2) - 1;
                                        else
                                            separation += get<0>(coords_2) - 1 - (*it_contig_ends_vector_list)[contig_idx_2-1];
                                    }
                                    else{
                                        // distance to contig end
                                        separation += (*it_contig_ends_vector_list)[contig_idx_2] - get<1>(coords_2);
                                    }

                                    // separation between the contigs
                                    separation += round(get<2>( (*it_ol_map).second )); // + get<3>( (*it_ol_map).second );

                                    valid = (separation > max_inter_feature_separation);

                                    if(valid){
                                        *it_pair_presence_pair_innerlist = pbb(false, true);

                                        pair_updated_separation_innerlist.push_back(separation);
                                        distant_pairs_count++;
                                        cout<< "\t\t";
                                        cout<< get<0>(*it_pair_tuples_innerlist) << "|";
                                        cout<< get<2>(*it_pair_tuples_innerlist) << "|";
                                        cout<< ( ( get<3>(*it_pair_tuples_innerlist) )?"1|":"0|" );
                                        cout<< ( ( get<4>(*it_pair_tuples_innerlist) )?"1|,":"0|," );
                                        cout<<":"<<assembly_idx<<":(via scaffolding)"<<contig_idx_1<<","<<contig_idx_2<<":"<<separation<<"\n";
                                    }
                                }
                            }
                        }
                    }

                    if(valid){
                        // distant pair found:
                        // (*it_outstrings_presence) += "2,";
                        opstr_presence += "2,";
                        new_update_count++;
                    }
                    else{
                        // pair does not exist:
                        // (*it_outstrings_presence) += "0,";
                        opstr_presence += "0,";
                    }
                }
                else{
                    // if(update_presence_pair.second){
                    //     // inital pair found was distant
                    //     opstr_presence += "2,";
                    // }
                    // else{
                    //     // near pair already found
                    //     // (*it_outstrings_presence) += "1,";
                    //     opstr_presence += "1,";
                    // }

                    // near pair already found
                    // (*it_outstrings_presence) += "1,";
                    opstr_presence += "1,";
                }


                it_pair_presence_pair_innerlist++;
                it_metablock_coords_map_list++;
                it_block_coords_map_list++;
                it_contig_ends_vector_list++;
                if(use_oriented_links)
                    it_ol_map_list++;

                // it_outstrings_presence++;
            }

            it_pair_tuples_innerlist++;
            it_current_pair_presence_pair_list++;
            opstr_presence.pop_back();
            opstr_presence += "\n";
            current_pair_updated_separation_list.push_back(pair_updated_separation_innerlist);
        }

        cout<<"\tCumulative distant pairs "<<partition_idx<<", "<<start_assembly_idx_str<<", "<<end_assembly_idx_str;
        cout<<" : "<<distant_pairs_count<<" "<<new_update_count<<"\n";

        //// Instead of loading all partitions and then saving the updates later, load a partition, process it, save it, and then go onto the next

        // pair_updated_separation_outerlist.push_back(current_pair_updated_separation_list);  


        // Save updated net presence
        /*opstr = "";
        for(it_outstrings_presence = outstrings_presence.begin(); it_outstrings_presence != outstrings_presence.end(); it_outstrings_presence++){
            opstr += (*it_outstrings_presence)
            (*it_outstrings_presence) = "";
            opstr.pop_back();
            opstr += "\n";
        }*/

        outFile.open((filename_prefix + to_string(partition_idx) + "_updated_presence").c_str(), ios::out);
        outFile << opstr_presence;
        outFile.close();

        cout<<"\tUpdated presence saved: "<<filename_prefix<<partition_idx<<"\n";


        // Load existing separation -- merge distant pair separation recorded and save

        it_pair_tuples_innerlist = pair_tuples_innerlist.begin();
        it_current_pair_presence_pair_list = current_pair_presence_pair_list.begin();
        // it_outstrings_separation = outstrings_separation.begin();
        it_current_pair_updated_separation_list = current_pair_updated_separation_list.begin();

        cout<<"\t\t"<<pair_tuples_innerlist.size()<<"\t\t"<<current_pair_presence_pair_list.size();
        cout<<"\t\t"<<outstrings_separation.size()<<"\t\t"<<current_pair_updated_separation_list.size()<<"\n";
        opstr = "";

        inFile.open((filename_prefix + to_string(partition_idx) + "_separation").c_str(), ios::binary);
        while(it_pair_tuples_innerlist != pair_tuples_innerlist.end() && it_current_pair_presence_pair_list != current_pair_presence_pair_list.end()){

            it_pair_presence_pair_innerlist = (*it_current_pair_presence_pair_list).begin();
            it_pair_updated_separation_innerlist = (*it_current_pair_updated_separation_list).begin();

            if(start_assembly_idx == 0){
                opstr += get<0>(*it_pair_tuples_innerlist) + "|";
                opstr += get<2>(*it_pair_tuples_innerlist) + "|";
                opstr += ( ( get<3>(*it_pair_tuples_innerlist) )?"1|":"0|" );
                opstr += ( ( get<4>(*it_pair_tuples_innerlist) )?"1|,":"0|," );
            }

            cout<< get<0>(*it_pair_tuples_innerlist) << "|";
            cout<< get<2>(*it_pair_tuples_innerlist) << "|";
            cout<< ( ( get<3>(*it_pair_tuples_innerlist) )?"1|":"0|" );
            cout<< ( ( get<4>(*it_pair_tuples_innerlist) )?"1|,":"0|," );

            for(assembly_idx = start_assembly_idx;
                    assembly_idx <= end_assembly_idx && it_pair_presence_pair_innerlist != (*it_current_pair_presence_pair_list).end();
                    assembly_idx++, it_pair_presence_pair_innerlist++){

                /*if( (*it_pair_presence_pair_innerlist).first ){
                    // Existing near pair
                    inFile.read((char*)(&separation), sizeof(separation));
                    (*it_outstrings_separation) += to_string(separation) + ",";
                }
                else{
                    if( (*it_pair_presence_pair_innerlist).second ){
                        // Detected distant pair
                        (*it_outstrings_separation) += to_string( *it_pair_updated_separation_innerlist ) + ",";
                        it_pair_updated_separation_innerlist++;
                    }
                    else{
                        // Pair does not exist
                        (*it_outstrings_separation) += "0,";
                    }
                }*/

                inFile.read((char*)(&separation), sizeof(separation));

                cout<<"\t"<<separation;

                if(!(*it_pair_presence_pair_innerlist).first && !(*it_pair_presence_pair_innerlist).second){
                    if(separation >0){
                        cout<<"ERROR!!! FEAUTRE_PAIR_ERROR!!! PAIR ABSENT SEPARATION PRESENT!!! "<<assembly_idx<<" ";
                        cout<<"("<<(*it_pair_presence_pair_innerlist).first<<","<<(*it_pair_presence_pair_innerlist).second<<") ";
                        cout<<separation<<"\n";
                    }
                }
                
                if( (*it_pair_presence_pair_innerlist).second ){
                    // Detected distant pair
                    cout<<"->";
                    cout<< *it_pair_updated_separation_innerlist;

                    opstr += to_string( *it_pair_updated_separation_innerlist ) + ",";
                    it_pair_updated_separation_innerlist++;
                }
                else{
                    // Pair either exists as a near pair or does not exists
                    opstr += to_string(separation) + ",";
                }
            }
            cout<<"\n";

            it_pair_tuples_innerlist++;
            it_current_pair_presence_pair_list++;
            opstr.pop_back();
            opstr += "\n";
            it_current_pair_updated_separation_list++;
        }
        inFile.close();

        /*opstr = "";
        for(it_outstrings_separation = outstrings_presence.begin(); it_outstrings_separation != outstrings_presence.end(); it_outstrings_separation++){
            opstr += (*it_outstrings_separation);
            (*it_outstrings_separation) = "";
            opstr.pop_back();
            opstr += "\n";
        }*/

        outFile.open((filename_prefix + to_string(partition_idx) + "_updated_separation").c_str(), ios::out);
        outFile << opstr;
        outFile.close();

        cout<<"\tUpdated separation saved: "<<filename_prefix<<partition_idx<<"\n\n";
    }

    cout<<"Total distant pairs found between "<<start_assembly_idx_str<<" and "<<end_assembly_idx_str<<" : "<<distant_pairs_count<<"\n";

    cout<<"load_features_and_update ended: "<<start_assembly_idx_str<<" "<<end_assembly_idx_str<<" "<<partition_count<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_distant_feature_pair_updating.cpp "<<argc<<"\n";
    string start_assembly_idx_str = argv[1];
    string end_assembly_idx_str = argv[2];
    short partition_count = static_cast<short>(stoi(argv[3]));
    long max_inter_feature_separation = stol(argv[4]);
    string filtered_feature_dir = argv[5];
    string paired_feature_dir = argv[6];
    string contig_len_dir = argv[7];
    bool use_oriented_links = strncmp(argv[8],"1",1)==0;
    string oriented_links_dir;

    if(use_oriented_links){
        oriented_links_dir = argv[9];
    }

    load_features_and_update(filtered_feature_dir, paired_feature_dir, start_assembly_idx_str, end_assembly_idx_str, partition_count,
                            max_inter_feature_separation, contig_len_dir, use_oriented_links, oriented_links_dir);

    return 0;
}
