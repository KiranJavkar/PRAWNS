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


void load_and_aggregate_features(string filtered_feature_dir, unsigned long assembly_count, unsigned long assemblies_per_partition,
                                string partition_idx_str){
    cout<<"load_and_aggregate_features started: "<<filtered_feature_dir<<" "<<assembly_count<<" "<<partition_idx_str<<"\n";

    ulli partitioned_assembly_feature_count, block_idx;//, feature_no;
    short feature_tuple_string_length, split_pos;
    
    unsigned long start_assembly_idx = 0, end_assembly_idx = assemblies_per_partition-1, assembly_idx_offset, missing_entries;
    // if(end_assembly_idx >= assembly_count)
    //     end_assembly_idx = assembly_count - 1;

    tup_uubb current_coords;
    tup_uub current_feature_coords, null_coords = make_tuple(0,0,false);
    
    string feature_tuple_string, filename, opstr_coords, opstr_bool, split_no_str, prev_tuple_string, opstr_block_idx;

    map<string, list< tup_uub > > metablock_tuple_idx_coords_map;
    map<ulli, list< tup_uub > > block_idx_coords_map;

    map<string, vector< bool > > metablock_presence_map;
    map<ulli, vector< bool > > block_presence_map;

    map<string, list< tup_uub > >::iterator it_metablock_coords_map;
    map<ulli, list< tup_uub > >::iterator it_block_coords_map;

    map<string, vector< bool > >::iterator it_metablock_presence_map;
    map<ulli, vector< bool > >::iterator it_block_presence_map;

    list<tup_uub>::iterator it_coords;
    vector<bool>::iterator it_presence;

    bool is_new_metablock;

    set<string> metablock_tuple_string_set;
    set<string>::iterator it_metablock_string_set;

    ifstream inFile;
    ofstream outFile;

    while(start_assembly_idx <= assembly_count){
        if(end_assembly_idx >= assembly_count)
            end_assembly_idx = assembly_count - 1;

        inFile.close();
        filename = filtered_feature_dir + to_string(start_assembly_idx) + "_" + to_string(end_assembly_idx) + "_" + partition_idx_str;
        cout<<"\t"<<filename<<"\n";
        inFile.open(filename.c_str(), ios::binary);
        if(!inFile){
            cout<<"ALERT!!! ERROR IN LOADING "<<filename<<"!!!";
            return;
        }

        for(assembly_idx_offset=start_assembly_idx; assembly_idx_offset<=end_assembly_idx; assembly_idx_offset++){
            // cout<<" "<<assembly_idx_offset<<" ";
            inFile.read((char*) (&partitioned_assembly_feature_count), sizeof(partitioned_assembly_feature_count));
            cout<<"\t\t"<<partitioned_assembly_feature_count;
            metablock_tuple_string_set.clear();
            
            while(partitioned_assembly_feature_count--){
                inFile.read((char*) (&feature_tuple_string_length), sizeof(feature_tuple_string_length));
                feature_tuple_string.resize(feature_tuple_string_length);
                inFile.read( &feature_tuple_string[0], feature_tuple_string_length);
                inFile.read((char*) (&current_coords), sizeof(current_coords));

                current_feature_coords = make_tuple( get<0>(current_coords), get<1>(current_coords), get<2>(current_coords) );

                if(partitioned_assembly_feature_count%500==0){
                    cout<<" "<<feature_tuple_string<<" ("<< get<0>(current_feature_coords) <<","<< get<1>(current_feature_coords) <<",";
                    cout<< get<2>(current_feature_coords) <<") ";
                }

                if( get<3>(current_coords) ){
                    // metablock

                    it_metablock_string_set = metablock_tuple_string_set.find(feature_tuple_string);
                    if(it_metablock_string_set != metablock_tuple_string_set.end()){
                        cout<<"ALERT!!! AGGREGATION_METABLOCK_STRING_REENCOUNTERED_ISSUE!!! "<<assembly_idx_offset<<" "<<feature_tuple_string<<"\n";
                    }
                    else
                        metablock_tuple_string_set.insert(feature_tuple_string);

                    it_metablock_coords_map = metablock_tuple_idx_coords_map.find(feature_tuple_string);
                    it_metablock_presence_map = metablock_presence_map.find(feature_tuple_string);

                    if(it_metablock_coords_map != metablock_tuple_idx_coords_map.end()){
                        // metablock tuple string encountered before: check if the list length is adjusted for missing values

                        if(it_metablock_presence_map == metablock_presence_map.end()){
                            cout<<"\nERROR!!! METABLOCK DETECTION ANAMOLY 1!!! "<<feature_tuple_string<<" "<<assembly_idx_offset<<"\n";
                            return;
                        }

                        if( (it_metablock_presence_map->second).size() >= assembly_idx_offset+1 ){
                            cout<<"\nERROR!!! METABLOCK PRESENCE ANAMOLY 2!!! "<<feature_tuple_string<<" "<<assembly_idx_offset<<" ";
                            cout<< (it_metablock_presence_map->second).size() << " " << prev_tuple_string<<"\n";
                            for(it_metablock_presence_map = metablock_presence_map.begin(); it_metablock_presence_map != metablock_presence_map.end(); it_metablock_presence_map++){
                                cout<<"\t"<< it_metablock_presence_map->first ;
                            }
                            cout<<'\n';
                            return;
                        }

                        // missing_entries = assembly_idx_offset - (it_metablock_coords_map->second).size();
                        missing_entries = assembly_idx_offset - (it_metablock_presence_map->second).size();
                        while(missing_entries--){
                            // (it_metablock_coords_map->second).push_back(null_coords);
                            (it_metablock_presence_map->second).push_back(false);
                        }

                        (it_metablock_coords_map->second).push_back(current_feature_coords);
                        (it_metablock_presence_map->second).push_back(true);
                    }
                    else{

                        if(it_metablock_presence_map != metablock_presence_map.end()){
                            cout<<"\nERROR!!! METABLOCK DETECTION ANAMOLY 2!!! "<<feature_tuple_string<<" "<<assembly_idx_offset<<"\n";
                            return;
                        }

                        list< tup_uub > new_feature_coords_list;
                        vector< bool > new_bool_vector;
                        missing_entries = assembly_idx_offset;
                        while(missing_entries--){
                            // new_feature_coords_list.push_back(null_coords);
                            new_bool_vector.push_back(false);
                        }

                        new_feature_coords_list.push_back(current_feature_coords);
                        new_bool_vector.push_back(true);

                        metablock_tuple_idx_coords_map[feature_tuple_string] = new_feature_coords_list;
                        metablock_presence_map[feature_tuple_string] = new_bool_vector;
                    }
                    prev_tuple_string = feature_tuple_string;
                }
                else{
                    // block
                    block_idx = strtoulli(feature_tuple_string);
                    it_block_coords_map = block_idx_coords_map.find(block_idx);
                    it_block_presence_map = block_presence_map.find(block_idx);

                    if(it_block_coords_map != block_idx_coords_map.end()){
                        // block_idx encountered before: check if the list length is adjusted for missing values
                        // missing_entries = assembly_idx_offset - (it_block_coords_map->second).size();
                        missing_entries = assembly_idx_offset - (it_block_presence_map->second).size();
                        while(missing_entries--){
                            // (it_block_coords_map->second).push_back(null_coords);
                            (it_block_presence_map->second).push_back(false);
                        }

                        (it_block_coords_map->second).push_back(current_feature_coords);
                        (it_block_presence_map->second).push_back(true);
                    }
                    else{
                        list< tup_uub > new_feature_coords_list;
                        vector< bool > new_bool_vector;
                        missing_entries = assembly_idx_offset;
                        while(missing_entries--){
                            // new_feature_coords_list.push_back(null_coords);
                            new_bool_vector.push_back(false);
                        }

                        new_feature_coords_list.push_back(current_feature_coords);
                        new_bool_vector.push_back(true);

                        block_idx_coords_map[block_idx] = new_feature_coords_list;
                        block_presence_map[block_idx] = new_bool_vector;
                    }
                }
            }

            cout<<"  "<<metablock_tuple_idx_coords_map.size()<<"  "<<block_idx_coords_map.size();
        }
        cout<<"\n";

        inFile.close();
        start_assembly_idx += assemblies_per_partition;
        end_assembly_idx += assemblies_per_partition;
        remove(filename.c_str());
    }

    cout<<"Metablocks in partition "<<partition_idx_str<<" : "<<metablock_tuple_idx_coords_map.size()<<"\t";
    cout<<"Blocks in partition "<<partition_idx_str<<" : "<<block_idx_coords_map.size()<<"\n";


    // Save partitioned feature matrix

    it_metablock_presence_map = metablock_presence_map.begin();
    for(it_metablock_coords_map = metablock_tuple_idx_coords_map.begin();
            it_metablock_coords_map != metablock_tuple_idx_coords_map.end() && it_metablock_presence_map != metablock_presence_map.end();
            it_metablock_coords_map++, it_metablock_presence_map++){

        split_pos = (it_metablock_coords_map->first).find("_");
        split_pos = (it_metablock_coords_map->first).find("_", split_pos+1);

        split_no_str = (it_metablock_coords_map->first).substr(split_pos+1);
        is_new_metablock = (split_no_str=="0");

        opstr_coords = it_metablock_coords_map->first + ",";
        if(is_new_metablock)
            opstr_bool = (it_metablock_coords_map->first).substr(0, split_pos) + ",";


        assembly_idx_offset = 0;
        // it_presence = (it_metablock_presence_map->second).begin();
        it_coords = (it_metablock_coords_map->second).begin();
        // for(it_coords = (it_metablock_coords_map->second).begin(); it_coords != (it_metablock_coords_map->second).end();
        //         it_coords++, assembly_idx_offset++){
        for(it_presence = (it_metablock_presence_map->second).begin(); it_presence != (it_metablock_presence_map->second).end();
                it_presence++, assembly_idx_offset++){

            if( *it_presence ){
                opstr_coords += to_string( get<0>(*it_coords) ) + "," + to_string( get<1>(*it_coords) ) + ",";
                opstr_coords += ( ( get<2>(*it_coords) )?"1,":"0," );

                if(is_new_metablock){
                    opstr_bool += "1,";
                    // opstr_bool += ( ( get<0>(*it_coords) == 0)?"0,":"1," );
                }

                it_coords++;
            }
            else{
                opstr_coords += "0,0,0,";
                opstr_bool += "0,";
            }
        }

        while(assembly_idx_offset++ < assembly_count){
            opstr_coords += "0,0,0,";
            opstr_bool += "0,";
            // if(is_new_metablock)
            //     opstr_bool += "0,";
        }

        opstr_coords.pop_back();
        opstr_coords += "\n";

        outFile.open((filtered_feature_dir + "metablock_coords_" + partition_idx_str).c_str(), ios::app);
        outFile << opstr_coords;
        outFile.close();

        opstr_bool.pop_back();
        opstr_bool += "\n";

        outFile.open((filtered_feature_dir + "metablock_presence_absence_" + partition_idx_str).c_str(), ios::app);
        outFile << opstr_bool;
        outFile.close();
    }

    cout<<"Metablocks saved\n";


    it_block_presence_map = block_presence_map.begin();
    opstr_block_idx = "";
    for(it_block_coords_map = block_idx_coords_map.begin();
            it_block_coords_map != block_idx_coords_map.end() && it_block_presence_map != block_presence_map.end();
            it_block_coords_map++, it_block_presence_map++){

        opstr_coords = to_string(it_block_coords_map->first) + ",";
        opstr_bool = to_string(it_block_coords_map->first) + ",";
        opstr_block_idx += to_string(it_block_coords_map->first) + "\n";

        assembly_idx_offset = 0;
        // it_presence = (it_block_presence_map->second).begin();
        it_coords = (it_block_coords_map->second).begin();
        // for(it_coords = (it_block_coords_map->second).begin(); it_coords != (it_block_coords_map->second).end();
        //         it_coords++, assembly_idx_offset++){
        for(it_presence = (it_block_presence_map->second).begin(); it_presence != (it_block_presence_map->second).end();
                it_presence++, assembly_idx_offset++){
            if( *it_presence){
                opstr_coords += to_string( get<0>(*it_coords) ) + "," + to_string( get<1>(*it_coords) ) + ",";
                opstr_coords += ( ( get<2>(*it_coords) )?"1,":"0," );
                opstr_bool += "1,";
                // opstr_bool += ( ( get<0>(*it_coords) == 0)?"0,":"1," );
            }
            else{
                opstr_coords += "0,0,0,";
                opstr_bool += "0,";
            }
        }

        while(assembly_idx_offset++ < assembly_count){
            opstr_coords += "0,0,0,";
            opstr_bool += "0,";
        }

        opstr_coords.pop_back();
        opstr_coords += "\n";

        outFile.open((filtered_feature_dir + "block_coords_" + partition_idx_str).c_str(), ios::app);
        outFile << opstr_coords;
        outFile.close();

        opstr_bool.pop_back();
        opstr_bool += "\n";

        outFile.open((filtered_feature_dir + "block_presence_absence_" + partition_idx_str).c_str(), ios::app);
        outFile << opstr_bool;
        outFile.close();
    }

    outFile.open((filtered_feature_dir + "filtered_block_idx_" + partition_idx_str).c_str(), ios::out);
    outFile << opstr_block_idx;
    outFile.close();

    cout<<"load_and_aggregate_features ended: "<<filtered_feature_dir<<" "<<assembly_count<<" "<<partition_idx_str<<"\n";
}


void load_aggregate_and_filter_paired_features( string paired_feature_dir, unsigned long assembly_count, unsigned long assemblies_per_partition,
                                                string partition_idx_str, unsigned long min_presence_count){
    cout<<"load_aggregate_and_filter_paired_features started: "<<paired_feature_dir<<" "<<assembly_count<<" "<<min_presence_count<<" "<<partition_idx_str<<"\n";

    ulli partitioned_assembly_feature_count, block_idx, feature_no, feature_count;
    short feature_1_tuple_string_length, feature_2_tuple_string_length, split_pos;
    
    unsigned long start_assembly_idx = 0, end_assembly_idx = assemblies_per_partition-1, assembly_idx_offset, missing_entries;
    long separation;

    tup_bbbbl current_paired_coords;
    
    string feature_string, feature_1_tuple_string, feature_2_tuple_string, filename, opstr_separation, opstr_bool;

    map<string, list< long > > paired_feature_separation_map;
    map<string, vector< bool > > paired_feature_presence_map;

    map<string, list< long > >::iterator it_paired_feature_map;
    map<string, vector< bool > >::iterator it_pair_presence_map;

    list<long>::iterator it_separation;
    vector<bool>::iterator it_presence;

    ifstream inFile;
    ofstream outFile;

    while(start_assembly_idx <= assembly_count){
        if(end_assembly_idx >= assembly_count)
            end_assembly_idx = assembly_count - 1;

        inFile.close();
        filename = paired_feature_dir + to_string(start_assembly_idx) + "_" + to_string(end_assembly_idx) + "_" + partition_idx_str;
        // cout<<"\t"<<filename<<"\n";
        inFile.open(filename.c_str(), ios::binary);
        if(!inFile){
            cout<<"ERROR IN LOADING "<<filename<<"!!!";
            return;
        }

        for(assembly_idx_offset=start_assembly_idx; assembly_idx_offset<=end_assembly_idx; assembly_idx_offset++){
            inFile.read((char*) (&partitioned_assembly_feature_count), sizeof(partitioned_assembly_feature_count));
            // cout<<"\t\t"<<partitioned_assembly_feature_count;

            while(partitioned_assembly_feature_count--){
                // feature 1
                inFile.read((char*)(&feature_1_tuple_string_length), sizeof(feature_1_tuple_string_length));
                feature_1_tuple_string.resize(feature_1_tuple_string_length);
                inFile.read( &feature_1_tuple_string[0], feature_1_tuple_string_length);

                feature_string = feature_1_tuple_string + "|";

                // feature 2
                inFile.read((char*)(&feature_2_tuple_string_length), sizeof(feature_2_tuple_string_length));
                feature_2_tuple_string.resize(feature_2_tuple_string_length);
                inFile.read( &feature_2_tuple_string[0], feature_2_tuple_string_length);

                feature_string += feature_2_tuple_string + "|";

                // feature types, relative orientations, and separation
                inFile.read((char*) (&current_paired_coords), sizeof(current_paired_coords));

                feature_string += ( ( get<2>(current_paired_coords) )?"1|":"0|" );
                feature_string += ( ( get<3>(current_paired_coords) )?"1|":"0|" );

                separation = get<4>(current_paired_coords);

                // if(partitioned_assembly_feature_count%500 == 0)
                //     cout<< " " << feature_string << " ";


                it_paired_feature_map = paired_feature_separation_map.find(feature_string);
                it_pair_presence_map = paired_feature_presence_map.find(feature_string);

                if(it_paired_feature_map != paired_feature_separation_map.end()){
                    // paired feature string encountered before: check if the list length is adjusted for missing values
                    // missing_entries = assembly_idx_offset - (it_paired_feature_map->second).size();
                    missing_entries = assembly_idx_offset - (it_pair_presence_map->second).size();
                    while(missing_entries--){
                        // (it_paired_feature_map->second).push_back(0);
                        (it_pair_presence_map->second).push_back(false);
                    }

                    (it_paired_feature_map->second).push_back(separation);
                    (it_pair_presence_map->second).push_back(true);
                }
                else{
                    list< long > new_separation_list;
                    vector< bool > new_bool_vector;
                    missing_entries = assembly_idx_offset;
                    while(missing_entries--){
                        // new_separation_list.push_back(0);
                        new_bool_vector.push_back(false);
                    }

                    new_separation_list.push_back(separation);
                    new_bool_vector.push_back(true);

                    paired_feature_separation_map[feature_string] = new_separation_list;
                    paired_feature_presence_map[feature_string] = new_bool_vector;
                }
            }
        }
        // cout<<"\n";

        inFile.close();
        start_assembly_idx += assemblies_per_partition;
        end_assembly_idx += assemblies_per_partition;
        remove(filename.c_str());
    }

    cout<<"Paired features (unfiltered) in partition "<<partition_idx_str<<" : "<<paired_feature_separation_map.size()<<"\t";


    // Filter and Save partitioned paired feature matrix

    it_pair_presence_map = paired_feature_presence_map.begin();

    for(it_paired_feature_map = paired_feature_separation_map.begin();
            it_paired_feature_map != paired_feature_separation_map.end() && it_pair_presence_map != paired_feature_presence_map.end();){

        // Filtering
        // if( (it_pair_presence_map->second) < min_presence_count){
        //     it_pair_presence_map = paired_feature_presence_map.erase(it_pair_presence_map);
        //     it_paired_feature_map = paired_feature_separation_map.erase(it_paired_feature_map);
        //     continue;
        // }


        opstr_separation = it_paired_feature_map->first + ",";
        opstr_bool = it_paired_feature_map->first + ",";

        assembly_idx_offset = 0;
        feature_count = 0;
        // it_presence = (it_pair_presence_map->second).begin();
        it_separation = (it_paired_feature_map->second).begin();

        // for(it_separation = (it_paired_feature_map->second).begin();
        //         it_separation != (it_paired_feature_map->second).end() && it_presence != (it_pair_presence_map->second).end();
        //         it_separation++, it_presence++, assembly_idx_offset++){

        for(it_presence = (it_pair_presence_map->second).begin();
                /*it_separation != (it_paired_feature_map->second).end() &&*/ it_presence != (it_pair_presence_map->second).end();
                /*it_separation++,*/ it_presence++, assembly_idx_offset++){

            // opstr_separation += to_string( *it_separation ) + ",";
            if( *it_presence ){
                if(it_separation == (it_paired_feature_map->second).end()){
                    cout<<"ALERT!!!ERROR IN PARSING PAIRED FEATURES LIST!!! "<< it_paired_feature_map->first <<"\n";
                    continue;
                }
                opstr_separation += to_string( *it_separation ) + ",";
                opstr_bool += "1,";
                feature_count++;
                it_separation++;
            }
            else{
                opstr_separation += "0,";
                opstr_bool += "0,";
            }
        }

        if(feature_count < min_presence_count){
            // Filter out lowly present pairs
            it_paired_feature_map = paired_feature_separation_map.erase(it_paired_feature_map);
            it_pair_presence_map = paired_feature_presence_map.erase(it_pair_presence_map);
            continue;
        }

        while(assembly_idx_offset++ < assembly_count){
            opstr_separation += "0,";
            opstr_bool += "0,";
        }

        opstr_separation.pop_back();
        opstr_separation += "\n";

        outFile.open((paired_feature_dir + "intrapair_separation_" + partition_idx_str).c_str(), ios::app);
        outFile << opstr_separation;
        outFile.close();

        opstr_bool.pop_back();
        opstr_bool += "\n";

        outFile.open((paired_feature_dir + "pair_presence_absence_" + partition_idx_str).c_str(), ios::app);
        outFile << opstr_bool;
        outFile.close();

        it_paired_feature_map++;
        it_pair_presence_map++;
    }

    cout<<"Paired features (filtered) in partition "<<partition_idx_str<<" : "<<paired_feature_separation_map.size()<<"\n";

    cout<<"load_aggregate_and_filter_paired_features ended: "<<paired_feature_dir<<" "<<assembly_count<<" "<<min_presence_count<<" "<<partition_idx_str<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_feature_aggregation_and_filtering.cpp "<<argc<<"\n";
    string partition_idx_str = argv[1];
    unsigned long assembly_count = stoul(argv[2]);
    unsigned long assemblies_per_partition = stoul(argv[3]);
    string filtered_feature_dir = argv[4];
    string paired_feature_dir = argv[5];
    unsigned long min_presence_count = stoul(argv[6]);

    load_and_aggregate_features(filtered_feature_dir, assembly_count, assemblies_per_partition, partition_idx_str);

    load_aggregate_and_filter_paired_features(paired_feature_dir, assembly_count, assemblies_per_partition, partition_idx_str, min_presence_count);

    return 0;
}