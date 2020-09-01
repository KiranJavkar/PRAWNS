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

    if(inFile.fail()){
        cerr<<oll_filename<<" cannot be opened\n";
        return oriented_links_map;
    }

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


map< unsigned long, vector<tup_uub> > load_required_block_coordinates(list<ulli> merged_sorted_block_list, list<unsigned long> merged_sorted_assembly_list,
                                                            string assembly_blocks_dir){
    cout<<"load_required_block_coordinates started: "<<merged_sorted_block_list.size()<<" "<<merged_sorted_assembly_list.size()<<" ";
    cout<<assembly_blocks_dir<<"\n";
    map< unsigned long, vector<tup_uub> > assembly_block_coords_map;
    list<string> lines, lsplit;
    list<string>::iterator it_coords;
    list<ulli>::iterator it_blocks;
    list<unsigned long>::iterator it_assembly;
    ulli current_block_idx, diff, start_pos, end_pos;
    bool orientation;
    string filename;

    /*for(it_blocks = merged_sorted_block_list.begin(); it_blocks != merged_sorted_block_list.end(); it_blocks++)
        cout<< *it_blocks <<" ";
    cout<<"\n\n";

    for(it_assembly = merged_sorted_assembly_list.begin(); it_assembly != merged_sorted_assembly_list.end(); it_assembly++)
        cout<< *it_assembly << " ";
    cout<<"\n\n";*/

    for(it_assembly = merged_sorted_assembly_list.begin(); it_assembly != merged_sorted_assembly_list.end(); it_assembly++){
        // cout<< (*it_assembly) <<"  ";
        vector<tup_uub> current_required_coordinates;
        // current_required_coordinates.reserve(merged_sorted_block_list.size());

        filename = assembly_blocks_dir + to_string(*it_assembly);
        lines = get_lines(filename.c_str());
        lsplit = split(lines.front(), ',');
        it_coords = lsplit.begin();
        current_block_idx = 0;

        // short show = 0;

        for(it_blocks = merged_sorted_block_list.begin(); it_blocks != merged_sorted_block_list.end(); it_blocks++){
            // cout<<"\n"<< *it_blocks << " " << current_block_idx << " ";

            diff = ( (*it_blocks) - current_block_idx )*3;

            // cout<< diff << " ";

            if(diff>0)
                advance(it_coords, diff);

            // cout<< *it_coords << "\n";

            start_pos = strtoulli( *it_coords++ );
            end_pos = strtoulli( *it_coords++ );
            orientation = ( *it_coords++ == "1" );

            current_required_coordinates.push_back( make_tuple(start_pos, end_pos, orientation) );

            /*if(show < 10){
                cout<< "[" << *it_blocks << " " << current_block_idx << " : (" << start_pos << "," << end_pos << "," << (orientation?"1":"0") << "] ";
                show++;
            }*/

            current_block_idx = (*it_blocks) + 1; // +1 due to the 3 it_coords++ operations
        }
        // cout<< current_required_coordinates.size() << "\n\n";

        assembly_block_coords_map[ *it_assembly ] = current_required_coordinates;
    }

    /*short count = 0;
    for(auto it = assembly_block_coords_map[ merged_sorted_assembly_list.front() ].begin();
            it != assembly_block_coords_map[ merged_sorted_assembly_list.front() ].end() && count < 10; it++){
        if( get<0>( *it ) > 0){
            cout << "(" << get<0>( *it ) << "," << get<1>( *it ) << "," << (get<2>( *it )?"1":"0") << ") \n";
            count++;
        }
    }*/

    cout<<"load_required_block_coordinates ended\n";
    return assembly_block_coords_map;
}


list<tup_uub> get_sorted_assembly_component_blocks( list<ulli> component_block_idx_list, list<ulli> merged_sorted_block_list,
                                                    vector<tup_uub> assembly_block_coords_vector){
    list<ulli>::iterator it_merged_blocks, it_current_blocks;
    list<tup_uub> block_locator_list; // (block_start_pos, block_idx, block_orientation)
    it_merged_blocks = merged_sorted_block_list.begin();
    it_current_blocks = component_block_idx_list.begin();
    ulli pos = 0;

    while(it_current_blocks != component_block_idx_list.end()){
        while(*it_current_blocks != *it_merged_blocks){
            it_merged_blocks++;
            pos++;
        }
        block_locator_list.push_back( make_tuple( get<0>( assembly_block_coords_vector[pos] ), *it_current_blocks,
                                                    get<2>( assembly_block_coords_vector[pos] ) ) );
        it_current_blocks++;
    }
    block_locator_list.sort();
    // cout << "get_sorted_assembly_component_blocks : "<<block_locator_list.size()<<" "<<assembly_block_coords_vector.size()<<"\n";
    // short count = 0;
    // for(auto it = block_locator_list.begin(); it != block_locator_list.end() && count<5; it++, count++)
    //     cout<< "(" << get<0>(*it) << "," << get<1>(*it) << "," << ((get<2>(*it))?"1":"0") << ")\n";
    return block_locator_list;
}


void check_coords_vector(vector<tup_uub> assembly_block_coords_vector){
    cout<<"check_coords_vector: "<< assembly_block_coords_vector.size() << " ";
    short count = 0;
    for(auto it = assembly_block_coords_vector.begin(); it != assembly_block_coords_vector.end() && count <10 ; it++)
        if( (get<0>(*it)) > 0){
            cout<< "(" << (get<0>(*it)) << "," << (get<1>(*it)) << ") ";
            count++;
        }
    cout<<"\n";
}



list<tup_uuu> get_sorted_assembly_component_block_coords( list<ulli> sorted_block_idx_list, list<ulli> merged_sorted_block_list,
                                                    vector<tup_uub> assembly_block_coords_vector){
    list<ulli>::iterator it_merged_blocks, it_current_blocks;
    list<tup_uuu> block_locator_list; // (block_start_pos, block_idx, block_end_pos)
    it_merged_blocks = merged_sorted_block_list.begin();
    it_current_blocks = sorted_block_idx_list.begin();
    ulli pos = 0, block_idx, start_pos, end_pos;

    while(it_current_blocks != sorted_block_idx_list.end()){
        while(*it_current_blocks != *it_merged_blocks){
            it_merged_blocks++;
            pos++;
        }
        start_pos = get<0>( assembly_block_coords_vector[pos] );
        end_pos = get<1>( assembly_block_coords_vector[pos] );
        block_idx = *it_current_blocks;
        // cout<<"("<<start_pos<<","<<block_idx<<","<<end_pos<<") ";
        block_locator_list.push_back( make_tuple( start_pos, block_idx, end_pos) );
        it_current_blocks++;
    }
    // cout<<"\n";
    block_locator_list.sort();

    // short count = 0;
    // for(auto it = block_locator_list.begin(); it != block_locator_list.end() && count<5; it++, count++)
    //     cout<< "(" << get<0>(*it) << "," << get<1>(*it) << "," << ((get<2>(*it))?"1":"0") << ")\n";
    return block_locator_list;
}


/*void generate_and_update_block_pair_counts( list<tup_uub> block_locator_list, vector<ulli> contig_ends_vector,
                                            map< tup_uubb, unsigned long >& block_pair_map, map< ulli, unsigned long >& block_map,
                                            bool new_component){
    unsigned int current_contig_idx = 0;
    bool new_contig_block = true;
    tup_uubb block_pair_tuple;
    tup_uub previous_block;
    map< tup_uubb, unsigned long >::iterator it_pair_map;
    map< ulli, unsigned long >::iterator it_single_map;

    for(list<tup_uub>::iterator it = block_locator_list.begin(); it != block_locator_list.end(); it++){
        if(get<0>( *it )==0)
            continue;

        while( get<0>( *it ) > contig_ends_vector[current_contig_idx]){
            current_contig_idx++;
            new_contig_block = true;
        }


        it_single_map = block_map.find( get<1>( *it ) );
        if(new_component && it_single_map == block_map.end())
            block_map[ get<1>( *it ) ] = 1;
        else
            (it_single_map->second) = (it_single_map->second) + 1;


        if( !new_contig_block && get<0>( *it ) > 0 && get<0>( previous_block ) > 0){
            if( get<1>(previous_block) < get<1>( *it ) )
                block_pair_tuple = make_tuple(get<1>(previous_block), get<1>( *it ), get<2>(previous_block), get<2>( *it ));
            else
                block_pair_tuple = make_tuple(get<1>( *it ), get<1>( previous_block ), !get<2>( *it ), !get<2>( previous_block ));

            it_pair_map = block_pair_map.find(block_pair_tuple);

            // We just need to look for core blocks
            //      ==> new tuples are to be added only if tuples from first assembly for that component are being loaded
            if(new_component && it_pair_map == block_pair_map.end())
                block_pair_map[block_pair_tuple] = 1;
            else
                (it_pair_map->second) = (it_pair_map->second) + 1;
        }
        new_contig_block = false;
        previous_block = *it;
    }
    // short count = 0;
    // for(it_pair_map = block_pair_map.begin(); it_pair_map != block_pair_map.end() && count<5; it_pair_map++, count++){
    //     cout << "("<< get<0>(it_pair_map->first) << "," << get<1>(it_pair_map->first) << ",";
    //     cout << ((get<2>(it_pair_map->first))?"1":"0") << "," << ((get<3>(it_pair_map->first))?"1":"0") << ")";
    //     cout<<":"<<it_pair_map->second<<"\n";
    // }
}*/


void generate_and_update_block_pair_counts( list<tup_uub> block_locator_list, vector<ulli> contig_ends_vector,
                                            map< tup_uubb, unsigned long >& block_pair_map, bool new_component){
    unsigned int current_contig_idx = 0;
    bool new_contig_block = true;
    tup_uubb block_pair_tuple;
    tup_uub previous_block;
    map< tup_uubb, unsigned long >::iterator map_it;

    for(list<tup_uub>::iterator it = block_locator_list.begin(); it != block_locator_list.end(); it++){
        if(get<0>( *it )==0)
            continue;

        while( get<0>( *it ) > contig_ends_vector[current_contig_idx]){
            current_contig_idx++;
            new_contig_block = true;
        }

        if( !new_contig_block && get<0>( *it ) > 0 && get<0>( previous_block ) > 0){
            if( get<1>(previous_block) < get<1>( *it ) )
                block_pair_tuple = make_tuple(get<1>(previous_block), get<1>( *it ), get<2>(previous_block), get<2>( *it ));
            else
                block_pair_tuple = make_tuple(get<1>( *it ), get<1>( previous_block ), !get<2>( *it ), !get<2>( previous_block ));

            map_it = block_pair_map.find(block_pair_tuple);

            // We just need to look for core blocks
            //      ==> new tuples are to be added only if tuples from first assembly for that component are being loaded
            if(new_component /*&& map_it == block_pair_map.end()*/)
                block_pair_map[block_pair_tuple] = 1;
            else if (map_it != block_pair_map.end())
                (map_it->second) = (map_it->second) + 1;
        }
        new_contig_block = false;
        previous_block = *it;
    }
    // short count = 0;
    // for(map_it = block_pair_map.begin(); map_it != block_pair_map.end() && count<5; map_it++, count++){
    //     cout << "("<< get<0>(map_it->first) << "," << get<1>(map_it->first) << ",";
    //     cout << ((get<2>(map_it->first))?"1":"0") << "," << ((get<3>(map_it->first))?"1":"0") << ")";
    //     cout<<":"<<map_it->second<<"\n";
    // }
}


/*void retain_core_block_pairs(map< tup_uubb, unsigned long >& block_pair_map, map< ulli, unsigned long >& block_map,
                            unsigned long component_assembly_count){
    set<ulli> core_blocks_set;
    set<ulli>::iterator set_it;

    cout << "retain_core_block_pairs: "<<block_pair_map.size()<<" ";

    for(map< tup_uubb, unsigned long >::iterator it_pair_map = block_pair_map.begin(); it_pair_map != block_pair_map.end(); ){
        if(it_pair_map->second == component_assembly_count){
            core_blocks_set.insert( get<0>( it_pair_map->first) );
            core_blocks_set.insert( get<1>( it_pair_map->first) );
            it_pair_map++;
        }
        else
            it_pair_map = block_pair_map.erase(it_pair_map);
    }

    cout<<block_pair_map.size()<<" ";

    for(map< ulli, unsigned long >::iterator it_single_map = block_map.begin(); it_single_map != block_map.end(); ){
        set_it = core_blocks_set.find( it_single_map->first );
        if(it_single_map->second == component_assembly_count && set_it == core_blocks_set.end()){
            // The block is by itself a core block, but adjoining blocks in the vicinity do not form core tuple
            block_pair_map[ make_tuple( it_single_map->first, it_single_map->first, true, true) ] = component_assembly_count;
        }
        it_single_map = block_map.erase(it_single_map);
    }

    cout<<block_pair_map.size()<<"\n";
    // cout<<"\t("<<max_val<<")\t";
    // short count = 0;
    // for(map< tup_uubb, unsigned long >::iterator it_pair_map = block_pair_map.begin(); it_pair_map != block_pair_map.end() && count<5;
    //         it_pair_map++, count++){
    //     cout << "("<< get<0>(it_pair_map->first) << "," << get<1>(it_pair_map->first) << ",";
    //     cout << ((get<2>(it_pair_map->first))?"1":"0") << "," << ((get<3>(it_pair_map->first))?"1":"0") << ")";
    //     cout<<":"<<it_pair_map->second<<"\n";
    // }
}*/


void retain_core_block_pairs(map< tup_uubb, unsigned long >& block_pair_map, unsigned long component_assembly_count){
    // unsigned long max_val = 0;
    for(map< tup_uubb, unsigned long >::iterator map_it = block_pair_map.begin(); map_it != block_pair_map.end(); ){
        // max_val = (max_val > map_it->second)?max_val:map_it->second;
        if(map_it->second == component_assembly_count)
            map_it++;
        else
            map_it = block_pair_map.erase(map_it);
    }
    // cout<<"\t("<<max_val<<")\t";
    // short count = 0;
    // for(map< tup_uubb, unsigned long >::iterator map_it = block_pair_map.begin(); map_it != block_pair_map.end() && count<5;
    //         map_it++, count++){
    //     cout << "("<< get<0>(map_it->first) << "," << get<1>(map_it->first) << ",";
    //     cout << ((get<2>(map_it->first))?"1":"0") << "," << ((get<3>(map_it->first))?"1":"0") << ")";
    //     cout<<":"<<map_it->second<<"\n";
    // }
}


// This is not a struct for metablock
// This is just to locate the chains of core blocks:
//      some of these could then be merged - some of these could span across contig boundaries in some assemblies
// The collection of the subsequent merged and those that could not be merged chains forms the metablocks
struct Chain{
    pulliulli chain_idx;
    unsigned long reference_assembly_idx;
    pulliulli ends_blocks;
    list< pulliulli > ends_coords_list;
    list< bool > orientations_list;
};


list<ulli> get_core_block_list(map< tup_uubb, unsigned long > block_pair_map){
    cout<<"get_core_block_list started: "<<block_pair_map.size()<<"\n";
    set<ulli> core_blocks_set;
    for(map< tup_uubb, unsigned long >::iterator map_it = block_pair_map.begin(); map_it != block_pair_map.end(); map_it++){
        core_blocks_set.insert( get<0>( map_it->first) );
        core_blocks_set.insert( get<1>( map_it->first) );
    }
    cout<<core_blocks_set.size()<<" ";

    list<ulli> core_blocks_list(core_blocks_set.begin(), core_blocks_set.end());
    // list<ulli> core_blocks_list;
    // for(set<ulli>::iterator it = core_blocks_set.begin(); it != core_blocks_set.end(); it++){
    //     core_blocks_list.push_back( *it );
    // }
    core_blocks_list.sort();
    cout<<core_blocks_list.size()<<" ";
    return core_blocks_list;
}


// Had a bug: it is possible that in some genome, the chaining of blocks from the reference genome, is now split across 2 contigs
// The other contig may still be placed next to the previous one, but in reverse orientation
// If only the end coords and blocks from the reference genome chaining is used, this causes issues
// Maintain a set of just the block indices tuples (orientations can be discarded now as they are checked and verified before)
// Get the location sorted list of blocks, just as before
// For each current block, check if a core block tuple is formed with its preceding block - using this block tuple set
// If yes, extend the current chain
list< Chain > merge_block_pairs_into_chains(map< tup_uubb, unsigned long > block_pair_map, list<ulli> merged_sorted_block_list,
                                            vector<tup_uub> assembly_block_coords_vector, ulli comp_idx, short max_neighbour_separation,
                                            unsigned long reference_assembly_idx, unsigned long entry_location ){
    // cout<<"merge_block_pairs_into_chains started: "<<block_pair_map.size()<<" "<<merged_sorted_block_list.size()<<" ";
    // cout<<assembly_block_coords_vector.size()<<" "<<comp_idx<<" "<<max_neighbour_separation<<"\n";
    list< Chain > component_chains_list;
    set<ulli> core_blocks_set;
    set<pulliulli> core_block_tuples_set;
    set<pulliulli>::iterator it_pair_presence;
    ulli previous_block_idx;

    for(map< tup_uubb, unsigned long >::iterator map_it = block_pair_map.begin(); map_it != block_pair_map.end(); map_it++){
        core_blocks_set.insert( get<0>( map_it->first) );
        core_blocks_set.insert( get<1>( map_it->first) );

        if(get<0>( map_it->first) < get<1>( map_it->first))
            core_block_tuples_set.insert( make_pair( get<0>( map_it->first) , get<1>( map_it->first) ) );
        else
            core_block_tuples_set.insert( make_pair( get<1>( map_it->first) , get<0>( map_it->first) ) );
    }
    // cout<<core_blocks_set.size()<<" ";

    list<ulli> core_blocks_list(core_blocks_set.begin(), core_blocks_set.end());
    core_blocks_list.sort();

    // cout<<core_blocks_list.size()<<" ";

    list<tup_uuu> location_sorted_blocks_list = get_sorted_assembly_component_block_coords(core_blocks_list, merged_sorted_block_list,
                                                                                            assembly_block_coords_vector);
    // cout<<location_sorted_blocks_list.size()<<"\n";

    Chain current_chain;
    ulli current_chain_no = 0;
    list<tup_uuu>::iterator it_block_loc = location_sorted_blocks_list.begin();
    list< pulliulli > reset_ends_coords_list;
    list<bool> default_orientation_list;
    for(unsigned long assembly_idx=0; assembly_idx<entry_location; assembly_idx++){
        reset_ends_coords_list.push_back( make_pair(0,0) );
        default_orientation_list.push_back(false);
    }
    default_orientation_list.push_back(true);

    current_chain.chain_idx = make_pair(comp_idx, current_chain_no);
    // current_chain.ends_blocks = make_pair( get<1>( *it_block_loc ), get<1>( *it_block_loc ));
    current_chain.ends_coords_list = reset_ends_coords_list;
    current_chain.orientations_list = default_orientation_list;
    current_chain.reference_assembly_idx = reference_assembly_idx;
    
    pulliulli ends_coords = make_pair( get<0>( *it_block_loc ), get<2>( *it_block_loc ));
    pulliulli ends_blocks = make_pair( get<1>( *it_block_loc ), get<1>( *it_block_loc ));
    unsigned long assembly_idx;

    previous_block_idx = get<1>( *it_block_loc );

    cout << (current_chain.chain_idx).first << "," << (current_chain.chain_idx).second << " : " << get<1>( *it_block_loc );
    if(ends_coords.first==0 || ends_coords.second==0)
        cout<<"\tNULL_CORE_BLOCK_ISSUE!!! " << reference_assembly_idx << "("<<ends_coords.first<<","<<ends_coords.second<<")\n";

    it_block_loc++;

    while(it_block_loc != location_sorted_blocks_list.end()){

        if(previous_block_idx < get<1>( *it_block_loc ))
            it_pair_presence = core_block_tuples_set.find( make_pair( previous_block_idx , get<1>( *it_block_loc ) ) );
        else
            it_pair_presence = core_block_tuples_set.find( make_pair( get<1>( *it_block_loc ) , previous_block_idx ) );

        if( ( it_pair_presence != core_block_tuples_set.end() ) && 
                ( ends_coords.second + max_neighbour_separation + 1 >= get<0>( *it_block_loc ) ) ){
            cout << "-" << get<1>( *it_block_loc );
            ends_blocks.second = get<1>( *it_block_loc );
            ends_coords.second = get<2>( *it_block_loc );
            if(ends_coords.first==0 || ends_coords.second==0)
                cout<<"\tNULL_CORE_BLOCK_ISSUE!!! " << reference_assembly_idx << "("<<ends_coords.first<<","<<ends_coords.second<<")\n";
        }
        else{
            current_chain.ends_blocks = ends_blocks;
            (current_chain.ends_coords_list).push_back( ends_coords );
            component_chains_list.push_back(current_chain);
            cout << " :: [" << (current_chain.ends_blocks).first << ",";
            cout << (current_chain.ends_blocks).second << "] => ";
            cout << " :: [" << ((current_chain.ends_coords_list).front()).first << ",";
            cout << ((current_chain.ends_coords_list).front()).second << "]\n";

            current_chain_no++;

            current_chain.chain_idx = make_pair(comp_idx, current_chain_no);
            // current_chain.ends_blocks = make_pair( get<1>( *it_block_loc ), get<1>( *it_block_loc ));
            ends_blocks = make_pair( get<1>( *it_block_loc ), get<1>( *it_block_loc ));
            current_chain.orientations_list = default_orientation_list;
            ends_coords = make_pair( get<0>( *it_block_loc ), get<2>( *it_block_loc ));
            // (current_chain.ends_coords_list).clear();
            current_chain.ends_coords_list = reset_ends_coords_list;
            cout << (current_chain.chain_idx).first << "," << (current_chain.chain_idx).second << " : " << get<1>( *it_block_loc );
            if(ends_coords.first==0 || ends_coords.second==0)
                cout<<"\tNULL_CORE_BLOCK_ISSUE!!! " << reference_assembly_idx << "("<<ends_coords.first<<","<<ends_coords.second<<")\n";
        }

        previous_block_idx = get<1>( *it_block_loc );

        it_block_loc++;
    }
    current_chain.ends_blocks = ends_blocks;
    (current_chain.ends_coords_list).push_back(ends_coords);
    component_chains_list.push_back(current_chain);
    // cout << " :: [" << (current_chain.ends_blocks).first << ",";
    // cout << (current_chain.ends_blocks).second << "] => ";
    // cout << " :: [" << ((current_chain.ends_coords_list).front()).first << ",";
    // cout << ((current_chain.ends_coords_list).front()).second << "]\n";

    cout<<"\nmerge_block_pairs_into_chains "<<comp_idx<<" : "<<component_chains_list.size()<<"\n";

    return component_chains_list;
}


void update_chain_coordinates_lists(list< Chain >& all_component_chains_list, list<ulli> chain_ends_block_list,
                                    list<ulli> merged_sorted_block_list, vector<tup_uub> assembly_block_coords_vector,
                                    unsigned long current_assembly_idx/*, unsigned long entry_location*/,
                                    vector< bool > component_presence_vector, ulli start_comp_idx){

    list<tup_uuu> location_sorted_blocks_list = get_sorted_assembly_component_block_coords(chain_ends_block_list, merged_sorted_block_list,
                                                                                            assembly_block_coords_vector);
    // cout<< location_sorted_blocks_list.size() << " ";

    map< ulli, pulliulli > block_coords_map;
    list<tup_uuu>::iterator it_block_loc;

    ulli diff;

    for(it_block_loc = location_sorted_blocks_list.begin(); it_block_loc != location_sorted_blocks_list.end(); it_block_loc++)
        block_coords_map[ get<1>( *it_block_loc) ] = make_pair( get<0>( *it_block_loc) , get<2>( *it_block_loc) );

    // cout << block_coords_map.size() << "\n";

    pulliulli block_1_coords, block_2_coords;

    for(list< Chain >::iterator it_chain = all_component_chains_list.begin(); it_chain != all_component_chains_list.end(); it_chain++){
        /*if( (*it_chain).reference_assembly_idx == current_assembly_idx ){
            if( (*it_chain).ends_coords_list.front().first >= (*it_chain).ends_coords_list.front().second ){
                cout << "CHAIN_ISSUE!!! : (" << (*it_chain).chain_idx.first << "," << (*it_chain).chain_idx.second << ") ";
                cout << (*it_chain).ends_coords_list.front().first << " " << (*it_chain).ends_coords_list.front().second << "\n";
            }
            continue;
        }

        else*/ 
        if((*it_chain).reference_assembly_idx < current_assembly_idx){

            if(component_presence_vector[ (*it_chain).chain_idx.first - start_comp_idx] == false){
                (*it_chain).ends_coords_list.push_back( make_pair(0, 0) );
                (*it_chain).orientations_list.push_back( false );
                continue;
            }

            block_1_coords = block_coords_map[ (*it_chain).ends_blocks.first ];
            block_2_coords = block_coords_map[ (*it_chain).ends_blocks.second ];
            // // Append only if the chain exists in the assembly
            // if(block_1_coords.first > 0 && block_2_coords.first > 0){
            //     if(block_1_coords.first < block_2_coords.first){
            //         // forward strand orientation
            //         (*it_chain).ends_coords_list.push_back( make_pair( block_1_coords.first, block_2_coords.second) );
            //         (*it_chain).orientations_list.push_back( true );
            //     }
            //     else{
            //         // reverse strand orientation
            //         (*it_chain).ends_coords_list.push_back( make_pair( block_2_coords.first, block_1_coords.second) );
            //         (*it_chain).orientations_list.push_back( false );
            //     }
            // }

            // This will result in all valid list entries to be appropriately placed
            // (lengths may vary depending on later assemblies that do not contain the corresponding blocks)
            /*diff = entry_location - (*it_chain).orientations_list.size();
            while(diff>0){
                (*it_chain).ends_coords_list.push_back( make_pair( 0, 0) );
                (*it_chain).orientations_list.push_back( false );
                diff--;
            }*/
            if(block_1_coords.first==0 || block_2_coords.first==0 || block_1_coords.second==0 || block_2_coords.second==0){
                cout<<"CHAIN_ENDS_BLOCKS_NULL_COORDS_ISSUE!!!: ("<<(*it_chain).chain_idx.first<<","<<(*it_chain).chain_idx.second<<") : ";
                cout<<current_assembly_idx<<" ("<<block_1_coords.first<<","<<block_1_coords.second<<") ";
                cout<<"("<<block_2_coords.first<<","<<block_2_coords.second<<")\n";
            }

            if(block_1_coords.first < block_2_coords.first){
                // forward strand orientation
                (*it_chain).ends_coords_list.push_back( make_pair( block_1_coords.first, block_2_coords.second) );
                (*it_chain).orientations_list.push_back( true );
            }
            else{
                // reverse strand orientation
                (*it_chain).ends_coords_list.push_back( make_pair( block_2_coords.first, block_1_coords.second) );
                (*it_chain).orientations_list.push_back( false );
            }
        }
    }
}


void filter_chain_coordinates_lists(Chain& current_chain, list<unsigned long> assembly_idx_innerlist,
                                    list<unsigned long> merged_sorted_assembly_list){
    // cout<<"filter_chain_coordinates_lists started: ("<<current_chain.chain_idx.first<<","<<current_chain.chain_idx.second<<") ";
    // cout<<current_chain.orientations_list.size()<<" "<<assembly_idx_innerlist.size()<<" ";
    // short count=0;
    // for(auto it = assembly_idx_innerlist.begin(); it != assembly_idx_innerlist.end() && count<10; it++, count++)
    //     cout<< *it << " ";
    // cout<<"\n";

    // unsigned long assembly_idx, assembly_count;
    list<unsigned long>::iterator it_assembly_innerlist = assembly_idx_innerlist.begin();
    list<unsigned long>::iterator it_merged_assembly = merged_sorted_assembly_list.begin();
    list< pulliulli >::iterator it_ends_coords = current_chain.ends_coords_list.begin();
    list< bool >::iterator it_orient = current_chain.orientations_list.begin();

    // assembly_idx = 0;
    // assembly_count = current_chain.orientations_list.size();

    /*while(it_assembly_innerlist != assembly_idx_innerlist.end()){
        if(assembly_idx != (*it_assembly_innerlist) ){
            it_ends_coords = current_chain.ends_coords_list.erase(it_ends_coords);
            it_orient = current_chain.orientations_list.erase(it_orient);
        }
        else{
            it_assembly_innerlist++;
            it_ends_coords++;
            it_orient++;
        }
        assembly_idx++;
    }
    while(it_orient != current_chain.orientations_list.end()){
        it_ends_coords = current_chain.ends_coords_list.erase(it_ends_coords);
        it_orient = current_chain.orientations_list.erase(it_orient);
    }*/

    list<pulliulli> filtered_ends_coords_list;
    list<bool> filtered_orientations_list;

    for(it_assembly_innerlist = assembly_idx_innerlist.begin(); it_assembly_innerlist != assembly_idx_innerlist.end();
            it_merged_assembly++, it_ends_coords++, it_orient++){
        if( (*it_assembly_innerlist) == (*it_merged_assembly) ){
            filtered_ends_coords_list.push_back( *it_ends_coords );
            filtered_orientations_list.push_back( *it_orient );
            if( it_ends_coords->first==0 || it_ends_coords->second==0 ){
                cout<<"CHAIN_COORDS_NULL_ISSUE!!! ("<<current_chain.chain_idx.first<<","<<current_chain.chain_idx.second<<") ";
                cout<<"<"<<current_chain.reference_assembly_idx << "> : ";
                cout<<(*it_assembly_innerlist)<<" "<<(it_ends_coords->first)<<" "<<(it_ends_coords->second)<<" "<<(*it_orient)<<"\n";
            }
            it_assembly_innerlist++;
        }
    }

    current_chain.ends_coords_list = filtered_ends_coords_list;
    current_chain.orientations_list = filtered_orientations_list;
    // cout<<"filter_chain_coordinates_lists ended: ("<<current_chain.chain_idx.first<<","<<current_chain.chain_idx.second<<") ";
    // cout<<current_chain.orientations_list.size()<<" "<<assembly_idx_innerlist.size()<<"\n";
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


struct Metablock{
    ulli comp_idx;
    bool is_contig_split;

    // "non-split" fields
    pulliulli ends_blocks;
    list< pulliulli > ends_coords_list;
    list< bool > orientations_list;

    // "split" fields
    list< pulliulli > contig_ends_blocks_list;
    list< list< pulliulli > > contig_ends_coords_list;
    list< list< bool > > contig_orientations_list;
};


Metablock init_chain_as_metablock(Chain current_chain){
    Metablock current_metablock;
    current_metablock.comp_idx = current_chain.chain_idx.first;
    current_metablock.is_contig_split = false;
    current_metablock.ends_blocks = current_chain.ends_blocks;
    current_metablock.ends_coords_list = current_chain.ends_coords_list;
    current_metablock.orientations_list = current_chain.orientations_list;

    return current_metablock;
}


list<unsigned int> get_chain_coords_contig_index_list(list<pulliulli> ends_coords_list, list< unsigned long> assembly_idx_innerlist,
                                                    map< unsigned long, vector<ulli> > mapped_contig_ends){
    list<unsigned int> contig_idx_list;
    vector<ulli> contig_ends_vector;
    list<pulliulli>::iterator it_ends_coords = ends_coords_list.begin();

    for(list<unsigned long>::iterator it_assembly_innerlist = assembly_idx_innerlist.begin(); it_assembly_innerlist != assembly_idx_innerlist.end();
            it_assembly_innerlist++, it_ends_coords++){
        contig_ends_vector = mapped_contig_ends[ *it_assembly_innerlist ];
        contig_idx_list.push_back( get_contig_index_by_coordinate( contig_ends_vector, (*it_ends_coords).first, 0, contig_ends_vector.size() - 1 ) );
    }

    return contig_idx_list;
}


/*list<ulli> get_chain_end_coord_list(Chain current_chain, bool is_preceding){
    list<ulli> end_coord_list;
    list<bool> orientations_list = current_chain.orientations_list;
    list<bool>::iterator it_orient = orientations_list.begin();

    for(list<pulliulli>::iterator it_ends_coords = current_chain.ends_coords_list.begin(); it_ends_coords != current_chain.ends_coords_list.end();
            it_ends_coords++, it_orient++){
        if(is_preceding){
            if( (*it_orient) )
                end_coord_list.push_back( (*it_ends_coords).second );
            else
                end_coord_list.push_back( (*it_ends_coords).first );
        }
        else{
            if( (*it_orient) )
                end_coord_list.push_back( (*it_ends_coords).first );
            else
                end_coord_list.push_back( (*it_ends_coords).second );
        }
    }

    return end_coord_list;
}*/


list<ulli> get_chain_end_coord_list(Chain current_chain, bool is_preceding, list<bool>assembly_contig_orientation_list){
    list<ulli> end_coord_list;
    list<bool> orientations_list = current_chain.orientations_list;
    list<bool>::iterator it_orient = orientations_list.begin();
    list<bool>::iterator it_contig_orient = assembly_contig_orientation_list.begin();

    for(list<pulliulli>::iterator it_ends_coords = current_chain.ends_coords_list.begin(); it_ends_coords != current_chain.ends_coords_list.end();
            it_ends_coords++, it_orient++, it_contig_orient++){
        if( !(is_preceding ^ (*it_contig_orient) ) ){
            if( (*it_orient) )
                end_coord_list.push_back( (*it_ends_coords).second );
            else
                end_coord_list.push_back( (*it_ends_coords).first );
        }
        else{
            if( (*it_orient) )
                end_coord_list.push_back( (*it_ends_coords).first );
            else
                end_coord_list.push_back( (*it_ends_coords).second );
        }
    }

    return end_coord_list;
}


// is_preceding ==> in forward orientation, if true: gives 'end' coordinate, else: gives 'start' coordinate; and vice versa
/*list<ulli> get_chain_end_coord_list(list<pulliulli> ends_coords_list, bool is_preceding, list<bool>assembly_contig_orientation_list){
    list<ulli> end_coord_list;
    // list<bool> orientations_list = current_chain.orientations_list;
    // list<bool>::iterator it_orient = orientations_list.begin();
    list<bool>::iterator it_orient = assembly_contig_orientation_list.begin();

    for(list<pulliulli>::iterator it_ends_coords = ends_coords_list.begin(); it_ends_coords != ends_coords_list.end();
            it_ends_coords++, it_orient++){
        if(is_preceding){
            if( (*it_orient) )
                end_coord_list.push_back( (*it_ends_coords).second );
            else
                end_coord_list.push_back( (*it_ends_coords).first );
        }
        else{
            if( (*it_orient) )
                end_coord_list.push_back( (*it_ends_coords).first );
            else
                end_coord_list.push_back( (*it_ends_coords).second );
        }
    }

    return end_coord_list;
}*/


list< Metablock > generate_metablocks(list< Chain > all_component_chains_list, list< list<unsigned long> > component_assembly_idx_outerlist,
                                    map< unsigned long, vector<ulli> > mapped_contig_ends, short max_neighbour_separation,
                                    bool use_oriented_links, map< unsigned long, map< puiui, tup_bbffs > > mapped_oriented_links_map,
                                    unsigned long iter_diff){
    cout<<"generate_metablocks started: "<< all_component_chains_list.size() << " " << component_assembly_idx_outerlist.size() << "\n";
    list< Metablock > metablock_list;

    Metablock current_metablock;
    bool current_contig_split, required_orientation, metablock_contig_orientation, chain_contig_orientation;
    bool chain_effective_orientation;

    list< Chain >::iterator it_chain = all_component_chains_list.begin();
    list< list< unsigned long > >::iterator it_assembly_outerlist = component_assembly_idx_outerlist.begin();
    advance(it_assembly_outerlist, iter_diff);

    list< unsigned int > metablock_end_contig_indices, chain_contig_indices;
    list< unsigned int >::iterator it_metablock_contig, it_chain_contig;
    list< ulli > metablock_start_coord_list, metablock_end_coord_list, chain_other_end_coord_list, chain_end_coord_list;
    list< ulli >::iterator it_metablock_start_coord, it_metablock_coord, it_chain_coord, it_chain_other_end_coord;
    list<bool> orientations_list, assembly_contig_orientation_list, metablock_last_orientations_list, orientations_innerlist;
    list<bool>::iterator it_orient, it_assembly_contig_orient, it_metablock_last_orient, it_metablock_orient;
    map< puiui, tup_bbffs > oriented_links_map;
    map< puiui, tup_bbffs >::iterator it_ol_map;
    tup_bbffs contig_orientation_tuple;
    long diff_metablock_contig, diff_chain_contig, diff_effective;
    // list< list< bool > > split_contig_orientations;
    // list< list< bool > >::iterator it_split_contig;
    unsigned long assembly_iter;

    list< list< pulliulli > >::iterator it_update_coords_outerlist;
    list< list< bool > >::iterator it_update_orientations_outerlist;
    list<pulliulli>::iterator it_ends_coords, it_metablock_ends_coords;
    list< pulliulli > ends_coords_innerlist;


    current_metablock = init_chain_as_metablock(*it_chain);
    metablock_end_contig_indices = get_chain_coords_contig_index_list(current_metablock.ends_coords_list, ( *it_assembly_outerlist ),
                                                                    mapped_contig_ends);

    for(assembly_iter = 0; assembly_iter < ( *it_assembly_outerlist ).size(); assembly_iter++)
        assembly_contig_orientation_list.push_back(true);
    it_assembly_contig_orient = assembly_contig_orientation_list.begin();

    // metablock_start_coord_list = get_chain_end_coord_list( *it_chain, !(*it_assembly_contig_orient) );
    // metablock_end_coord_list = get_chain_end_coord_list( *it_chain, (*it_assembly_contig_orient) );
    // metablock_start_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, false, (*it_chain).ends_coords_list.orientations_list);
    // metablock_end_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, true, (*it_chain).ends_coords_list.orientations_list);
    metablock_start_coord_list = get_chain_end_coord_list( *it_chain, false, assembly_contig_orientation_list);
    metablock_end_coord_list = get_chain_end_coord_list( *it_chain, true, assembly_contig_orientation_list);
    list< unsigned long >::iterator it_assembly_innerlist = ( *it_assembly_outerlist ).begin();
    current_contig_split = false;
    metablock_last_orientations_list = (*it_chain).orientations_list;
    it_metablock_last_orient = metablock_last_orientations_list.begin();
    cout<< current_metablock.comp_idx << " : " << (*it_chain).chain_idx.second << " ";
    // cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") ";
    short count = 0;
    auto it_bool_auto = (*it_chain).orientations_list.begin();
    for(auto it_auto = (*it_chain).ends_coords_list.begin(); it_auto != (*it_chain).ends_coords_list.end() && count<5;
            it_auto++, count++, it_bool_auto++)
        cout << "(" << (*it_auto).first << "," << (*it_auto).second << "," << ( (*it_bool_auto)?"1":"0" ) << ") ";

    it_chain++;
    bool is_chain_extension_valid;


    while(it_chain != all_component_chains_list.end()){
        if(metablock_last_orientations_list.size() != ( *it_assembly_outerlist ).size()){
            cout<<"\nLIST_SIZE_MISMATCH_ISSUE!!! : "<<( *it_assembly_outerlist ).size()<<" "<<metablock_last_orientations_list.size()<<" ";
            cout<<current_metablock.orientations_list.size()<<" "<<current_metablock.contig_orientations_list.size()<<"\n";
        }

        if( (*it_chain).chain_idx.first != current_metablock.comp_idx ){
            metablock_list.push_back(current_metablock);
            iter_diff = (*it_chain).chain_idx.first - current_metablock.comp_idx;
            // it_assembly_outerlist++;
            advance(it_assembly_outerlist, iter_diff);

            cout<< "\n";

            current_metablock = init_chain_as_metablock(*it_chain);
            metablock_end_contig_indices = get_chain_coords_contig_index_list(current_metablock.ends_coords_list, ( *it_assembly_outerlist ),
                                                                            mapped_contig_ends);
            // Once new metablock is to be formed, irrespective of use of oriented links, we assume the contig to be in the forward orientation.
            // This is to be done even if a part of the current contig was taken in flipped orientation, since the chains are sorted by start indices
            // for the 1st assembly containing the component, which are generated assuming the contigs to be in forward orientation
            // + the chains are sorted (by design) on comp_idx.

            assembly_contig_orientation_list.clear();
            for(assembly_iter = 0; assembly_iter < ( *it_assembly_outerlist ).size(); assembly_iter++)
                assembly_contig_orientation_list.push_back(true);
            it_assembly_contig_orient = assembly_contig_orientation_list.begin();

            // metablock_start_coord_list = get_chain_end_coord_list( *it_chain, !(*it_assembly_contig_orient) );
            // metablock_end_coord_list = get_chain_end_coord_list( *it_chain, (*it_assembly_contig_orient) );
            // metablock_start_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, false,
            //                                                         (*it_chain).ends_coords_list.orientations_list);
            // metablock_end_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, true,
            //                                                         (*it_chain).ends_coords_list.orientations_list);
            metablock_start_coord_list = get_chain_end_coord_list( *it_chain, false, assembly_contig_orientation_list);
            metablock_end_coord_list = get_chain_end_coord_list( *it_chain, true, assembly_contig_orientation_list);
            it_assembly_innerlist = ( *it_assembly_outerlist ).begin();
            current_contig_split = false;

            cout<<"\nNEW_COMPONENT\n";

            // Commented out as it anyways needs to be done at the end
            // metablock_last_orientations_list = (*it_chain).orientations_list;
            // it_metablock_last_orient = metablock_last_orientations_list.begin();

            cout<< current_metablock.comp_idx << " : " << (*it_chain).chain_idx.second << " ";
            // cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") ";
            count = 0;
            it_bool_auto = (*it_chain).orientations_list.begin();
            for(auto it_auto = (*it_chain).ends_coords_list.begin(); it_auto != (*it_chain).ends_coords_list.end() && count<5;
                    it_auto++, count++, it_bool_auto++)
                cout << "(" << (*it_auto).first << "," << (*it_auto).second << "," << ( (*it_bool_auto)?"1":"0" ) << ") ";

        }
        else{
            // check if end coord for the adjoining chain and current metablock
            //      are within the required max_neighbour_separation distance apart
            //      lie on the same contig
            //          if oriented links is available and some lie on separate contigs, check if contig orientation aligns for the pair formation
            //      if the requirements are fulfilled by all assemblies, update the metablock
            //          in case of split contig, append the new chain coords
            //          else, if the metablock hasn't been spanned across contigs so far (is_contig_split == false), update the 'non-split' fields
            //          else, update the last list entries from the 'split' lists fields

            it_assembly_contig_orient = assembly_contig_orientation_list.begin();
            it_metablock_last_orient = metablock_last_orientations_list.begin();

            chain_contig_indices = get_chain_coords_contig_index_list((*it_chain).ends_coords_list, ( *it_assembly_outerlist ),
                                                                            mapped_contig_ends);
            // chain_end_coord_list = get_chain_end_coord_list( *it_chain, !(*it_assembly_contig_orient) );
            // chain_end_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, false, (*it_chain).ends_coords_list.orientations_list);
            chain_end_coord_list = get_chain_end_coord_list( *it_chain, false, assembly_contig_orientation_list);

            // it_metablock_contig = metablock_end_contig_indices.begin();
            // it_chain_contig = chain_contig_indices.begin();
            // cout<<"\n";
            // while(it_metablock_contig != metablock_end_contig_indices.end() && it_chain_contig != chain_contig_indices.end()){
            //     cout << "{" << *it_metablock_contig << "," << *it_chain_contig << "} ";
            //     it_metablock_contig++;
            //     it_chain_contig++;
            // }
            // cout<<"\n";


            it_metablock_contig = metablock_end_contig_indices.begin();
            it_chain_contig = chain_contig_indices.begin();
            it_metablock_start_coord = metablock_start_coord_list.begin();
            it_metablock_coord = metablock_end_coord_list.begin();
            it_chain_coord = chain_end_coord_list.begin();

            is_chain_extension_valid = true;
            it_assembly_innerlist = ( *it_assembly_outerlist ).begin();
            // orientations_list = (*it_chain).orientations_list;
            // it_orient = orientations_list.begin();
            it_orient = (*it_chain).orientations_list.begin();

            if(use_oriented_links){
                // chain_other_end_coord_list = get_chain_end_coord_list( *it_chain, (*it_assembly_contig_orient) );
                // chain_other_end_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, true, 
                //                                                         (*it_chain).ends_coords_list.orientations_list);
                chain_other_end_coord_list = get_chain_end_coord_list( *it_chain, true, assembly_contig_orientation_list);
                it_chain_other_end_coord = chain_other_end_coord_list.begin();
            }

            while(is_chain_extension_valid && it_metablock_contig != metablock_end_contig_indices.end()){
                if( (*it_metablock_contig) == (*it_chain_contig) ){
                    if( /*(*it_orient)*/ (*it_metablock_last_orient) ){
                        if( (*it_metablock_start_coord) < (*it_chain_coord) ){
                            if ( ( (*it_metablock_coord) + max_neighbour_separation + 1 ) < (*it_chain_coord) ){
                                is_chain_extension_valid = false;
                                cout<<"\nLARGER_CHAIN_SEPARATION_1: "<< (*it_metablock_coord) << " " << (*it_chain_coord) << "\n";
                            }
                        }
                        else{
                            is_chain_extension_valid = false;
                            cout << "\nINVERTED_CHAIN_SEPARATION_1: "<< (*it_metablock_start_coord) << " " << (*it_metablock_coord) << " ";
                            cout << (*it_chain_coord) << "\n";
                        }
                    }
                    else{
                        if( (*it_metablock_start_coord) > (*it_chain_coord) ){
                            if ( ( (*it_chain_coord) + max_neighbour_separation + 1 ) < (*it_metablock_coord) ){
                                is_chain_extension_valid = false;
                                cout<<"\nLARGER_CHAIN_SEPARATION_2: "<< (*it_metablock_coord) << " " << (*it_chain_coord) << "\n";
                            }
                        }
                        else{
                            is_chain_extension_valid = false;
                            cout <<"\nINVERTED_CHAIN_SEPARATION_2: "<< (*it_metablock_start_coord) << " " << (*it_metablock_coord) << " ";
                            cout << (*it_chain_coord) << "\n";
                        }
                    }

                    // Biological relevance: genomic sequences have been identified which if present in a particular orientation in the upstream
                    // region of certain gene/s increase the expression levels which is, in turn, causal for a particular phenotype; and
                    // when these sequences are flipped around at the same location, decrease the expressions of the corresponding genes, thereby
                    // resulting in the other phenotype

                    if(( *it_assembly_innerlist ) == ( *it_assembly_outerlist ).front()){
                        required_orientation = !( (*it_metablock_last_orient) ^ (*it_orient) );
                    }
                    else{
                        // This is to check for the possibility of a small region, that by itself forms a chain, but can be flipped in different assemblies
                        if(is_chain_extension_valid && required_orientation != !( (*it_metablock_last_orient) ^ (*it_orient) )){
                            is_chain_extension_valid = false;
                            cout<<"\nINCONSISTENT_CHAIN_RELATIVE_ORIENTATION\n";
                        }
                    }

                }
                else{
                    if(!use_oriented_links){
                        is_chain_extension_valid = false;
                        cout<<"\nSEPARATE_CONTIGS_AND_NO_ORIENTED_LINKS\n";
                    }
                    else{
                        oriented_links_map = mapped_oriented_links_map[ (*it_assembly_innerlist) ];
                        it_ol_map = oriented_links_map.find( puiui( (*it_metablock_contig) , (*it_chain_contig) ) );

                        // If metablock is a part of the first of the 2 contigs from the oriented links row
                        //      -------[___METABLOCK___]---    ---[Chain]------
                        //          ==> the part of the contig that is upstream of the metablock should have been evaluated already
                        //          if metablock is oriented on the forward strand:
                        //              then the contig end linking the contigs should be 'E' (false)   B---[---METABLOCK--->]-----E
                        //          else the contig end linking the contigs should be 'B' (true)        E---[---METABLOCK--->]-----B
                        //      ==> for the oriented contig pairing to be valid, ( metablock orientation ^ contig end ) == false

                        if( it_ol_map != oriented_links_map.end() ){
                            contig_orientation_tuple = it_ol_map->second;
                            metablock_contig_orientation = get<0>(contig_orientation_tuple);
                            chain_contig_orientation = get<1>(contig_orientation_tuple);

                            // if ( metablock orientation ^ contig end ) == true ==> invalid extension
                            // else
                            //  check for orientation of the contig containing the chain
                            //  if the existing metablock and current chain pair has been checked before (i.e. this is not the first comparison made)
                            //      check if the contig orientations are consistent with the required relative orientation of the chain w.r.t. existing metablock
                            //  else, set the current identified relative orientation as the required orientation
                            //  check for the distances of the metablock and the chain from the corresponding contig ends (and subtract the mean distance)
                            //  verify if the effective separation between the metablock and the chain is within the max_neighbour_separation
                            //  if all these conditions are satisfied, set current_contig_split to true,
                            //      if the chain contig is flipped , i.e. corresponding contig end in oriented links is 'E' (false),
                            //          set (*it_assembly_contig_orient) = false

                            if( metablock_contig_orientation ^ *it_metablock_last_orient ){
                                is_chain_extension_valid = false;
                                cout<<"\nINVALID_CURRENT_METABLOCK_CONTIG_ORIENTATION_1\n";
                            }
                        }
                        else{
                            it_ol_map = oriented_links_map.find( puiui( (*it_chain_contig) , (*it_metablock_contig) ) );

                            // If the metablock is a part of the second of the 2 contigs from the oriented links row
                            //      ---------[Chain]---   ----[___METABLOCK___]------
                            //          ==> as before, the part of the metablock contig upstream of the current metablock has been evaluated
                            //              (right of metablock in the diagram)
                            //          if metablock is oriented on the forward strand:
                            //              then the contig end linking the contigs should again be 'E':  E---[<---METABLOCK---]-----B
                            //          else, the contig end linking the contigs should be 'B' (true):    B---[<---METABLOCK---]-----E
                            //       ==> as before, for the oriented contig pairing to be valid, ( metablock orientation ^ contig end ) == false

                            if( it_ol_map != oriented_links_map.end() ){
                                contig_orientation_tuple = it_ol_map->second;
                                metablock_contig_orientation = get<1>(contig_orientation_tuple);
                                chain_contig_orientation = get<0>(contig_orientation_tuple);

                                // this, therefore, follows identically to the previously mentioned checks for chain and its contig orientations

                                if( metablock_contig_orientation ^ *it_metablock_last_orient ){
                                    is_chain_extension_valid = false;
                                    cout<<"\nINVALID_CURRENT_METABLOCK_CONTIG_ORIENTATION_2\n";
                                }
                            }
                            else{
                                is_chain_extension_valid = false;
                                cout<<"\nCONTIGS_NOT_LINKED\n";
                            }
                        }

                        if(is_chain_extension_valid){
                            // as the checks to be done are the same irrespective of whether the metablock is a part of the first contig or second
                            // all the required checks can just be done here instead of writing the same steps twice above

                            chain_effective_orientation = !( chain_contig_orientation ^ (*it_orient) );

                            if(( *it_assembly_innerlist ) == ( *it_assembly_outerlist ).front()){
                                required_orientation = !( (*it_metablock_last_orient) ^ chain_effective_orientation );
                            }
                            else{
                                // This is to check for the possibility of a small region, that by itself forms a chain, but can be flipped in different assemblies
                                if(required_orientation != !( (*it_metablock_last_orient) ^ chain_effective_orientation )){
                                    is_chain_extension_valid = false;
                                    cout<<"\nINCONSISTENT_CHAIN_RELATIVE_ORIENTATION\n";
                                }
                            }

                            if(!is_chain_extension_valid)
                                break;

                            if( (*it_metablock_last_orient) ){
                                diff_metablock_contig = mapped_contig_ends[ (*it_assembly_innerlist) ][ (*it_metablock_contig) ] - (*it_metablock_coord);
                            }
                            else{
                                if( (*it_metablock_contig) > 0)
                                    diff_metablock_contig = mapped_contig_ends[ (*it_assembly_innerlist) ][ (*it_metablock_contig) - 1] - (*it_metablock_coord);
                                else
                                    diff_metablock_contig = (*it_metablock_coord) - 1;
                            }

                            if(chain_contig_orientation){
                                if( (*it_chain_contig) > 0){
                                    // diff_chain_contig = mapped_contig_ends[ (*it_assembly_innerlist) ][ (*it_chain_contig) - 1] - (*it_chain_other_end_coord);
                                    diff_chain_contig = mapped_contig_ends[ (*it_assembly_innerlist) ][ (*it_chain_contig) - 1] - (*it_chain_coord);
                                }
                                else{
                                    // diff_chain_contig = (*it_chain_other_end_coord) - 1;
                                    diff_chain_contig = (*it_chain_coord) - 1;
                                }
                            }
                            else{
                                // diff_chain_contig = mapped_contig_ends[ (*it_assembly_innerlist) ][ (*it_chain_contig) ] - (*it_chain_coord);
                                diff_chain_contig = mapped_contig_ends[ (*it_assembly_innerlist) ][ (*it_chain_contig) ] - (*it_chain_other_end_coord);
                            }

                            // Effective distance between linked contigs is taken to be sum of bases from metablock and chain ends to the respective contigs
                            //  plus the mean separation (would be negative in case of overlap) between the contigs
                            //  (haven't made use of the std deviation so far)

                            diff_effective = diff_metablock_contig + diff_chain_contig + (long)( get<2>(contig_orientation_tuple));

                            if(diff_effective > max_neighbour_separation){
                                is_chain_extension_valid = false;
                                cout<<"\nLARGE_INTERCONTIG_INTERCHAIN_SEPARATION: "<< (*it_metablock_coord) << " " << (*it_chain_coord);
                                cout<<" " << (*it_chain_other_end_coord) << " " << diff_metablock_contig << " " << diff_chain_contig << " ";
                                cout<<diff_effective<<" "<< ( get<2>(contig_orientation_tuple)) << " " << ( get<3>(contig_orientation_tuple)) <<"\n";
                            }
                            else{
                                if(chain_contig_orientation)
                                    *it_assembly_contig_orient = true; // previous value for current assembly contig orientation can be false
                                else
                                    *it_assembly_contig_orient = false;
                                current_contig_split = true;
                            }
                        }
                    }
                }

                // incremement all list pointers

                it_metablock_contig++;
                it_chain_contig++;
                it_metablock_start_coord++;
                it_metablock_coord++;
                it_chain_coord++;
                it_assembly_innerlist++;
                it_orient++;
                it_assembly_contig_orient++;
                it_metablock_last_orient++;

                if(use_oriented_links)
                    it_chain_other_end_coord++;
            }

            // if is_chain_extension_valid == true && current_contig_split == true
            //  set metablock's is_contig_split = True
            //  the ends blocks and end coords of the current chain are appended to the list the corresponding metablock lists
            //      if these lists were empty (is_contig_split was previously false),
            //          the non split fields are appended to these lists and then the relevant chain fields are appended
            //      if is_contig_split was previously false
            //          set the metablock's contig_orientations_list to have true as the first entry in each inner list
            //              (#inner lists == ( *it_assembly_outerlist ).size)
            //      append the assembly_contig_orientation_list entries to the corresponding inner lists of contig_orientations_list


            // if is_chain_extension_valid == true && current_contig_split == false
            //  if is_contig_split == false, update the required single field entries corresponding to the non-split fields
            //  else, update the required pair entry (requires taking into account the assembly contig orientation) from corresponding split lists

            if(is_chain_extension_valid){

                if(current_contig_split){
                    // Metablock spans across contig boundaries (can happen only if use_oriented_links == true)

                    cout << " <split>  ";

                    if(current_metablock.is_contig_split){

                        current_metablock.contig_ends_blocks_list.push_back( (*it_chain).ends_blocks );

                        it_ends_coords = ((*it_chain).ends_coords_list).begin();
                        for(it_update_coords_outerlist = (current_metablock.contig_ends_coords_list).begin();
                                it_update_coords_outerlist != (current_metablock.contig_ends_coords_list).end();
                                it_update_coords_outerlist++, it_ends_coords++){
                            ( *it_update_coords_outerlist ).push_back( (*it_ends_coords) );
                        }

                        it_assembly_contig_orient = assembly_contig_orientation_list.begin();
                        for(it_update_orientations_outerlist = (current_metablock.contig_orientations_list).begin();
                                it_update_orientations_outerlist != (current_metablock.contig_orientations_list).end();
                                it_update_orientations_outerlist++, it_assembly_contig_orient++){
                            // Metablock is on the reverse strand if *it_assembly_contig_orient == false
                            ( *it_update_orientations_outerlist ).push_back( (*it_assembly_contig_orient) );
                        }
                    }
                    else{

                        current_metablock.is_contig_split = true;
                        current_metablock.contig_ends_blocks_list.push_back( current_metablock.ends_blocks );
                        current_metablock.contig_ends_blocks_list.push_back( (*it_chain).ends_blocks );

                        // append the ends coordinates and orientations to the respective lists

                        it_metablock_ends_coords = (current_metablock.ends_coords_list).begin();
                        for(it_ends_coords = ((*it_chain).ends_coords_list).begin(); it_ends_coords != ((*it_chain).ends_coords_list).end();
                                it_ends_coords++, it_metablock_ends_coords++){

                            ends_coords_innerlist.clear();
                            ends_coords_innerlist.push_back( (*it_metablock_ends_coords) );
                            ends_coords_innerlist.push_back( (*it_ends_coords) );
                            (current_metablock.contig_ends_coords_list).push_back(ends_coords_innerlist);
                        }

                        
                        it_metablock_orient = (current_metablock.orientations_list).begin();
                        for(it_assembly_contig_orient = assembly_contig_orientation_list.begin();
                                it_assembly_contig_orient != assembly_contig_orientation_list.end();
                                it_assembly_contig_orient++, it_metablock_orient++){

                            orientations_innerlist.clear();
                            orientations_innerlist.push_back( (*it_metablock_orient) );
                            orientations_innerlist.push_back( (*it_assembly_contig_orient) );
                            (current_metablock.contig_orientations_list).push_back(orientations_innerlist);
                        }
                    }

                    metablock_start_coord_list = get_chain_end_coord_list( *it_chain, false, assembly_contig_orientation_list);

                }
                else{

                    if(current_metablock.is_contig_split){

                        bool are_ends_blocks_updated = false;

                        it_ends_coords = ((*it_chain).ends_coords_list).begin();
                        it_update_coords_outerlist = (current_metablock.contig_ends_coords_list).begin();

                        for(it_metablock_last_orient = metablock_last_orientations_list.begin();
                                it_metablock_last_orient != metablock_last_orientations_list.end();
                                it_metablock_last_orient++, it_update_coords_outerlist++, it_ends_coords++){

                            if(!are_ends_blocks_updated){
                                if( (*it_metablock_last_orient) ){
                                    // if current region of metablock is on forward strand,
                                    //      then the end block (right-most on forward strand) for current region is updated to
                                    //          the end block (right-most again) of the chain to be merged
                                    ((current_metablock.contig_ends_blocks_list).back()).second = ((*it_chain).ends_blocks).second;
                                }
                                else{
                                    ((current_metablock.contig_ends_blocks_list).back()).first = ((*it_chain).ends_blocks).first;
                                }
                                are_ends_blocks_updated = true;
                            }

                            if( (*it_metablock_last_orient) ){
                                (( *it_update_coords_outerlist ).back()).second = ( *it_ends_coords ).second;
                            }
                            else{
                                (( *it_update_coords_outerlist ).back()).first = ( *it_ends_coords ).first;
                            }
                        }
                    }
                    else{

                        bool are_ends_blocks_updated = false;

                        it_ends_coords = ((*it_chain).ends_coords_list).begin();
                        it_metablock_ends_coords = (current_metablock.ends_coords_list).begin();
                        cout<<"\n"<<(current_metablock.orientations_list).size() << " " << metablock_last_orientations_list.size() << "\n";

                        it_orient = (*it_chain).orientations_list.begin();

                        // it_metablock_last_orient = metablock_last_orientations_list.begin();
                        // for(it_metablock_orient = (current_metablock.orientations_list).begin();
                        //         it_metablock_orient != (current_metablock.orientations_list).end();
                        //         it_metablock_orient++, it_metablock_ends_coords++, it_ends_coords++, it_metablock_last_orient++){
                        it_metablock_orient = (current_metablock.orientations_list).begin();
                        for(it_metablock_last_orient = metablock_last_orientations_list.begin();
                                it_metablock_last_orient != metablock_last_orientations_list.end();
                                it_metablock_last_orient++, it_metablock_ends_coords++, it_ends_coords++, it_metablock_orient++, it_orient++){

                            if(!are_ends_blocks_updated){
                                // if( (*it_metablock_orient) ){
                                if( (*it_metablock_last_orient) ){
                                    // if metablock is on forward strand,
                                    //      then the end block (right-most on forward strand) for current region is updated to
                                    //          the end block (right-most again) of the chain to be merged
                                    (current_metablock.ends_blocks).second = ((*it_chain).ends_blocks).second;
                                }
                                else{
                                    (current_metablock.ends_blocks).first = ((*it_chain).ends_blocks).first;
                                }
                                are_ends_blocks_updated = true;
                            }

                            // if( (*it_metablock_orient) ){
                            ulli prev_coords_first = ( *it_metablock_ends_coords ).first;
                            ulli prev_coords_second = ( *it_metablock_ends_coords ).second;
                            if( (*it_metablock_last_orient) ){
                                ( *it_metablock_ends_coords ).second = ( *it_ends_coords ).second;
                            }
                            else{
                                ( *it_metablock_ends_coords ).first = ( *it_ends_coords ).first;
                            }

                            if( ( *it_metablock_ends_coords ).first > ( *it_metablock_ends_coords ).second ){
                                cout<<"ALERT!!!! ";
                                cout << current_metablock.comp_idx << " " ;
                                cout << prev_coords_first << " " << prev_coords_second << " " << (*it_metablock_contig) << " =>? ";
                                cout << ( *it_metablock_ends_coords ).first << " " << ( *it_metablock_ends_coords ).second << " ";
                                cout << (*it_chain_contig) << "\t\t" << (*it_metablock_orient) << " " << (*it_metablock_last_orient) << "\t\t";
                                cout << ( *it_ends_coords ).first << " " << ( *it_ends_coords ).second << " " << (*it_orient) << "\t\t";
                                cout << ( *it_assembly_outerlist ).size() << " " << metablock_last_orientations_list.size() << " ";
                                cout << (*it_chain).orientations_list.size() << "\n";
                            }
                        }
                    }
                }

                cout<< " -- " << (*it_chain).chain_idx.second << " ";
                // cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") ";
                count = 0;
                it_bool_auto = (*it_chain).orientations_list.begin();
                for(auto it_auto = (*it_chain).ends_coords_list.begin(); it_auto != (*it_chain).ends_coords_list.end() && count<5; it_auto++, count++, it_bool_auto++)
                    cout << "(" << (*it_auto).first << "," << (*it_auto).second << "," << ( (*it_bool_auto)?"1":"0" ) << ") ";


                metablock_end_contig_indices = chain_contig_indices;
                
                metablock_end_coord_list = get_chain_end_coord_list( *it_chain, true, assembly_contig_orientation_list);
            }
            else{
                // Metablock cannot be extended despite current chain's component being the same as the current metablock's component
                // Following steps are same as those for new metablock construction when chain's comp_idx is different
                // EXCEPT incrementation of it_assembly_outerlist as the chain is present in the same assemblies as that of the previous metablock

                metablock_list.push_back(current_metablock);
                // THIS IS THE ONLY DIFFERENCE
                // iter_diff = (*it_chain).chain_idx.first - current_metablock.comp_idx;
                // // it_assembly_outerlist++;
                // advance(it_assembly_outerlist, iter_diff);
                cout<< "\n";

                current_metablock = init_chain_as_metablock(*it_chain);
                metablock_end_contig_indices = get_chain_coords_contig_index_list(current_metablock.ends_coords_list, ( *it_assembly_outerlist ),
                                                                                mapped_contig_ends);
                // Once new metablock is to be formed, irrespective of use of oriented links, we assume the contig to be in the forward orientation.
                // This is to be done even if a part of the current contig was taken in flipped orientation, since the chains are sorted by start indices
                // for the 1st assembly containing the component, which are generated assuming the contigs to be in forward orientation
                // + the chains are sorted (by design) on comp_idx.

                assembly_contig_orientation_list.clear();
                for(assembly_iter = 0; assembly_iter < ( *it_assembly_outerlist ).size(); assembly_iter++)
                    assembly_contig_orientation_list.push_back(true);
                it_assembly_contig_orient = assembly_contig_orientation_list.begin();

                // metablock_start_coord_list = get_chain_end_coord_list( *it_chain, !(*it_assembly_contig_orient) );
                // metablock_end_coord_list = get_chain_end_coord_list( *it_chain, (*it_assembly_contig_orient) );
                // metablock_start_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, false,
                //                                                         (*it_chain).ends_coords_list.orientations_list);
                // metablock_end_coord_list = get_chain_end_coord_list( (*it_chain).ends_coords_list, true,
                //                                                         (*it_chain).ends_coords_list.orientations_list);
                metablock_start_coord_list = get_chain_end_coord_list( *it_chain, false, assembly_contig_orientation_list);
                metablock_end_coord_list = get_chain_end_coord_list( *it_chain, true, assembly_contig_orientation_list);
                it_assembly_innerlist = ( *it_assembly_outerlist ).begin();
                current_contig_split = false;

                // Commented out as it anyways needs to be done at the end
                // metablock_last_orientations_list = (*it_chain).orientations_list;
                // it_metablock_last_orient = metablock_last_orientations_list.begin();

                cout<< current_metablock.comp_idx << " : " << (*it_chain).chain_idx.second << " ";
                // cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") ";
                count = 0;
                it_bool_auto = (*it_chain).orientations_list.begin();
                for(auto it_auto = (*it_chain).ends_coords_list.begin(); it_auto != (*it_chain).ends_coords_list.end() && count<5; it_auto++, count++, it_bool_auto++)
                    cout << "(" << (*it_auto).first << "," << (*it_auto).second << "," << ( (*it_bool_auto)?"1":"0" ) << ") ";
            }


            // metablock is updated
        }

        current_contig_split = false; // reset contig split boolean for next chain
        // Orientations of last chain region form metablock's last merged region's orientation
        metablock_last_orientations_list = (*it_chain).orientations_list;
        it_metablock_last_orient = metablock_last_orientations_list.begin();
        it_chain++;
    }

    // append the last metablock constructed to the metablock_list
    metablock_list.push_back(current_metablock);
    cout<< "\n";

    return metablock_list;
}


/*
struct Metablock{
    ulli comp_idx;
    bool is_contig_split;

    // "non-split" fields
    pulliulli ends_blocks;
    list< pulliulli > ends_coords_list;
    list< bool > orientations_list;

    // "split" fields
    list< pulliulli > contig_ends_blocks_list;
    list< list< pulliulli > > contig_ends_coords_list;
    list< list< bool > > contig_orientations_list;
};
*/
void save_metablocks(list< Metablock > metablock_list, list< list<unsigned long> > component_assembly_idx_outerlist, ulli start_comp_idx,
                    unsigned long assembly_count, unsigned int assemblies_per_partition, string metablock_outfile_prefix){
    cout<<"save_metablocks started: "<<metablock_outfile_prefix<<" "<<metablock_list.size()<<"\n";

    list< list<unsigned long> >::iterator it_assembly_outerlist = component_assembly_idx_outerlist.begin();
    list<unsigned long>::iterator it_assembly_innerlist;
    ulli comp_idx = start_comp_idx, current_metablock_no = 0;
    unsigned long assembly_idx, partition_count = ceil(assembly_count*1.0/assemblies_per_partition), diff ;
    int current_comp_metablock = 0, split_no;
    bool append_idx_tup;

    list< tup_uib > metablock_idx_tuple_list;
    list< list< tup_uubu > > assembly_metablock_coords_outerlist;

    list< list< tup_uubu > > metablock_1st_assembly_coords_outerlist;

    list< pulliulli >::iterator it_ends_coords;
    list< bool >::iterator it_orient;
    list< list< pulliulli > >::iterator it_ends_coords_list;
    list< list< bool > >::iterator it_orient_list;
    list< list< tup_uubu > >::iterator it_assembly_metablock_coords_list;
    list< list< tup_uubu > >::iterator it_metablock_1st_assembly_coords_list;

    for(assembly_idx=0; assembly_idx<assembly_count; assembly_idx++){
        list< tup_uubu > current_metablock_coords_list, current_metablock_coords_list2;
        assembly_metablock_coords_outerlist.push_back(current_metablock_coords_list);
        metablock_1st_assembly_coords_outerlist.push_back(current_metablock_coords_list2);
    }

    for(list<Metablock>::iterator it_metablock = metablock_list.begin(); it_metablock != metablock_list.end(); it_metablock++, current_metablock_no++){

        while( (*it_metablock).comp_idx != comp_idx){
            comp_idx++;
            it_assembly_outerlist++;
            current_comp_metablock = 0;
        }

        it_assembly_innerlist = (*it_assembly_outerlist).begin();
        assembly_idx = 0;


        if( (*it_metablock).is_contig_split ){

            it_ends_coords_list = (*it_metablock).contig_ends_coords_list.begin();
            it_orient_list = (*it_metablock).contig_orientations_list.begin();
            it_assembly_metablock_coords_list = assembly_metablock_coords_outerlist.begin();
            it_metablock_1st_assembly_coords_list = metablock_1st_assembly_coords_outerlist.begin();

            if((*it_metablock).contig_ends_coords_list.size() != (*it_metablock).contig_orientations_list.size()){
                cout<<"SPLIT_CONTIG_METABLOCK_LIST_ISSUE_1!!! : ";
                cout<<(*it_metablock).contig_ends_coords_list.size()<<" "<<(*it_metablock).contig_orientations_list.size()<<"\n";
            }
            append_idx_tup = true;

            while(it_orient_list != (*it_metablock).contig_orientations_list.end()){
                // in case of a metablock split over multiple contigs
                //      there would be equivalent multiple repeated entries in the metablock_idx_tuple_list

                // metablock_idx_tuple_list.push_back( tup_uib( comp_idx, current_comp_metablock, (*it_metablock).is_contig_split ) );

                diff = (*it_assembly_innerlist) - assembly_idx;
                if(diff>0){
                    assembly_idx += diff;
                    advance(it_assembly_metablock_coords_list, diff);
                    if(append_idx_tup)
                        advance(it_metablock_1st_assembly_coords_list, diff);
                }
                // while(assembly_idx < (*it_assembly_innerlist)){
                //     assembly_idx++;
                //     it_assembly_metablock_coords_list++;
                // }

                it_ends_coords = (*it_ends_coords_list).begin();
                it_orient = (*it_orient_list).begin();
                split_no = 0;

                while(it_orient != (*it_orient_list).end()){

                    if((*it_ends_coords_list).size() != (*it_orient_list).size()){
                        cout<<"SPLIT_CONTIG_METABLOCK_LIST_ISSUE_2!!! : ";
                        cout<<(*it_ends_coords_list).size()<<" "<<(*it_orient_list).size()<<"\n";
                    }

                    if(append_idx_tup){
                        metablock_idx_tuple_list.push_back( tup_uib( comp_idx, current_comp_metablock, (*it_metablock).is_contig_split ) );
                        (*it_metablock_1st_assembly_coords_list).push_back( tup_uubu( (it_ends_coords->first), (it_ends_coords->second), (*it_orient),
                                                                                    current_metablock_no + split_no ) );
                    }

                    (*it_assembly_metablock_coords_list).push_back( tup_uubu( (it_ends_coords->first), (it_ends_coords->second), (*it_orient),
                                                                            current_metablock_no + split_no ) );
                    it_ends_coords++;
                    it_orient++;
                    split_no++;
                }

                it_assembly_innerlist++;
                it_ends_coords_list++;
                it_orient_list++;
                assembly_idx++;
                it_assembly_metablock_coords_list++;
                append_idx_tup = false;
            }

            current_metablock_no += (split_no-1);

        }
        else{
            it_ends_coords = (*it_metablock).ends_coords_list.begin();
            it_orient = (*it_metablock).orientations_list.begin();
            it_assembly_metablock_coords_list = assembly_metablock_coords_outerlist.begin();
            it_metablock_1st_assembly_coords_list = metablock_1st_assembly_coords_outerlist.begin();

            metablock_idx_tuple_list.push_back( tup_uib( comp_idx, current_comp_metablock, (*it_metablock).is_contig_split ) );
            append_idx_tup = true;

            while(it_orient != (*it_metablock).orientations_list.end()){

                diff = (*it_assembly_innerlist) - assembly_idx;
                if(diff>0){
                    assembly_idx += diff;
                    advance(it_assembly_metablock_coords_list, diff);
                    if(append_idx_tup)
                        advance(it_metablock_1st_assembly_coords_list, diff);
                }
                // while(assembly_idx < (*it_assembly_innerlist)){
                //     // outstring += "0,0,0,";
                //     assembly_idx++;
                //     it_assembly_metablock_coords_list++;
                // }

                // outstring += to_string(it_ends_coords->first) + "," + to_string(it_ends_coords->second) + ",";
                // outstring += (*it_orient)?"1,":"0,";

                (*it_assembly_metablock_coords_list).push_back( tup_uubu( (it_ends_coords->first), (it_ends_coords->second), (*it_orient),
                                                                        current_metablock_no ) );
                if(append_idx_tup){
                    (*it_metablock_1st_assembly_coords_list).push_back( tup_uubu( (it_ends_coords->first), (it_ends_coords->second), (*it_orient),
                                                                        current_metablock_no ) );
                    append_idx_tup = false;
                }

                it_assembly_innerlist++;
                it_ends_coords++;
                it_orient++;
                assembly_idx++;
                it_assembly_metablock_coords_list++;
            }

            // while(assembly_idx++ < assembly_count){
            //     // outstring += "0,0,0,";
            //     it_assembly_metablock_coords_list++;
            // }
        }

        // outstring.pop_back();
        // outstring += "\n";
        current_comp_metablock++;
    }

    cout<<assembly_metablock_coords_outerlist.size()<<" "<<metablock_idx_tuple_list.size()<<"\n";

    ofstream outFile;
    outFile.open((metablock_outfile_prefix + "metablock_indices").c_str(), ios::binary);
    ulli metablock_count = metablock_idx_tuple_list.size();
    outFile.write((char*) (&metablock_count), sizeof(metablock_count));

    for(list< tup_uib >::iterator it_idx = metablock_idx_tuple_list.begin(); it_idx != metablock_idx_tuple_list.end(); it_idx++){
        tup_uib &current_index = (*it_idx);
        outFile.write((char*) (&current_index), sizeof(current_index));
    }
    outFile.close();

    ulli start_assembly_idx = 0;
    it_assembly_metablock_coords_list = assembly_metablock_coords_outerlist.begin();
    list< tup_uubu >::iterator it_assembly_metablock_coords_innerlist;

    for(unsigned long partition_idx=0; partition_idx<partition_count; partition_idx++){

        outFile.open((metablock_outfile_prefix + to_string(partition_idx)).c_str(), ios::binary);

        for(assembly_idx = start_assembly_idx; assembly_idx<assembly_count && assembly_idx<(start_assembly_idx+assemblies_per_partition);
                assembly_idx++, it_assembly_metablock_coords_list++){
            metablock_count = (*it_assembly_metablock_coords_list).size();
            outFile.write((char*) (&metablock_count), sizeof(metablock_count));
            for(it_assembly_metablock_coords_innerlist = (*it_assembly_metablock_coords_list).begin();
                    it_assembly_metablock_coords_innerlist != (*it_assembly_metablock_coords_list).end(); it_assembly_metablock_coords_innerlist++){
                tup_uubu &current_coords = (*it_assembly_metablock_coords_innerlist);
                outFile.write((char*) (&current_coords), sizeof(current_coords));
            }
        }
        start_assembly_idx += assemblies_per_partition;

        outFile.close();
    }


    start_assembly_idx = 0;
    it_metablock_1st_assembly_coords_list = metablock_1st_assembly_coords_outerlist.begin();

    for(unsigned long partition_idx=0; partition_idx<partition_count; partition_idx++){

        outFile.open((metablock_outfile_prefix + "fo_" + to_string(partition_idx)).c_str(), ios::binary);

        for(assembly_idx = start_assembly_idx; assembly_idx<assembly_count && assembly_idx<(start_assembly_idx+assemblies_per_partition);
                assembly_idx++, it_metablock_1st_assembly_coords_list++){
            metablock_count = (*it_metablock_1st_assembly_coords_list).size();
            outFile.write((char*) (&metablock_count), sizeof(metablock_count));
            for(it_assembly_metablock_coords_innerlist = (*it_metablock_1st_assembly_coords_list).begin();
                    it_assembly_metablock_coords_innerlist != (*it_metablock_1st_assembly_coords_list).end(); it_assembly_metablock_coords_innerlist++){
                tup_uubu &current_coords = (*it_assembly_metablock_coords_innerlist);
                outFile.write((char*) (&current_coords), sizeof(current_coords));
            }
        }
        start_assembly_idx += assemblies_per_partition;

        outFile.close();
    }

    cout<<"save_metablocks ended: "<<metablock_outfile_prefix<<" "<<metablock_list.size()<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_metablock_constructor.cpp "<<argc<<"\n";
    string comp_group_filename = argv[1];
    string assembly_blocks_dir = argv[2];
    short max_neighbour_separation = static_cast<short>(stoi(argv[3]));
    string contig_len_dir  = argv[4];
    unsigned long assembly_count = stoul(argv[5]);
    unsigned int assemblies_per_partition = stoi(argv[6]);
    string metablock_outdir = argv[7];
    bool use_oriented_links = strncmp(argv[8],"1",1)==0;
    string oriented_links_dir;

    if(use_oriented_links){
        oriented_links_dir = argv[9];
    }


    list<string> lines = get_lines(comp_group_filename.c_str());
    // cout<<lines.size()<<"\n";
    list<string>::iterator line_it = lines.begin();
    ulli start_comp_idx = strtoulli(*line_it);
    line_it++;

    list< list<ulli> > component_block_idx_outerlist, merged_block_outerlist;
    list< list<unsigned long> > component_assembly_idx_outerlist, merged_assembly_outerlist;
    list< list<ulli> >::iterator it_block_outerlist;
    list< list<unsigned long> >::iterator it_assembly_outerlist;
    list<ulli> block_idx_innerlist;
    list<unsigned long> assembly_idx_innerlist;
    list<ulli>::iterator it_block_innerlist;
    list<unsigned long>::iterator it_assembly_innerlist;

    list<string> lsplit = split(comp_group_filename, '_');
    string substring = lsplit.back();
    lsplit = split(substring, '.');
    string group_no_str = lsplit.front();
    string metablock_outfile_prefix = metablock_outdir + group_no_str + "_";
    cout<<metablock_outfile_prefix<<"\n";

    list<string>::iterator it_idx;
    unsigned long assembly_idx;
    ulli comp_idx = start_comp_idx, merge_iter = 0;
    short count;


    // Load block and assembly indices

    while(line_it != lines.end()){
        block_idx_innerlist.clear();
        assembly_idx_innerlist.clear();
        // cout<< *line_it <<"\n";
        lsplit = split(*line_it, ',');

        for(it_idx = lsplit.begin(); it_idx != lsplit.end(); it_idx++)
            block_idx_innerlist.push_back(strtoulli(*it_idx));

        component_block_idx_outerlist.push_back(block_idx_innerlist);
        merged_block_outerlist.push_back(block_idx_innerlist);
        // short show = 0;
        // for(it_block_innerlist = block_idx_innerlist.begin(); it_block_innerlist != block_idx_innerlist.end() && show <10; it_block_innerlist++, show++)
        //     cout<< *it_block_innerlist <<" ";
        // cout<<"\n\n";
        line_it++;
        
        assembly_idx=0;
        lsplit = split(*line_it, ',');

        for(it_idx = lsplit.begin(); it_idx != lsplit.end(); it_idx++, assembly_idx++)
            if((*it_idx) == "1")
                assembly_idx_innerlist.push_back(assembly_idx);

        component_assembly_idx_outerlist.push_back(assembly_idx_innerlist);
        merged_assembly_outerlist.push_back(assembly_idx_innerlist);
        line_it++;
    }


    // Maintain a sorted list of all block indices across the corresponding assemblies from all components to be checked here
    // Allows loading of the corresponding block indices instead of redundant access to same files

    list<ulli> inner_list_1, inner_list_2, merged_sorted_block_list;
    list<unsigned long> inner_assembly_list_1, inner_assembly_list_2, merged_sorted_assembly_list;

    while(!merged_block_outerlist.empty()){
        inner_list_1.clear();
        inner_list_2.clear();
        inner_list_1 = merged_block_outerlist.front();
        merged_block_outerlist.pop_front();
        if(!merged_block_outerlist.empty()){
            inner_list_2 = merged_block_outerlist.front();
            merged_block_outerlist.pop_front();
            inner_list_1.merge(inner_list_2);
            merged_block_outerlist.push_back(inner_list_1);
        }
    }
    merged_sorted_block_list = inner_list_1;
    inner_list_1.clear();
    inner_list_2.clear();

    while(!merged_assembly_outerlist.empty()){
        inner_assembly_list_1.clear();
        inner_assembly_list_2.clear();
        inner_assembly_list_1 = merged_assembly_outerlist.front();
        merged_assembly_outerlist.pop_front();

        if(!merged_assembly_outerlist.empty()){
            inner_assembly_list_2 = merged_assembly_outerlist.front();
            merged_assembly_outerlist.pop_front();
            it_assembly_innerlist = inner_assembly_list_1.begin();
            while(it_assembly_innerlist!=inner_assembly_list_1.end() && !inner_assembly_list_2.empty()){
                if(*it_assembly_innerlist < inner_assembly_list_2.front())
                    it_assembly_innerlist++;
                else if(*it_assembly_innerlist == inner_assembly_list_2.front()){
                    it_assembly_innerlist++;
                    inner_assembly_list_2.pop_front();
                }
                else /*(*it_assembly_innerlist > inner_assembly_list_2.front())*/{
                    inner_assembly_list_1.insert(it_assembly_innerlist, inner_assembly_list_2.front());
                    inner_assembly_list_2.pop_front();
                }
            }
            while(!inner_assembly_list_2.empty()){
                inner_assembly_list_1.push_back(inner_assembly_list_2.front());
                inner_assembly_list_2.pop_front();
            }
            merged_assembly_outerlist.push_back(inner_assembly_list_1);
        }
    }
    merged_sorted_assembly_list = inner_assembly_list_1;
    cout<<merged_sorted_block_list.size()<<" "<<merged_sorted_assembly_list.size()<<"\n\n";
    inner_assembly_list_1.clear();
    inner_assembly_list_2.clear();


    vector< vector< bool > > assembly_component_presence_outervector;
    vector< bool > current_presence_vector;
    for(assembly_idx=0; assembly_idx < merged_sorted_assembly_list.size(); assembly_idx++)
        assembly_component_presence_outervector.push_back(current_presence_vector);

    vector< vector< bool > >::iterator it_component_presence_outervector;
    // ulli components_in_assembly;
    
    for(it_assembly_outerlist=component_assembly_idx_outerlist.begin(); it_assembly_outerlist!=component_assembly_idx_outerlist.end();
            it_assembly_outerlist++){

        list<unsigned long>::iterator it_merged_assembly;
        it_component_presence_outervector = assembly_component_presence_outervector.begin();

        it_assembly_innerlist = (*it_assembly_outerlist).begin();
        // components_in_assembly = 0;

        for(it_merged_assembly = merged_sorted_assembly_list.begin();
                it_merged_assembly != merged_sorted_assembly_list.end() && it_assembly_innerlist != (*it_assembly_outerlist).end();
                it_merged_assembly++, it_component_presence_outervector++){

            if( (*it_merged_assembly) == (*it_assembly_innerlist) ){
                (*it_component_presence_outervector).push_back(true);
                it_assembly_innerlist++;
                // components_in_assembly++;
            }
            else{
                (*it_component_presence_outervector).push_back(false);
            }
        }
        while(it_merged_assembly != merged_sorted_assembly_list.end()){
            (*it_component_presence_outervector).push_back(false);
            it_merged_assembly++;
            it_component_presence_outervector++;
        }
    }



    /*count = 25;
    for(it_block_innerlist = merged_sorted_block_list.begin(); it_block_innerlist != merged_sorted_block_list.end() && count>0;
            it_block_innerlist++, count--)
        cout<< *it_block_innerlist << " ";
    cout<<"\n";
    count = 25;
    for(it_assembly_innerlist = merged_sorted_assembly_list.begin(); it_assembly_innerlist != merged_sorted_assembly_list.end() && count>0;
            it_assembly_innerlist++, count--)
        cout<< *it_assembly_innerlist << " ";
    cout<<"\n";*/


    // Load the assemblywise coordinates of the required blocks from all components
    // If oriented links are available, load the contig length and oriented links maps for the corresponding assemblies

    map< unsigned long, vector<tup_uub> > assembly_block_coords_map = load_required_block_coordinates(merged_sorted_block_list,
                                                                                        merged_sorted_assembly_list, assembly_blocks_dir);
    cout<<assembly_block_coords_map.size()<<" "<<assembly_block_coords_map[merged_sorted_assembly_list.front()].size()<<"\n";

    vector<ulli> contig_ends_vector;
    map< puiui, tup_bbffs > oriented_links_map;
    map< unsigned long, vector<ulli> > mapped_contig_ends;
    map< unsigned long, map< puiui, tup_bbffs > > mapped_oriented_links_map;
    string assembly_idx_str;

    for(it_assembly_innerlist = merged_sorted_assembly_list.begin(); it_assembly_innerlist != merged_sorted_assembly_list.end();
            it_assembly_innerlist++){
        assembly_idx_str = to_string( *it_assembly_innerlist);
        contig_ends_vector = get_contig_ends_vector(contig_len_dir + assembly_idx_str);
        mapped_contig_ends[ *it_assembly_innerlist ] = contig_ends_vector;
        if(use_oriented_links){
            oriented_links_map = load_oriented_links_map(oriented_links_dir + assembly_idx_str);
            mapped_oriented_links_map[ *it_assembly_innerlist ] = oriented_links_map;
        }

        // cout << (*it_assembly_innerlist) << " : " << assembly_block_coords_map[ *it_assembly_innerlist ].size() << "\n";
    }
    cout<<"Contig ends (and Oriented links loaded): "<<mapped_contig_ends.size()<<" "<<mapped_oriented_links_map.size()<<"\n";

//     return 0;
// }

    // Metablock construction module

    //  Find and merge core blocks to form chains

    it_block_outerlist = component_block_idx_outerlist.begin();
    it_assembly_outerlist = component_assembly_idx_outerlist.begin();
    unsigned long reference_assembly_idx, component_assembly_count;

    map< tup_uubb, unsigned long > block_pair_map;
    // map< ulli, unsigned long > block_map;
    list<tup_uub> block_locator_list;
    bool new_component;
    comp_idx = start_comp_idx;
    list< Chain > component_chains_list, all_component_chains_list;
    list< Chain >::iterator it_chain;
    list<tup_uuu> location_sorted_blocks_list;
    list<ulli> chain_ends_block_list, core_blocks_list;
    // short count = 10;
    ulli discard_count = 0, cummulative_discard_len = 0, max_discard_len = 0, entry_location = 0;
    // list< list<ulli> > core_blocks_outerlist;
    // vector< bool > retained_comp_idx_vector;
    // list< pululli > comp_reference_assembly_idx_list;

    while(it_block_outerlist!=component_block_idx_outerlist.end() /*&& count-- > 0*/){
        // locate ordered core tuples froms blocks of the current component and form chains
        // merge chains (with oriented links if available) and form the metablocks' coordinates

        block_pair_map.clear();
        // block_map.clear();

        reference_assembly_idx = (*it_assembly_outerlist).front();
        component_assembly_count = (*it_assembly_outerlist).size();
        new_component = true;

        for(it_assembly_innerlist = (*it_assembly_outerlist).begin(); it_assembly_innerlist != (*it_assembly_outerlist).end();
                it_assembly_innerlist++){
            // if(new_component)
            //     reference_assembly_idx = *it_assembly_innerlist;
            // cout<< *it_assembly_innerlist <<" ";
            block_locator_list =  get_sorted_assembly_component_blocks( *it_block_outerlist, merged_sorted_block_list,
                                                                        assembly_block_coords_map[ *it_assembly_innerlist ]);
            // cout<< block_locator_list.size() <<" ";
            generate_and_update_block_pair_counts(block_locator_list, mapped_contig_ends[ *it_assembly_innerlist ], block_pair_map, new_component);
            // generate_and_update_block_pair_counts(block_locator_list, mapped_contig_ends[ *it_assembly_innerlist ], block_pair_map,
                                                // block_map, new_component);
            new_component = false;
            // cout<<block_pair_map.size()<<"\n";
        }

        // cout<<comp_idx<<" "<<component_assembly_count<<" "<<(*it_block_outerlist).size()<<"\t"<<block_pair_map.size()<<" => ";
        retain_core_block_pairs(block_pair_map, component_assembly_count);
        // retain_core_block_pairs(block_pair_map, block_map, component_assembly_count);
        // cout<<block_pair_map.size()<<"\n";

        if(block_pair_map.size() == 0){
            // cout<< comp_idx << " : " << "DISCARDED\n\n";
            discard_count++;
            cummulative_discard_len += (*it_block_outerlist).size();
            max_discard_len = (max_discard_len>(*it_block_outerlist).size())?max_discard_len:(*it_block_outerlist).size();
            // retained_comp_idx_vector.push_back(false);
        }
        else{
            cout<<"\n"<<block_pair_map.size()<<" "<<reference_assembly_idx<<" "<<(*it_assembly_outerlist).front()<<" ";
            cout<<component_assembly_count << " " << (*it_assembly_outerlist).size() << " ";
            cout<<assembly_block_coords_map[ (*it_assembly_outerlist).front() ].size()<<"\n";
            // component_chains_list = merge_block_pairs_into_chains(block_pair_map, merged_sorted_block_list,
            //                                             assembly_block_coords_map[reference_assembly_idx], comp_idx, max_neighbour_separation);
            entry_location = 0;
            for(it_assembly_innerlist = merged_sorted_assembly_list.begin();
                    it_assembly_innerlist != merged_sorted_assembly_list.end() && *it_assembly_innerlist != (*it_assembly_outerlist).front();
                    it_assembly_innerlist++)
                entry_location++;

            component_chains_list = merge_block_pairs_into_chains(block_pair_map, merged_sorted_block_list,
                                                                assembly_block_coords_map[ (*it_assembly_outerlist).front() ],
                                                                comp_idx, max_neighbour_separation, (*it_assembly_outerlist).front(), entry_location);
            cout<<comp_idx<<" : "<<component_chains_list.size() << " " << (*it_assembly_outerlist).size() <<"\n";
            // core_blocks_list = get_core_block_list(block_pair_map);
            // core_blocks_outerlist.push_back(core_blocks_list);
            // retained_comp_idx_vector.push_back(true);
            // pululli ref_comp_idx_pair = make_pair(reference_assembly_idx, comp_idx);
            // comp_reference_assembly_idx_list.push_back( make_pair(reference_assembly_idx, comp_idx) );
            // for(auto it = core_blocks_list.begin(); it != core_blocks_list.end(); it++)
            //     cout << *it << " ";
            // cout<<"\n\n";
            // /*location_sorted_blocks_list =*/ get_sorted_assembly_component_block_coords(core_blocks_list, merged_sorted_block_list,
            //                                                                         assembly_block_coords_map[reference_assembly_idx]);
            // check_coords_vector( assembly_block_coords_map[ (*it_assembly_outerlist).front() ] );

            all_component_chains_list.insert(all_component_chains_list.end(), component_chains_list.begin(), component_chains_list.end());
            cout<<all_component_chains_list.size()<<"\n";

            for(it_chain = component_chains_list.begin(); it_chain != component_chains_list.end(); it_chain++){
                chain_ends_block_list.push_back((*it_chain).ends_blocks.first);
                chain_ends_block_list.push_back((*it_chain).ends_blocks.second);
                // cout << "(" << (*it_chain).chain_idx.first << "," << (*it_chain).chain_idx.second << ") : ";
                // cout << "(" << (*it_chain).ends_blocks.first << "," << (*it_chain).ends_blocks.second << ") : ";
                // cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") : ";
                // cout << (((*it_chain).orientations_list.front()) ? "1" : "0") << "\n";
            }
            component_chains_list.clear();
        }

        comp_idx++;
        it_block_outerlist++;
        it_assembly_outerlist++;
    }
    cout << merged_sorted_block_list.size() << " -> " << all_component_chains_list.size() << " + |" << cummulative_discard_len << "|";
    cout << " " << chain_ends_block_list.size() << "\n";

    chain_ends_block_list.sort();

    count = 0;
    for(it_chain = all_component_chains_list.begin(); it_chain != all_component_chains_list.end() && count<25; it_chain++, count++){
        // if((*it_chain).ends_blocks.first == (*it_chain).ends_blocks.second)
        //     cout<<"#########################\n";
        cout << "<" << (*it_chain).reference_assembly_idx << "> : ";
        cout << "(" << (*it_chain).chain_idx.first << "," << (*it_chain).chain_idx.second << ") : ";
        cout << "(" << (*it_chain).ends_blocks.first << "," << (*it_chain).ends_blocks.second << ") : ";
        cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") : ";
        cout << (((*it_chain).orientations_list.front()) ? "1" : "0") << "\n";
        // if((*it_chain).ends_blocks.first == (*it_chain).ends_blocks.second)
        //     cout<<"#########################\n";
    }
    cout<<"\n";

    count = 0;
    for(it_block_innerlist = chain_ends_block_list.begin(); it_block_innerlist != chain_ends_block_list.end() && count<25;
            it_block_innerlist++, count++)
        cout << *it_block_innerlist << " ";
    cout<<"\n\n";


    //  Update chain coordinates across all assemblies containing the chains

    // ulli entry_location = 0;
    it_component_presence_outervector = assembly_component_presence_outervector.begin();

    for(it_assembly_innerlist = merged_sorted_assembly_list.begin(); it_assembly_innerlist != merged_sorted_assembly_list.end();
            it_assembly_innerlist++/*, entry_location++*/){
        update_chain_coordinates_lists(all_component_chains_list, chain_ends_block_list, merged_sorted_block_list,
                                        assembly_block_coords_map[ (*it_assembly_innerlist) ], *it_assembly_innerlist/*, entry_location*/,
                                        (*it_component_presence_outervector), start_comp_idx);
        it_component_presence_outervector++;
    }
    
    count = 0;
    for(it_chain = all_component_chains_list.begin(); it_chain != all_component_chains_list.end(); it_chain++){
        if(count<25 || (*it_chain).orientations_list.size() > merged_sorted_assembly_list.size() ){
            cout << "<" << (*it_chain).reference_assembly_idx << "> : ";
            cout << "(" << (*it_chain).chain_idx.first << "," << (*it_chain).chain_idx.second << ") : ";
            cout << "(" << (*it_chain).ends_blocks.first << "," << (*it_chain).ends_blocks.second << ") : ";
            cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") : ";
            cout << "(" << ((*it_chain).ends_coords_list.back()).first << "," << ((*it_chain).ends_coords_list.back()).second << ") : ";
            cout << (*it_chain).ends_coords_list.size() << " ";
            cout << (((*it_chain).orientations_list.front()) ? "1" : "0") << " " << (((*it_chain).orientations_list.back()) ? "1" : "0") << " ";
            cout << (*it_chain).orientations_list.size() << "\n";
            count++;
        }
    }
    cout<<"\n";

    
    comp_idx = start_comp_idx;
    it_assembly_outerlist = component_assembly_idx_outerlist.begin();
    for(it_chain = all_component_chains_list.begin(); it_chain != all_component_chains_list.end(); it_chain++){
        while((*it_chain).chain_idx.first != comp_idx){
            comp_idx++;
            it_assembly_outerlist++;
        }
        filter_chain_coordinates_lists( (*it_chain), *it_assembly_outerlist, merged_sorted_assembly_list );
    }

    cout<<"Filtered coordinates\n";
    
    count = 0;
    for(it_chain = all_component_chains_list.begin(); it_chain != all_component_chains_list.end() && count<25; it_chain++, count++){
        cout << "<" << (*it_chain).reference_assembly_idx << "> : ";
        cout << "(" << (*it_chain).chain_idx.first << "," << (*it_chain).chain_idx.second << ") : ";
        cout << "(" << (*it_chain).ends_blocks.first << "," << (*it_chain).ends_blocks.second << ") : ";
        cout << "(" << ((*it_chain).ends_coords_list.front()).first << "," << ((*it_chain).ends_coords_list.front()).second << ") : ";
        cout << "(" << ((*it_chain).ends_coords_list.back()).first << "," << ((*it_chain).ends_coords_list.back()).second << ") : ";
        cout << (*it_chain).ends_coords_list.size() << " ";
        cout << (((*it_chain).orientations_list.front()) ? "1" : "0") << " " << (((*it_chain).orientations_list.back()) ? "1" : "0") << " ";
        cout << (*it_chain).orientations_list.size() << "\n";
    }
    cout<<"\n";

    comp_idx = start_comp_idx;
    it_assembly_outerlist = component_assembly_idx_outerlist.begin();
    for(it_chain = all_component_chains_list.begin(); it_chain != all_component_chains_list.end(); it_chain++){
        // cout << (*it_chain).orientations_list.size() << " " << (*it_assembly_outerlist).size() <<"\n";
        while((*it_chain).chain_idx.first != comp_idx){
            comp_idx++;
            it_assembly_outerlist++;
        }
        if( (*it_chain).orientations_list.size() != (*it_assembly_outerlist).size() ){
            cout << "CHAIN_LIST_SIZE_MISMATCH_ISSUE!!! : <" << (*it_chain).reference_assembly_idx << "> : ";
            cout << "(" << (*it_chain).chain_idx.first << "," << (*it_chain).chain_idx.second << ") : ";
            cout << (*it_chain).orientations_list.size() << " " << (*it_assembly_outerlist).size() <<"\n";
        }
    }

    unsigned long iter_diff = all_component_chains_list.front().chain_idx.first - start_comp_idx;
    list< Metablock > metablock_list = generate_metablocks(all_component_chains_list, component_assembly_idx_outerlist, mapped_contig_ends,
                                                        max_neighbour_separation, use_oriented_links, mapped_oriented_links_map, iter_diff);

    cout << "\n" << merged_sorted_block_list.size() << " -> " << all_component_chains_list.size() << " + |" << cummulative_discard_len << "| --> ";
    cout<< all_component_chains_list.size() << " -> " << metablock_list.size( )<< "\n\n";

    cout << discard_count << "\t" << cummulative_discard_len << "\t" << ((float)cummulative_discard_len)/discard_count << "\t" << max_discard_len << "\n";

    save_metablocks(metablock_list, component_assembly_idx_outerlist, start_comp_idx, assembly_count, assemblies_per_partition, metablock_outfile_prefix);

    // for(it_assembly_innerlist = merged_sorted_assembly_list.begin(); it_assembly_innerlist != merged_sorted_assembly_list.end();
    //         it_assembly_innerlist++){
    //     cout << (*it_assembly_innerlist) << " : " << assembly_block_coords_map[ *it_assembly_innerlist ].size() << "\n";
    // }

    // cout<<"\n";
    // cout << core_blocks_outerlist.size() << " " << comp_reference_assembly_idx_list.size() << " " << retained_comp_idx_vector.size() << "\n";

    //for(auto it = comp_reference_assembly_idx_list.begin(); it != comp_reference_assembly_idx_list.end(); it++)
        //cout<< it->second << " , " << it->first /*<< " : " << assembly_block_coords_map[ it->first ].size()*/ << "\n";

    remove(comp_group_filename.c_str());

    return 0;
}