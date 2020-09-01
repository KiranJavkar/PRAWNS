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


map<ulli, string> get_metablock_tuple_idx_map(string metablock_dir, unsigned long group_count, vector<ulli>& count_offset){
    cout<<"get_metablock_tuple_idx_map started: "<<metablock_dir<<" "<<group_count<<"\n";
    map<ulli, string> metablock_tuple_idx_map;
    ulli current_idx = 0, metablock_count;
    tup_uib current_tuple;
    short split_no = 0 /*, count*/;
    string metablock_tuple_string;

    ifstream inFile;

    set<string> metablock_tuple_string_set;
    set<string>::iterator it_metablock_string_set;

    for(unsigned long group_no=0; group_no<group_count; group_no++){
        // cout<<"\t"<<group_no<<"\t";
        inFile.open( (metablock_dir + to_string(group_no) + "_metablock_indices").c_str(), ios::binary);

        inFile.read((char*) (&metablock_count), sizeof(metablock_count));
        // marks the starting mapping index for metablocks from a new group
        count_offset.push_back(current_idx);
        split_no = 0;

        while(metablock_count--){
            inFile.read((char*) (&current_tuple), sizeof(current_tuple));

            if( get<2>(current_tuple) )
                split_no++;
            else
                split_no = 0;

            metablock_tuple_string = to_string( get<0>(current_tuple) ) + "_" + to_string( get<1>(current_tuple) ) + "_" + to_string(split_no);

            // metablock_tuple_idx_map[current_idx++] = to_string( get<0>(current_tuple) ) + "_" + to_string( get<1>(current_tuple) ) + "_" + to_string(split_no);

            it_metablock_string_set = metablock_tuple_string_set.find(metablock_tuple_string);

            if(it_metablock_string_set != metablock_tuple_string_set.end()){
                cout<<"ALERT!!! METABLOCK_STRING_REENCOUNTERED_ISSUE!!! "<<current_idx<<" "<<metablock_tuple_string<<"\n";
            }
            else
                metablock_tuple_string_set.insert(metablock_tuple_string);

            metablock_tuple_idx_map[current_idx++] = metablock_tuple_string;

            /*if(count<10){
                cout << metablock_tuple_idx_map[current_idx-1] << "\n";
                count++;
            }*/
        }

        inFile.close();

        // cout<<metablock_tuple_idx_map.size()<<"\n";
        // count = 0;
    }

    cout<<"get_metablock_tuple_idx_map ended: "<<metablock_dir<<" "<<group_count<<" "<<count_offset.size();
    cout<<current_idx<<" "<<metablock_tuple_idx_map.size()<<"\n";

    return metablock_tuple_idx_map;
}


void init_partitioned_feature_lists(list< list< list<tup_uubub> > >& feature_partition_outerlist,
                                    list< list< list<tup_ububbbl> > >& paired_feature_partition_outerlist,
                                    short feature_partitions, unsigned long current_assembly_count){
    // tup_uubub: <start, end, strand, feature_idx, is_metablock>
    // tup_ububbbl: <feature_1_idx, is_1_metablock, feature_2_idx, is_2_metablock, rel_orientation_1, rel_orientation_2, inter_feature_separation>

    list< list<tup_uubub> > feature_partition_innerlist;
    list< list<tup_ububbbl> > paired_feature_partition_innerlist;

    list<tup_uubub> partitioned_assembly_features;
    list<tup_ububbbl> partitioned_assembly_paired_features;

    for(unsigned long assembly_idx_offset = 0; assembly_idx_offset < current_assembly_count; assembly_idx_offset++){
        feature_partition_innerlist.push_back(partitioned_assembly_features);
        paired_feature_partition_innerlist.push_back(partitioned_assembly_paired_features);
    }

    for(short partition_idx = 0; partition_idx < feature_partitions; partition_idx++){
        feature_partition_outerlist.push_back(feature_partition_innerlist);
        paired_feature_partition_outerlist.push_back(paired_feature_partition_innerlist);
    }
}


void add_feature_to_partitioned_list(list< list< list<tup_uubub> > >& feature_partition_outerlist, tup_uubub feature,
                                    unsigned long assembly_idx_offset, short partition_number){
    list< list< list<tup_uubub> > >::iterator it_feature_partition_outerlist = feature_partition_outerlist.begin();
    advance(it_feature_partition_outerlist, partition_number);

    list< list<tup_uubub> >:: iterator it_feature_partition_innerlist = (*it_feature_partition_outerlist).begin();
    advance(it_feature_partition_innerlist, assembly_idx_offset);

    (*it_feature_partition_innerlist).push_back(feature);
}


void add_feature_to_partitioned_list(list< list< list<tup_ububbbl> > >& paired_feature_partition_outerlist, tup_ububbbl paired_feature,
                                    unsigned long assembly_idx_offset, short smaller_partition_number){
    list< list< list<tup_ububbbl> > >::iterator it_paired_feature_partition_outerlist = paired_feature_partition_outerlist.begin();
    advance(it_paired_feature_partition_outerlist, smaller_partition_number);

    list< list<tup_ububbbl> >:: iterator it_paired_feature_partition_innerlist = (*it_paired_feature_partition_outerlist).begin();
    advance(it_paired_feature_partition_innerlist, assembly_idx_offset);

    (*it_paired_feature_partition_innerlist).push_back(paired_feature);
}


// uubub: start, end, orientation, block_idx/metablock_mapping_idx, is_metablock
// list< list< tup_uubub > > get_filtered_features(string metablock_dir, unsigned long group_count, string assembly_blocks_dir,
//                                                 unsigned long start_assembly_idx, unsigned long end_assembly_idx, vector<ulli> count_offset,
//                                                 string partition_idx_str, unsigned long min_sv_size, string contig_len_dir, vector<ulli>& max_metablock_idx){
void generate_filtered_features(string metablock_dir, unsigned long group_count, string assembly_blocks_dir, unsigned long start_assembly_idx,
                                unsigned long end_assembly_idx, vector<ulli> count_offset, string partition_idx_str, unsigned long min_sv_size,
                                list< list< list<tup_uubub> > >& feature_partition_outerlist,
                                list< list< list<tup_ububbbl> > >& paired_feature_partition_outerlist,
                                ulli blocks_per_partition, ulli metablocks_per_partition, short feature_partitions, ulli max_inter_feature_separation,
                                bool use_oriented_links, string oriented_links_dir, string contig_len_dir, vector<ulli>& max_metablock_idx){
    cout<<"generate_filtered_features started: "<<metablock_dir<<" "<<group_count<<" "<<partition_idx_str<<" "<<start_assembly_idx<<" "<<end_assembly_idx<<" ";
    cout<<count_offset.size()<<" "<<partition_idx_str<<" "<<min_sv_size<<" "<<feature_partition_outerlist.size()<<" ";
    cout<<paired_feature_partition_outerlist.size()<<" "<<blocks_per_partition<<" "<<metablocks_per_partition<<" "<<feature_partitions<<" ";
    cout<<max_inter_feature_separation<<" "<<use_oriented_links<<" "<<oriented_links_dir<<" "<<contig_len_dir<<" "<<max_metablock_idx.size()<<"\n";

    tup_uubu current_coords;
    ulli metablock_count, current_metablock_size, prev_contig_end;
    unsigned long assembly_idx, assembly_idx_offset=0;
    list< list< tup_uubub > > feature_coords_outerlist, block_coords_outerlist;
    list< list< tup_uubub > >::iterator it_coords_outerlist, it_block_coords_outerlist;
    list< tup_uubub > block_coords_innerlist;
    list< tup_uubub >::iterator it_coords_innerlist, it_block_coords_innerlist;
    vector<ulli>metablock_coverage, cumulative_length, max_coord_vec, max_metablock_size, min_metablock_size, metablock_counts, contig_ends_vector;
    vector< tup_uub >largest_metablock_vector;
    unsigned int contig_idx;
    map< unsigned long, vector<ulli> > mapped_contig_ends;
    map< unsigned int, ptupuubbtupuubb> contig_ends_feature_map;
    map< unsigned int, ptupuubbtupuubb>::iterator it_contig_ends_feature_map;
    tup_uubb contig_start_feature_tuple, contig_end_feature_tuple; // < distance_from_start, feature_idx, is_metablock, strand >
    map< puiui, tup_bbffs > oriented_links_map;
    map< puiui, tup_bbffs >::iterator it_ol_map;

    string filename;
    list<string> lines, lsplit;
    list<string>::iterator it_coords;
    ulli current_block_idx, diff, start_pos, end_pos, discard_count, previous_contig_end, distance, feature_count, paired_feature_count;
    bool orientation, is_new_contig;
    long separation;
    list< tup_uubub > previous_features;
    list< tup_uubub >::iterator it_prev_feat;

    short current_feature_partition_no /*, previous_feature_partition_no*/, count = 0;
    list<short> previous_feature_partition_numbers;
    list<short>::iterator it_prev_partition_no;

    // list< list< list<tup_uubub> > > feature_partition_outerlist;
    // list< list< list<tup_ububbbl> > > paired_feature_partition_outerlist;
    tup_uubub previous_feature_tuple;
    tup_ububbbl paired_feature;

    // init_partitioned_feature_lists( feature_partition_outerlist, paired_feature_partition_outerlist,
    //                                 feature_partitions, end_assembly_idx - start_assembly_idx + 1);

    // Load metablocks - unsorted by design

    for(assembly_idx=start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++, assembly_idx_offset++){
        list< tup_uubub > current_coords_list;
        feature_coords_outerlist.push_back(current_coords_list);
        metablock_coverage.push_back(0);
        cumulative_length.push_back(0);
        max_coord_vec.push_back(0);
        max_metablock_size.push_back(0);
        min_metablock_size.push_back(100);
        metablock_counts.push_back(0);
        largest_metablock_vector.push_back( tup_uub(0,0,false) );

        contig_ends_vector = get_contig_ends_vector(contig_len_dir + to_string(assembly_idx));
        mapped_contig_ends[assembly_idx_offset] = contig_ends_vector;
    }

    // cout<<"CONTIG ENDS VECTORS LOADED\n";

    ifstream inFile;

    for(unsigned long group_no=0; group_no<group_count; group_no++){
        inFile.open( (metablock_dir + to_string(group_no) + "_" + partition_idx_str).c_str(), ios::binary);
        // cout<<metablock_dir + to_string(group_no) + "_" + partition_idx_str<<"\n";
        assembly_idx_offset = 0;

        for(it_coords_outerlist = feature_coords_outerlist.begin(); it_coords_outerlist != feature_coords_outerlist.end();
                it_coords_outerlist++, assembly_idx_offset++){

            inFile.read((char*) (&metablock_count), sizeof(metablock_count));
            // cout<<group_no<<" "<<metablock_count<<"\n";

            while(metablock_count--){
                inFile.read((char*) (&current_coords), sizeof(current_coords));

                // if(count<5){
                //     cout << get<0>(current_coords) << " " << get<1>(current_coords) << " " << get<2>(current_coords);
                //     cout << " " << get<3>(current_coords) << " " << count_offset[group_no] << " ";
                //     cout<< get<3>(current_coords) + count_offset[group_no] << " 1\n";
                //     count++;
                // }
                
                (*it_coords_outerlist).push_back( tup_uubub( get<0>(current_coords), get<1>(current_coords), get<2>(current_coords),
                                                            get<3>(current_coords) + count_offset[group_no], true ) );

                // if(get<1>(current_coords) - get<0>(current_coords) + 1 >= min_sv_size)
                //     (*it_coords_outerlist).push_back( tup_uubub( get<0>(current_coords), get<1>(current_coords), get<2>(current_coords),
                //                                             get<3>(current_coords) + count_offset[group_no], true ) );

                if(get<0>(current_coords) > get<1>(current_coords)){
                    cout << "ALERT!!! " << assembly_idx_offset<< " " << get<0>(current_coords) << " " << get<1>(current_coords) << " ";
                    cout<< get<2>(current_coords) << " " << get<3>(current_coords) + count_offset[group_no] << " 1\n";
                }
                current_metablock_size = get<1>(current_coords) - get<0>(current_coords) + 1;
                metablock_coverage[assembly_idx_offset] += current_metablock_size;
                max_coord_vec[assembly_idx_offset] = (max_coord_vec[assembly_idx_offset] < (get<1>(current_coords)))?(get<1>(current_coords)):max_coord_vec[assembly_idx_offset];
                if(max_metablock_size[assembly_idx_offset] < current_metablock_size){
                    max_metablock_size[assembly_idx_offset] = current_metablock_size;
                    largest_metablock_vector[assembly_idx_offset] = tup_uub(get<0>(current_coords), get<1>(current_coords), get<2>(current_coords));
                    max_metablock_idx[assembly_idx_offset] = get<3>(current_coords) + count_offset[group_no];
                }
                if(min_metablock_size[assembly_idx_offset] > current_metablock_size)
                    min_metablock_size[assembly_idx_offset] = current_metablock_size;
            }
        }
        inFile.close();
        count = 0;
    }

    // cout<<"METABLOCKS LOADED\n";

    
    // Sort metablock lists by location

    assembly_idx_offset = 0;
    for(it_coords_outerlist = feature_coords_outerlist.begin(); it_coords_outerlist != feature_coords_outerlist.end();
            it_coords_outerlist++, assembly_idx_offset++){

        (*it_coords_outerlist).sort();
        metablock_counts[assembly_idx_offset] = (*it_coords_outerlist).size();

        // cout<< start_assembly_idx+assembly_idx_offset << " : " << (*it_coords_outerlist).size() << " " << max_metablock_size[assembly_idx_offset] << " ";
        // cout<< metablock_coverage[assembly_idx_offset]*1.0/((*it_coords_outerlist).size()) << " "<< metablock_coverage[assembly_idx_offset] << " ";
        // cout<< max_coord_vec[assembly_idx_offset] << " " << metablock_coverage[assembly_idx_offset]*1.0/max_coord_vec[assembly_idx_offset] << "\n";
    }

    // cout<<"METABLOCKS SORTED\n";


    // Load and sort individual blocks above required minimum size

    assembly_idx_offset = 0;
    for(assembly_idx=start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++, assembly_idx_offset++){
        block_coords_innerlist.clear();
        filename = assembly_blocks_dir + to_string(assembly_idx);
        lines = get_lines(filename.c_str());
        lsplit = split(lines.front(), ',');
        it_coords = lsplit.begin();
        current_block_idx = 0;
        discard_count = 0;

        while(it_coords != lsplit.end()){
            if( *it_coords == "0"){
                advance(it_coords, 3);
                current_block_idx++;
                continue;
            }
            start_pos = strtoulli( *it_coords++);
            end_pos = strtoulli( *it_coords++);
            orientation = ( *it_coords++ == "1" );
            max_coord_vec[assembly_idx_offset] = (max_coord_vec[assembly_idx_offset] < end_pos)?end_pos:max_coord_vec[assembly_idx_offset];

            if(end_pos - start_pos + 1 >= min_sv_size)
                block_coords_innerlist.push_back(tup_uubub( start_pos, end_pos, orientation, current_block_idx, false ));
            else
                discard_count++;

            current_block_idx++;
        }

        block_coords_innerlist.sort();

        block_coords_outerlist.push_back(block_coords_innerlist);
        cout<<"\t"<<assembly_idx<<" : "<<block_coords_innerlist.size()<<"\t"<<discard_count<<"\n";
    }

    // cout<<"LONGER BLOCKS LOADED\n";


    // Add screened blocks to combined feature lists

    it_block_coords_outerlist = block_coords_outerlist.begin();
    assembly_idx=start_assembly_idx;
    for(it_coords_outerlist = feature_coords_outerlist.begin(); it_coords_outerlist != feature_coords_outerlist.end();
            it_coords_outerlist++,it_block_coords_outerlist++, assembly_idx++){

        cout<<assembly_idx<<" : "<<(*it_coords_outerlist).size()<<" ("<< (*it_coords_outerlist).size() + (*it_block_coords_outerlist).size()<<")"<<" -> ";

        it_coords_innerlist = (*it_coords_outerlist).begin();
        // ulli block_counter=0;

        for(it_block_coords_innerlist = (*it_block_coords_outerlist).begin(); it_block_coords_innerlist != (*it_block_coords_outerlist).end();
                it_block_coords_innerlist++/*, block_counter++*/){
            // if(block_counter%1000==0){
            //     cout<<"\t\t"<<get<0>( *it_coords_innerlist)<<" "<<get<1>( *it_coords_innerlist)<<"        ";
            //     cout<<get<0>( *it_block_coords_innerlist)<<" "<<get<1>( *it_block_coords_innerlist)<<"\n";
            // }

            // while metablock end is less than current block end, increment the metablock iterator
            // note that if the metablock and block ends match, block is contained inside the metablock - this case is handled in the later if clause
            while( it_coords_innerlist != (*it_coords_outerlist).end() && get<1>( *it_coords_innerlist) < get<1>( *it_block_coords_innerlist) )
                it_coords_innerlist++;

            if(it_coords_innerlist == (*it_coords_outerlist).end()){
                // the remaining blocks are now necessarily ending after the last metablock
                (*it_coords_outerlist).insert( it_coords_innerlist, it_block_coords_innerlist, (*it_block_coords_outerlist).end() );
                it_block_coords_innerlist = (*it_block_coords_outerlist).end();
                it_block_coords_innerlist--;
            }
            else{

                // // if metablock start is greater than current block start (block preceeds completely or overlaps at start of metablock)
                // // insert the block before the metablock
                // // note that even after the insert, the iterator would point to the metablock location and not the newly inserted block
                // if( get<0>( *it_coords_innerlist) > get<0>( *it_block_coords_innerlist) )
                //     (*it_coords_outerlist).insert( it_coords_innerlist, *it_block_coords_innerlist );
                // // else if the metablock end is less than the current block end (current block is not completely contained inside the metablock)
                // // insert the block after the metablock
                // // Note - here the iterator would point to the metablock (if available) starting after the start of the newly inserted block
                // else if( get<1>( *it_coords_innerlist) < get<1>( *it_block_coords_innerlist) ){
                //     it_coords_innerlist++;
                //     if(it_coords_innerlist == (*it_coords_outerlist).end()){
                //         (*it_coords_outerlist).insert( it_coords_innerlist, it_block_coords_innerlist, (*it_block_coords_outerlist).end() );
                //         it_block_coords_innerlist = (*it_block_coords_outerlist).end();
                //         it_block_coords_innerlist--;
                //     }
                //     else{
                //         (*it_coords_outerlist).insert( it_coords_innerlist, *it_block_coords_innerlist );
                //     }
                // }

                // Current block should not be contained inside the metablock (metablock starts after block start)
                if( get<0>( *it_coords_innerlist) > get<0>( *it_block_coords_innerlist) )
                    (*it_coords_outerlist).insert( it_coords_innerlist, *it_block_coords_innerlist );
            }
        }

        cout<<(*it_coords_outerlist).size()<<"\n";
    }

    // cout<<"METABLOCKS AND LONGER BLOCKS MERGED AND SORTED\n";


    // Generate partitioned individual and paired features

    assembly_idx_offset = 0;
    assembly_idx=start_assembly_idx;
    for(it_coords_outerlist = feature_coords_outerlist.begin(); it_coords_outerlist != feature_coords_outerlist.end();
            it_coords_outerlist++, assembly_idx++, assembly_idx_offset++){
        ulli last_pos = 0;
        cumulative_length[assembly_idx_offset] = 0;

        contig_ends_feature_map.clear();
        is_new_contig = true;
        contig_ends_vector = mapped_contig_ends[assembly_idx_offset];
        contig_idx = 0;
        previous_contig_end = 0;
        feature_count = 0;
        paired_feature_count = 0;
        previous_features.clear();
        previous_feature_partition_numbers.clear();
        // it_prev_feat = previous_features.begin();

        for(it_coords_innerlist = (*it_coords_outerlist).begin(); it_coords_innerlist != (*it_coords_outerlist).end(); it_coords_innerlist++){
            if(get<0>(*it_coords_innerlist) > last_pos)
                last_pos = get<0>(*it_coords_innerlist);
            cumulative_length[assembly_idx_offset] += get<1>(*it_coords_innerlist)-last_pos+1;
            last_pos = get<1>(*it_coords_innerlist);

            while( get<3>( *it_coords_innerlist) > contig_ends_vector[contig_idx] ){
                // New contig
                // if current contig had at least 1 feature generated ( ==> is_new_contig == false)
                // set this feature as the contig end feature for the contig (if use_oriented_links == true)

                if(use_oriented_links && !is_new_contig){
                    // (distance_from_end, feature_2_idx, is_metablock, strand)
                    // distance = contig_ends_vector[contig_idx] - get<1>(previous_feature_tuple);
                    // contig_end_feature_tuple = make_tuple(  distance, get<3>(previous_feature_tuple),
                    //                                         get<4>(previous_feature_tuple), get<2>(previous_feature_tuple) );

                    distance = contig_ends_vector[contig_idx] - get<1>(previous_features.back());
                    contig_end_feature_tuple = make_tuple(  distance, get<3>(previous_features.back()),
                                                            get<4>(previous_features.back()), get<2>(previous_features.back()));

                    contig_ends_feature_map[contig_idx] = make_pair( contig_start_feature_tuple, contig_end_feature_tuple);
                }

                previous_contig_end = contig_ends_vector[contig_idx];
                contig_idx++;
                is_new_contig = true;
                previous_features.clear();
                previous_feature_partition_numbers.clear();
            }

            if( get<4>( *it_coords_innerlist ) ) // metablock
                current_feature_partition_no = (get<3>( *it_coords_innerlist ))/metablocks_per_partition;
            else // block
                current_feature_partition_no = (get<3>( *it_coords_innerlist ))/blocks_per_partition;

            // This would happen in case the metablock index is greater than the estimated total metablock count
            if(current_feature_partition_no >= feature_partitions)
                current_feature_partition_no = feature_partitions - 1;

            add_feature_to_partitioned_list(feature_partition_outerlist, *it_coords_innerlist, assembly_idx_offset, current_feature_partition_no);
            if(feature_count%10000 == 0){
                cout << "\t\t (" << get<0>(*it_coords_innerlist) << "," << get<1>(*it_coords_innerlist) << "," << get<2>(*it_coords_innerlist) << ",";
                cout << get<3>(*it_coords_innerlist) << "," << get<4>(*it_coords_innerlist) << ") : " << current_feature_partition_no << " ";
            }
            feature_count++;

            if(is_new_contig){
                if(use_oriented_links){
                    // (distance_from_start, feature_1_idx, is_metablock, strand)
                    distance = previous_contig_end - get<0>( *it_coords_innerlist ) - 1;
                    contig_start_feature_tuple = make_tuple(distance, get<3>( *it_coords_innerlist ),
                                                            get<4>( *it_coords_innerlist ), get<2>( *it_coords_innerlist ));
                }
            }
            else{
                /*if( get<3>( previous_feature_tuple ) <= get<3>( *it_coords_innerlist ) ){
                    // ( feature_1_idx, is_1_metablock, feature_2_idx, is_2_metablock, rel_orientation_1, rel_orientation_2, inter_feature_separation)
                    paired_feature = make_tuple(get<3>( previous_feature_tuple ), get<4>( previous_feature_tuple ),
                                                get<3>( *it_coords_innerlist ), get<4>( *it_coords_innerlist ),
                                                get<2>( previous_feature_tuple), get<2>( *it_coords_innerlist ),
                                                get<0>( *it_coords_innerlist ) - get<1>( previous_feature_tuple ) - 1 );
                    add_feature_to_partitioned_list(paired_feature_partition_outerlist, paired_feature, assembly_idx_offset,
                                                    previous_feature_partition_no);
                }
                else{
                    paired_feature = make_tuple(get<3>( *it_coords_innerlist ), get<4>( *it_coords_innerlist ),
                                                get<3>( previous_feature_tuple ), get<4>( previous_feature_tuple ),
                                                !(get<2>( *it_coords_innerlist )), !(get<2>( previous_feature_tuple)),
                                                get<0>( *it_coords_innerlist ) - get<1>( previous_feature_tuple ) - 1 );
                    add_feature_to_partitioned_list(paired_feature_partition_outerlist, paired_feature, assembly_idx_offset,
                                                    current_feature_partition_no);
                }  

                if(paired_feature_count%10000 == 0){
                    cout << "\t\t p(" << get<0>(paired_feature) << "," << get<1>(paired_feature) << "," << get<2>(paired_feature) << ",";
                    cout << get<3>(paired_feature) << "," << get<4>(paired_feature) << "," << get<5>(paired_feature) << ",";
                    cout << get<6>(paired_feature) << ") : ";
                    cout << (get<3>(previous_feature_tuple)<=get<3>(*it_coords_innerlist)?previous_feature_partition_no:current_feature_partition_no) <<" ";
                }
                paired_feature_count++;*/

                it_prev_feat = previous_features.begin();
                it_prev_partition_no = previous_feature_partition_numbers.begin();

                while(it_prev_feat != previous_features.end()){
                    separation = static_cast<long> (get<0>( *it_coords_innerlist ) - get<1>( *it_prev_feat ) - 1);
                    if(separation > max_inter_feature_separation){
                        it_prev_feat = previous_features.erase(it_prev_feat);
                        it_prev_partition_no = previous_feature_partition_numbers.erase(it_prev_partition_no);
                    }
                    else{
                        // if( ( get<3>( *it_prev_feat ) != get<3>( *it_coords_innerlist ) ) ||
                        //         ( ( get<3>( *it_prev_feat ) == get<3>( *it_coords_innerlist ) ) && ( get<4>( *it_prev_feat ) != get<4>( *it_coords_innerlist ) ) ) )

                        if( get<3>( *it_prev_feat ) <= get<3>( *it_coords_innerlist ) ){
                            // ( feature_1_idx, is_1_metablock, feature_2_idx, is_2_metablock, rel_orientation_1, rel_orientation_2, inter_feature_separation)
                            paired_feature = make_tuple(get<3>( *it_prev_feat ), get<4>( *it_prev_feat ),
                                                        get<3>( *it_coords_innerlist ), get<4>( *it_coords_innerlist ),
                                                        get<2>( *it_prev_feat), get<2>( *it_coords_innerlist ), separation);
                            add_feature_to_partitioned_list(paired_feature_partition_outerlist, paired_feature, assembly_idx_offset,
                                                            *it_prev_partition_no);
                        }
                        else{
                            paired_feature = make_tuple(get<3>( *it_coords_innerlist ), get<4>( *it_coords_innerlist ),
                                                        get<3>( *it_prev_feat ), get<4>( *it_prev_feat ),
                                                        !(get<2>( *it_coords_innerlist )), !(get<2>( *it_prev_feat)), separation);
                            add_feature_to_partitioned_list(paired_feature_partition_outerlist, paired_feature, assembly_idx_offset,
                                                            current_feature_partition_no);
                        }  

                        if(paired_feature_count%10000 == 0){
                            cout << "\t\t p(" << get<0>(paired_feature) << "," << get<1>(paired_feature) << "," << get<2>(paired_feature) << ",";
                            cout << get<3>(paired_feature) << "," << get<4>(paired_feature) << "," << get<5>(paired_feature) << ",";
                            cout << get<6>(paired_feature) << ") : ";
                            cout << (get<3>(*it_coords_innerlist)<=get<3>(*it_coords_innerlist)?(*it_prev_partition_no):current_feature_partition_no) <<" ";
                        }

                        if((get<0>(paired_feature)==0 && get<1>(paired_feature)==false) || (get<2>(paired_feature)==0 && get<3>(paired_feature)==false)){
                            cout << "\t\tALERT!!! CHECK PAIR!!! p(" << get<0>(paired_feature) << "," << get<1>(paired_feature) << "," << get<2>(paired_feature) << ",";
                            cout << get<3>(paired_feature) << "," << get<4>(paired_feature) << "," << get<5>(paired_feature) << ",";
                            cout << get<6>(paired_feature) << ") : ";
                            cout << (get<3>(*it_coords_innerlist)<=get<3>(*it_coords_innerlist)?(*it_prev_partition_no):current_feature_partition_no) <<" ";
                        }

                        paired_feature_count++;
                        it_prev_feat++;
                        it_prev_partition_no++;
                    }
                }
            }

            is_new_contig = false;
            // previous_feature_partition_no = current_feature_partition_no;
            // previous_feature_tuple = (*it_coords_innerlist);
            previous_features.push_back( *it_coords_innerlist );
            previous_feature_partition_numbers.push_back(current_feature_partition_no);
        }

        cout << "\n" << feature_count << " " << paired_feature_count << "\n";

        if(use_oriented_links){

            if(!is_new_contig){
                // last feature detected
                contig_end_feature_tuple = make_tuple( contig_ends_vector[contig_idx] - get<1>(previous_feature_tuple),
                                                get<3>(previous_feature_tuple), get<4>(previous_feature_tuple), get<2>(previous_feature_tuple));

                contig_ends_feature_map[contig_idx] = make_pair( contig_start_feature_tuple, contig_end_feature_tuple);
            }


            // Add paired features formed using contig ends features

            oriented_links_map = load_oriented_links_map(oriented_links_dir + to_string(assembly_idx));

            bool rel_orientation_1, rel_orientation_2;
            short partition_idx_1, partition_idx_2;
            ulli feature_1_idx, feature_2_idx;

            for(it_ol_map = oriented_links_map.begin(); it_ol_map != oriented_links_map.end(); it_ol_map++){
                //
                rel_orientation_1 = get<0>(it_ol_map->second);
                rel_orientation_2 = get<1>(it_ol_map->second);

                // contig_start_feature_tuple ==> left contig tuple
                // contig_end_feature_tuple ==> right contig tuple

                it_contig_ends_feature_map = contig_ends_feature_map.find( (it_ol_map->first).first );

                if(it_contig_ends_feature_map == contig_ends_feature_map.end() )
                    continue;

                // if(rel_orientation_1)
                //     contig_start_feature_tuple = contig_ends_feature_map[ (it_ol_map->first).first ].first;
                // else
                //     contig_start_feature_tuple = contig_ends_feature_map[ (it_ol_map->first).first ].second;

                if(rel_orientation_1)
                    contig_start_feature_tuple = (it_contig_ends_feature_map->second).first;
                else
                    contig_start_feature_tuple = (it_contig_ends_feature_map->second).second;


                it_contig_ends_feature_map = contig_ends_feature_map.find( (it_ol_map->first).second );

                if(it_contig_ends_feature_map == contig_ends_feature_map.end() )
                    continue;

                // if(rel_orientation_2)
                //     contig_end_feature_tuple = contig_ends_feature_map[ (it_ol_map->first).second ].first;
                // else
                //     contig_end_feature_tuple = contig_ends_feature_map[ (it_ol_map->first).second ].second;

                if(rel_orientation_2)
                    contig_end_feature_tuple = (it_contig_ends_feature_map->second).first;
                else
                    contig_end_feature_tuple = (it_contig_ends_feature_map->second).second;


                separation = get<0>(contig_start_feature_tuple) + get<0>(contig_end_feature_tuple);// + round(get<2>(it_ol_map->second)); //+ get<3>(it_ol_map->second);

                // Add the paired feature only if the interfeature separation is acceptable
                if(separation > max_inter_feature_separation && (separation + round(get<2>(it_ol_map->second))) > max_inter_feature_separation)
                    continue;

                feature_1_idx = get<1>( contig_start_feature_tuple );
                feature_2_idx = get<1>( contig_end_feature_tuple );


                if( get<3>( contig_start_feature_tuple ) ) // metablock
                    partition_idx_1 = feature_1_idx/metablocks_per_partition;
                else // block
                    partition_idx_1 = feature_1_idx/blocks_per_partition;

                // This would happen in case the metablock index is greater than the estimated total metablock count
                if(partition_idx_1 >= feature_partitions)
                    partition_idx_1 = feature_partitions - 1;


                if( get<3>( contig_end_feature_tuple ) ) // metablock
                    partition_idx_2 = feature_2_idx/metablocks_per_partition;
                else // block
                    partition_idx_2 = feature_2_idx/blocks_per_partition;

                // This would happen in case the metablock index is greater than the estimated total metablock count
                if(partition_idx_2 >= feature_partitions)
                    partition_idx_2 = feature_partitions - 1;


                if(rel_orientation_1){
                    // contig end: 'B' ==> toggle the actual feature strand
                    rel_orientation_1 = !( get<3>(contig_start_feature_tuple) );
                }
                else
                    rel_orientation_1 = get<3>(contig_start_feature_tuple);

                if(rel_orientation_2)
                    rel_orientation_2 = get<3>(contig_end_feature_tuple);
                else
                    rel_orientation_2 = !( get<3>(contig_end_feature_tuple) );

                if(feature_1_idx <= feature_2_idx){
                    paired_feature = make_tuple(feature_1_idx, get<2>( contig_start_feature_tuple ),
                                                feature_2_idx, get<1>( contig_end_feature_tuple ),
                                                rel_orientation_1, rel_orientation_2, separation);
                    add_feature_to_partitioned_list(paired_feature_partition_outerlist, paired_feature, assembly_idx_offset,
                                                    partition_idx_1);
                }
                else{
                    paired_feature = make_tuple(feature_2_idx, get<1>( contig_end_feature_tuple ),
                                                feature_1_idx, get<2>( contig_start_feature_tuple ),
                                                !rel_orientation_2, !rel_orientation_1, separation);
                    add_feature_to_partitioned_list(paired_feature_partition_outerlist, paired_feature, assembly_idx_offset,
                                                    partition_idx_2);
                }

                if(paired_feature_count%10000 == 0){
                    cout << "\t\t p(" << get<0>(paired_feature) << "," << get<1>(paired_feature) << "," << get<2>(paired_feature) << ",";
                    cout << get<3>(paired_feature) << "," << get<4>(paired_feature) << "," << get<5>(paired_feature) << ",";
                    cout << get<6>(paired_feature) << ") : ";
                    cout << (feature_1_idx<=feature_2_idx?partition_idx_1:partition_idx_2) <<" ";
                }


                if((get<0>(paired_feature)==0 && get<1>(paired_feature)==false) || (get<2>(paired_feature)==0 && get<3>(paired_feature)==false)){
                    cout << "\t\tALERT!!! CHECK PAIR (INTERCONTIG) !!! p(" << get<0>(paired_feature) << "," << get<1>(paired_feature) << "," << get<2>(paired_feature) << ",";
                    cout << get<3>(paired_feature) << "," << get<4>(paired_feature) << "," << get<5>(paired_feature) << ",";
                    cout << get<6>(paired_feature) << ") : ";
                    cout << (get<3>(*it_coords_innerlist)<=get<3>(*it_coords_innerlist)?(*it_prev_partition_no):current_feature_partition_no) <<" ";
                }

                paired_feature_count++;
            }
        }

        contig_idx = get_contig_index_by_coordinate(contig_ends_vector, get<0>(largest_metablock_vector[assembly_idx_offset]), 0,
                                                    contig_ends_vector.size() - 1 );
        prev_contig_end = (contig_idx>0)?(contig_ends_vector[contig_idx-1]):0;

        cout<< start_assembly_idx+assembly_idx_offset << " : " << (*it_coords_outerlist).size() << " ";
        cout<< min_metablock_size[assembly_idx_offset] << " " << max_metablock_size[assembly_idx_offset] << " ";
        cout<< "(" << get<0>(largest_metablock_vector[assembly_idx_offset]) << "," << get<1>(largest_metablock_vector[assembly_idx_offset]) << ",";
        cout<< ( ( get<2>(largest_metablock_vector[assembly_idx_offset]) )?"1":"0" ) << ") :: ";
        cout<< "[ " << contig_idx+1 << ", (" <<  get<0>(largest_metablock_vector[assembly_idx_offset]) - prev_contig_end << ",";
        cout<< get<1>(largest_metablock_vector[assembly_idx_offset]) - prev_contig_end << ") ] ";
        cout<< metablock_coverage[assembly_idx_offset]*1.0/metablock_counts[assembly_idx_offset] << " "<< cumulative_length[assembly_idx_offset] << " ";
        cout<< max_coord_vec[assembly_idx_offset] << " " << cumulative_length[assembly_idx_offset]*1.0/max_coord_vec[assembly_idx_offset] << "\t";
        cout<< feature_count << " " << paired_feature_count << "\n";
    }

    // return feature_coords_outerlist;
    cout<<"generate_filtered_features ended: "<<metablock_dir<<" "<<group_count<<" "<<partition_idx_str<<"\n";
}


void save_partitioned_features( list< list< list<tup_uubub> > > feature_partition_outerlist,
                                list< list< list<tup_ububbbl> > > paired_feature_partition_outerlist,
                                map<ulli, string> metablock_tuple_idx_map, string start_assembly_idx_str, string end_assembly_idx_str,
                                string filtered_feature_outdir, string paired_feature_outdir){
    cout<<"save_partitioned_features started: "<<start_assembly_idx_str<<" "<<end_assembly_idx_str<<"\n";

    string filename, feature_tuple_string;
    // metablock_tuple_string: <comp_idx>_<current_metablock_no>_<current_metablock_split_no>

    short partition_idx = 0, feature_tuple_string_length;
    ulli partitioned_assembly_feature_count, feature_idx;
    tup_uubb current_coords;
    tup_bbbbl current_paired_coords;

    list< list< list<tup_uubub> > >::iterator it_feature_partition_outerlist;
    list< list< list<tup_ububbbl> > >::iterator it_paired_feature_partition_outerlist = paired_feature_partition_outerlist.begin();

    list< list<tup_uubub> >:: iterator it_feature_partition_innerlist; 
    list< list<tup_ububbbl> >:: iterator it_paired_feature_partition_innerlist; // = (*it_paired_feature_partition_outerlist).begin();

    ofstream outFile;

    // Save partitioned metablocks and screen blocks

    for(it_feature_partition_outerlist = feature_partition_outerlist.begin(); it_feature_partition_outerlist != feature_partition_outerlist.end();
            it_feature_partition_outerlist++, partition_idx++){
        
        filename = filtered_feature_outdir + start_assembly_idx_str + "_" + end_assembly_idx_str + "_" + to_string(partition_idx);
        cout<<"\t"<<filename<<"\n";

        outFile.open(filename.c_str(), ios::binary);

        for(it_feature_partition_innerlist = (*it_feature_partition_outerlist).begin();
                it_feature_partition_innerlist != (*it_feature_partition_outerlist).end(); it_feature_partition_innerlist++){

            partitioned_assembly_feature_count = (*it_feature_partition_innerlist).size();
            outFile.write((char*) (&partitioned_assembly_feature_count), sizeof(partitioned_assembly_feature_count));

            for(auto it_coords = (*it_feature_partition_innerlist).begin(); it_coords != (*it_feature_partition_innerlist).end(); it_coords++){

                feature_idx = get<3>( *it_coords);
                
                if( get<4>( *it_coords) ){
                    // metablock
                    feature_tuple_string = metablock_tuple_idx_map[feature_idx];
                }
                else
                    feature_tuple_string = to_string(feature_idx);

                feature_tuple_string_length = feature_tuple_string.size();

                outFile.write((char*)(&feature_tuple_string_length), sizeof(feature_tuple_string_length));
                outFile.write( &feature_tuple_string[0], feature_tuple_string_length);

                // <start, end, orientation, is_metablock>
                current_coords = make_tuple( get<0>(*it_coords), get<1>(*it_coords), get<2>(*it_coords), get<4>(*it_coords) );
                outFile.write((char*) (&current_coords), sizeof(current_coords));
            }
        }

        outFile.close();
    }

    cout<<"Individual features saved\n";

    // Save partitioned paired features

    partition_idx = 0;
    for(it_paired_feature_partition_outerlist = paired_feature_partition_outerlist.begin();
            it_paired_feature_partition_outerlist != paired_feature_partition_outerlist.end(); it_paired_feature_partition_outerlist++, partition_idx++){
        
        filename = paired_feature_outdir + start_assembly_idx_str + "_" + end_assembly_idx_str + "_" + to_string(partition_idx);
        cout<<"\t"<<filename<<"\n";

        outFile.open(filename.c_str(), ios::binary);

        for(it_paired_feature_partition_innerlist = (*it_paired_feature_partition_outerlist).begin();
                it_paired_feature_partition_innerlist != (*it_paired_feature_partition_outerlist).end(); it_paired_feature_partition_innerlist++){

            partitioned_assembly_feature_count = (*it_paired_feature_partition_innerlist).size();
            outFile.write((char*) (&partitioned_assembly_feature_count), sizeof(partitioned_assembly_feature_count));

            for(auto it_coords = (*it_paired_feature_partition_innerlist).begin(); it_coords != (*it_paired_feature_partition_innerlist).end(); it_coords++){

                // feature 1 of the pair
                feature_idx = get<0>( *it_coords);
                
                if( get<1>( *it_coords) ){
                    // metablock
                    feature_tuple_string = metablock_tuple_idx_map[feature_idx];
                }
                else
                    feature_tuple_string = to_string(feature_idx);

                feature_tuple_string_length = feature_tuple_string.size();

                outFile.write((char*)(&feature_tuple_string_length), sizeof(feature_tuple_string_length));
                outFile.write( &feature_tuple_string[0], feature_tuple_string_length);


                // feature 2 of the pair
                feature_idx = get<2>( *it_coords);
                
                if( get<3>( *it_coords) ){
                    // metablock
                    feature_tuple_string = metablock_tuple_idx_map[feature_idx];
                }
                else
                    feature_tuple_string = to_string(feature_idx);

                feature_tuple_string_length = feature_tuple_string.size();

                outFile.write((char*)(&feature_tuple_string_length), sizeof(feature_tuple_string_length));
                outFile.write( &feature_tuple_string[0], feature_tuple_string_length);


                // <is_1_metablock, is_2_metablock, rel_orientation_1, rel_orientation_2, separation>
                current_paired_coords = make_tuple( get<1>(*it_coords), get<3>(*it_coords), get<4>(*it_coords), get<5>(*it_coords), get<6>(*it_coords) );
                outFile.write((char*) (&current_paired_coords), sizeof(current_paired_coords));
            }
        }

        outFile.close();
    }

    cout<<"save_partitioned_features ended: "<<start_assembly_idx_str<<" "<<end_assembly_idx_str<<"\n";
}


void locate_and_save_first_occurrence_features( string metablock_dir, unsigned long group_count, string assembly_blocks_dir,
                                                unsigned long start_assembly_idx, unsigned long end_assembly_idx, vector<ulli> count_offset,
                                                string partition_idx_str, map<ulli, string> metablock_tuple_idx_map, string filtered_feature_outdir){
    cout<<"locate_and_save_first_occurrence_features started: "<<metablock_dir<<" "<<group_count<<" "<<count_offset.size()<<" "<<start_assembly_idx<<"\n";
    string filename, metablock_tuple_string;
    unsigned long assembly_idx;

    list< list< tup_uubu > > first_occurrence_coords_outerlist;
    list< list< tup_uubu > >::iterator it_first_occurrence_coords_outerlist;
    list< tup_uubu >::iterator it_coords;

    for(assembly_idx=start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++){
        list< tup_uubu > coords_list;
        first_occurrence_coords_outerlist.push_back(coords_list);
    }

    ifstream inFile;
    ulli metablock_count;
    tup_uubu current_coords;

    for(unsigned long group_no=0; group_no<group_count; group_no++){
        it_first_occurrence_coords_outerlist = first_occurrence_coords_outerlist.begin();
        filename = (metablock_dir + to_string(group_no) + "_fo_" + partition_idx_str);
        inFile.open( filename.c_str(), ios::binary);

        for(assembly_idx = start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++, it_first_occurrence_coords_outerlist++){
            inFile.read((char*) (&metablock_count), sizeof(metablock_count));

            while(metablock_count--){
                inFile.read((char*) (&current_coords), sizeof(current_coords));
                get<3>(current_coords) += count_offset[group_no];
                (*it_first_occurrence_coords_outerlist).push_back(current_coords);
            }
        }

        inFile.close();
        remove(filename.c_str());
    }

    ofstream outFile;
    ulli metablock_idx, count;
    short metablock_tuple_string_length;
    tup_uub coords;
    
    assembly_idx = start_assembly_idx;
    for(it_first_occurrence_coords_outerlist = first_occurrence_coords_outerlist.begin();
            it_first_occurrence_coords_outerlist != first_occurrence_coords_outerlist.end() && assembly_idx<=end_assembly_idx;
            it_first_occurrence_coords_outerlist++, assembly_idx++){
        filename = (metablock_dir + "fo_" + to_string(assembly_idx));
        outFile.open(filename.c_str(), ios::binary);

        metablock_count = (*it_first_occurrence_coords_outerlist).size();
        outFile.write((char*) (&metablock_count), sizeof(metablock_count));
        count = 0;

        for(it_coords = (*it_first_occurrence_coords_outerlist).begin(); it_coords != (*it_first_occurrence_coords_outerlist).end(); it_coords++){

            coords = tup_uub( get<0>( *it_coords), get<1>( *it_coords), get<2>( *it_coords));
            metablock_idx = get<3>( *it_coords);
            metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];

            metablock_tuple_string_length = metablock_tuple_string.size();

            outFile.write((char*)(&metablock_tuple_string_length), sizeof(metablock_tuple_string_length));
            outFile.write( &metablock_tuple_string[0], metablock_tuple_string_length);
            outFile.write((char*) (&coords), sizeof(coords));

            if(count++%100==0){
                cout<<metablock_tuple_string<<" "<<metablock_tuple_string_length<<" (";
                cout<< get<0>( *it_coords) << "," << get<1>( *it_coords) << "," << get<2>( *it_coords) << ") ";
            }
        }
        cout<<"\n";
        outFile.close();
    }

    cout<<"locate_and_save_first_occurrence_features ended: "<<metablock_dir<<" "<<group_count<<" "<<count_offset.size()<<" "<<start_assembly_idx<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_feature_filtering_and_pairing.cpp "<<argc<<"\n";
    unsigned long group_count = stoul(argv[1]);
    string metablock_dir = argv[2];
    unsigned long min_sv_size = stoul(argv[3]);
    string assembly_blocks_dir = argv[4];
    string partition_idx_str = argv[5];
    string start_assembly_idx_str = argv[6];
    unsigned long start_assembly_idx = stoul(start_assembly_idx_str);
    string end_assembly_idx_str = argv[7];
    unsigned long end_assembly_idx = stoul(end_assembly_idx_str);
    ulli max_inter_feature_separation = strtoulli(argv[8]); 
    string contig_len_dir  = argv[9];
    short feature_partitions = static_cast<short>(stoi(argv[10]));
    ulli total_block_count = strtoulli(argv[11]);
    ulli total_metablock_count = strtoulli(argv[12]);
    string filtered_feature_outdir = argv[13];
    string paired_feature_outdir = argv[14];
    bool use_oriented_links = strncmp(argv[15],"1",1)==0;
    string oriented_links_dir;

    if(use_oriented_links){
        oriented_links_dir = argv[16];
    }

    ulli blocks_per_partition = ceil(total_block_count*1.0/feature_partitions);
    blocks_per_partition = ceil(blocks_per_partition/10.0)*10;
    ulli metablocks_per_partition = ceil(total_metablock_count*1.0/feature_partitions);
    metablocks_per_partition = ceil(metablocks_per_partition/10.0)*10;

    unsigned long assembly_idx;
    vector<ulli> count_offset, max_metablock_idx;

    map<ulli, string> metablock_tuple_idx_map = get_metablock_tuple_idx_map(metablock_dir, group_count, count_offset);

    // map< unsigned long, map< puiui, tup_bbffs > > mapped_oriented_links_map;
    // map< puiui, tup_bbffs > oriented_links_map;

    // contig ends map: contig_idx -> [ (distance_from_start, feature_1_idx, is_metablock, strand),
    //                                  (distance_from_end, feature_2_idx, is_metablock, strand) ]
    // list< map< unsigned int, ptupuubbtupuubb> > contig_ends_feature_map_list;

    for(auto it = count_offset.begin(); it != count_offset.end(); it++)
        cout<< *it << " ";
    cout<<"\n";

    for(assembly_idx=start_assembly_idx; assembly_idx<=end_assembly_idx; assembly_idx++){
        /*if(use_oriented_links){
            oriented_links_map.clear();
            oriented_links_map = load_oriented_links_map(oriented_links_dir + to_string(assembly_idx));
            mapped_oriented_links_map[ *it_assembly_innerlist ] = oriented_links_map;
        }*/
        max_metablock_idx.push_back(0);
    }

    list< list< list<tup_uubub> > > feature_partition_outerlist;
    list< list< list<tup_ububbbl> > > paired_feature_partition_outerlist;

    init_partitioned_feature_lists( feature_partition_outerlist, paired_feature_partition_outerlist,
                                    feature_partitions, end_assembly_idx - start_assembly_idx + 1);
    cout << feature_partition_outerlist.size() << " " << feature_partition_outerlist.front().size() << " ";
    cout << paired_feature_partition_outerlist.size() << " " << paired_feature_partition_outerlist.front().size() << "\n";

    // list< list< tup_uubub > > filtered_features_outerlist = get_filtered_features(metablock_dir, group_count, assembly_blocks_dir,
    //                                                                         start_assembly_idx, end_assembly_idx, count_offset, partition_idx_str,
    //                                                                         min_sv_size, contig_len_dir, max_metablock_idx);

    generate_filtered_features( metablock_dir, group_count, assembly_blocks_dir, start_assembly_idx, end_assembly_idx, count_offset,
                                partition_idx_str, min_sv_size, feature_partition_outerlist, paired_feature_partition_outerlist,
                                blocks_per_partition, metablocks_per_partition, feature_partitions, max_inter_feature_separation,
                                use_oriented_links, oriented_links_dir, contig_len_dir, max_metablock_idx);

    for(auto it=max_metablock_idx.begin(); it!=max_metablock_idx.end(); it++)
        cout<<metablock_tuple_idx_map[*it]<<" ";
    cout<<"\n";

    save_partitioned_features(  feature_partition_outerlist, paired_feature_partition_outerlist, metablock_tuple_idx_map,
                                start_assembly_idx_str, end_assembly_idx_str, filtered_feature_outdir, paired_feature_outdir);

    locate_and_save_first_occurrence_features(metablock_dir, group_count, assembly_blocks_dir, start_assembly_idx, end_assembly_idx,
                                            count_offset, partition_idx_str, metablock_tuple_idx_map, filtered_feature_outdir);

    return 0;
}