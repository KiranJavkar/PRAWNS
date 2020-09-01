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


list<ulli> get_block_indices(const char* FILENAME){
    ifstream input_file(FILENAME);
    list<ulli> block_idx_list;
    string text;
    while (getline( input_file, text )){
        block_idx_list.push_back(strtoulli(text));
    }
    return block_idx_list;
}


string generate_formatted_fasta(string header, string fasta_sequence, short line_width=60){
    ulli i, seq_length = fasta_sequence.size();
    string formatted_string = ">" + header + " " + "len=" + to_string(seq_length) + "\n";
    for(i=0; i<seq_length; i+=line_width){
        if(i+line_width < seq_length)
            formatted_string += fasta_sequence.substr(i, line_width) + "\n";
        else
            formatted_string += fasta_sequence.substr(i);
    }
    return formatted_string;
}


void save_fasta_from_buffers(list<string>& buffer, bool is_metablock, string fasta_outdir, string batch_no_str){
    string opstr = "";
    for(list<string>::iterator it = buffer.begin(); it != buffer.end(); it++)
        opstr += (*it) + "\n";

    ofstream outFile;
    if(is_metablock)
        outFile.open(fasta_outdir + "metablocks_" + batch_no_str + ".fasta", ios::app);
    else
        outFile.open(fasta_outdir + "blocks_" + batch_no_str + ".fasta", ios::app);

    outFile << opstr;
    outFile.close();

    buffer.clear();
}


// void extract_and_save_variant_fasta(unsigned long start_assembly_idx, list<string> fasta_filepaths, string metablocks_dir,
//                                     list<ulli> filtered_block_idx_list, string assembly_blocks_dir, unsigned long seq2write_batch,
//                                     string fasta_outdir, string batch_no_str){
//     cout<<"extract_and_save_variant_fasta started: "<<start_assembly_idx<<" "<<fasta_filepaths.size()<<" "<<batch_no_str;
//     cout<<" "<<filtered_block_idx_list.size()<<"\n";

//     string last_sequence, current_fasta, filename, metablock_tuple_string;
//     list<string> buffer_metablock, buffer_block, assembly_fasta_lines;
//     ulli prev_start, prev_end, current_start, current_end, pos, block_count, block_idx, metablock_count, metablock_idx;
//     list<string>::iterator it_assembly_fasta = fasta_filepaths.begin();
//     list<string>::iterator it_assembly_fasta_lines;
//     unsigned long assembly_idx = start_assembly_idx; //, end_assembly_idx = (start_assembly_idx + fasta_filepaths.size() - 1);
//     short count, metablock_tuple_string_length, read_line_len, loc_on_line;
//     tup_uub coords;
//     list<tup_uubub> feature_coords_list; // <start, end, strand, feature_idx, is_metablock>
//     // in case of metablock, the feature_idx is a placeholder mapping to the corresponding metablock tuple
//     map<ulli, string> metablock_tuple_idx_map;
//     list<ulli>::iterator it_filtered_block_idx_list;
//     list<tup_uubub>::iterator it_features;
//     bool iterate, being_read, continue_on_current_line;

//     ifstream inFile;

//     while(it_assembly_fasta != fasta_filepaths.end()){

//         buffer_metablock.clear();
//         buffer_block.clear();
//         feature_coords_list.clear();
//         metablock_tuple_idx_map.clear();
//         prev_start = 0;
//         prev_end = 0;
//         pos = 1;
//         loc_on_line = 0;

//         // Load all long blocks' coords with first occurrences in current assembly
//         // The merged block idx list as well as first occurrences list are sorted on block idx by design

//         filename = assembly_blocks_dir + "fo_" + to_string(assembly_idx);
//         cout<<"\t"<<filename;
//         inFile.open( filename.c_str(), ios::binary);
//         inFile.read((char*) (&block_count), sizeof(block_count));
//         cout<<" "<<block_count;
//         it_filtered_block_idx_list = filtered_block_idx_list.begin();
//         iterate = (it_filtered_block_idx_list != filtered_block_idx_list.end());

//         while(iterate && block_count--){
//             inFile.read((char*) (&block_idx), sizeof(block_idx));
//             inFile.read((char*) (&coords), sizeof(coords));
//             while(it_filtered_block_idx_list != filtered_block_idx_list.end() && (*it_filtered_block_idx_list) < block_idx)
//                 it_filtered_block_idx_list++;
//             if(it_filtered_block_idx_list == filtered_block_idx_list.end())
//                 iterate = false;
//             else{
//                 if((*it_filtered_block_idx_list)==block_idx){
//                     feature_coords_list.push_back( tup_uubub(get<0>(coords), get<1>(coords), get<2>(coords), block_idx, false) );
//                     it_filtered_block_idx_list++;
//                 }
//             }
//         }

//         inFile.close();
//         // remove(filename.c_str());

//         cout<<"\tFiltered blocks located: "<<assembly_idx<<" "<<feature_coords_list.size()<<"\n";

//         count = 0;
//         for(it_features=feature_coords_list.begin(); it_features!=feature_coords_list.end() && count<10; it_features++, count++){
//             cout<<"("<< get<0>(*it_features) << ","<< get<1>(*it_features) << ","<< get<2>(*it_features) << ",";
//             cout<< get<3>(*it_features) << ","<< get<4>(*it_features) << ") ";
//         }
//         cout<<"\n";


//         // Load all metablocks' coords with first occurrences in current assembly

//         filename = metablocks_dir + "fo_" + to_string(assembly_idx);
//         cout<<"\t"<<filename;
//         inFile.open( filename.c_str(), ios::binary);
//         inFile.read((char*) (&metablock_count), sizeof(metablock_count));
//         cout<<" "<<metablock_count;
//         metablock_idx = 0;

//         while(metablock_count--){

//             inFile.read((char*)(&metablock_tuple_string_length), sizeof(metablock_tuple_string_length));
//             metablock_tuple_string.resize(metablock_tuple_string_length);
//             inFile.read( &metablock_tuple_string[0], metablock_tuple_string_length);
//             inFile.read((char*) (&coords), sizeof(coords));

//             feature_coords_list.push_back( tup_uubub(get<0>(coords), get<1>(coords), get<2>(coords), metablock_idx, true) );
//             metablock_tuple_idx_map[metablock_idx++] = metablock_tuple_string;
//         }

//         inFile.close();
//         // remove(filename.c_str());

//         cout<<"\tMetablocks located and appended: "<<assembly_idx<<" "<<feature_coords_list.size()<<"\n";

//         // Sort the current assembly first occurrence features by their start positions
//         feature_coords_list.sort();

//         count = 0;
//         for(it_features=feature_coords_list.begin(); it_features!=feature_coords_list.end() && count<10; it_features++, count++){
//             cout<<"("<< get<0>(*it_features) << ","<< get<1>(*it_features) << ","<< get<2>(*it_features) << ",";
//             cout<< get<3>(*it_features) << ","<< get<4>(*it_features) << ") ";
//         }
//         cout<<"\n";


//         // Load assembly fasta sequence
//         // Note that coords start from location index 1 (and not 0)
//         assembly_fasta_lines = get_lines( (*it_assembly_fasta).c_str() );
//         it_features = feature_coords_list.begin();
//         iterate = (it_features != feature_coords_list.end());
//         being_read = false;
//         block_count = 0;
//         metablock_count = 0;

//         for(it_assembly_fasta_lines = assembly_fasta_lines.begin();
//                 it_features != feature_coords_list.end() && it_assembly_fasta_lines != assembly_fasta_lines.end();
//                 /*it_features++, it_assembly_fasta_lines++*/){

//             // Write buffers
//             if(metablock_count == seq2write_batch){
//                 save_fasta_from_buffers(buffer_metablock, true, fasta_outdir, batch_no_str);
//                 metablock_count = 0;
//             }
//             if(block_count == seq2write_batch){
//                 save_fasta_from_buffers(buffer_block, false, fasta_outdir, batch_no_str);
//                 block_count = 0;
//             }

//             current_start = get<0>(*it_features);
//             current_end = get<1>(*it_features);

//             if((*it_assembly_fasta_lines)[0] == '>'){
//                 last_sequence = "";
//                 prev_end = pos-1; // pos: next location to read on assembly (1 indexed)
//                 prev_start = prev_end;
//                 loc_on_line = 0; // points to the current location on the line (0 indexed) to be read
//                 if(being_read){
//                     cout<<"ALERT!!! ERROR IN PARSING AND LOCATING VARIANTS: OVERFLOWING CONTIG BOUNDARIES!!! "<< (*it_assembly_fasta_lines);
//                     cout<<" "<<current_fasta<<"\n";
//                 }
//                 it_assembly_fasta_lines++;
//                 continue;
//             }

//             // being_read = current_start<pos;
//             read_line_len = (*it_assembly_fasta_lines).size();

//             if( (pos + read_line_len - 1 - loc_on_line) < current_start){
//                 pos += (read_line_len - loc_on_line);
//                 it_assembly_fasta_lines++;
//                 loc_on_line = 0;
//                 continue;
//             }

//             cout<< *it_assembly_fasta_lines <<" "<< current_start<<" "<<current_end<<" "<<loc_on_line<<" "<<pos<<" ";
//             cout<<being_read<<" "<<prev_start<<" "<<prev_end<<" "<<buffer_metablock.size()<<" "<<buffer_block.size()<<"\n";

//             if(!being_read){

//                 if(current_start < pos){
//                     // current variant overlaps previous variant
//                     // (can be contained inside previous one as well since a block - that is a part of a metablock - is a variant to be located)
//                     current_fasta = last_sequence.substr(current_start - prev_start);

//                     if(current_end <= prev_end){
//                         // variant is contained inside the previous variant located
//                         // ** DO NOT UPDATE THE SAVED COORDS AND SEQUENCE FOR PREVIOUS VARIANT after saving the current one
//                         current_fasta = current_fasta.substr(0, current_end - current_start + 1);
//                         // Add Fasta to the corresponding buffer
//                         if(get<4>(*it_features)){
//                             // metablock
//                             metablock_idx = get<3>(*it_features);
//                             metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];
//                             buffer_metablock.push_back( generate_formatted_fasta(metablock_tuple_string + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
//                             metablock_count++;
//                         }
//                         else{
//                             block_idx = get<3>(*it_features);
//                             buffer_block.push_back( generate_formatted_fasta( to_string(block_idx) + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
//                             block_count++;
//                         }
//                         it_features++;
//                         continue;
//                     }
//                     else{
//                         // being_read = true;
//                         current_fasta += (*it_assembly_fasta_lines).substr(loc_on_line);
//                     }
//                 }
//                 else{
//                     current_fasta = (*it_assembly_fasta_lines).substr(loc_on_line + current_start - pos);
//                     loc_on_line += current_start - pos;
//                     pos = current_start;
//                 }

//                 if(pos + read_line_len - 1 - loc_on_line < current_end){
//                     // variant sequence is to be continued onto the next line
//                     pos += (read_line_len - loc_on_line); // start location for next line
//                     being_read = true;
//                     it_assembly_fasta_lines++;
//                     loc_on_line = 0;
//                     continue;
//                 }
//                 else{
//                     // variant end on the current line
//                     current_fasta = current_fasta.substr(0, current_end - current_start + 1);
//                     loc_on_line += current_end - pos + 1; // next location to read on line
//                     pos = current_end + 1; // next location to read on assembly
//                     being_read = false;

//                     if(get<4>(*it_features)){
//                         // metablock
//                         metablock_idx = get<3>(*it_features);
//                         metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];
//                         // Add to buffer
//                         buffer_metablock.push_back( generate_formatted_fasta(metablock_tuple_string + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
//                         metablock_count++;
//                     }
//                     else{
//                         block_idx = get<3>(*it_features);
//                         buffer_block.push_back( generate_formatted_fasta( to_string(block_idx) + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
//                         block_count++;
//                     }
//                     // Set the variant fasta sequence to last_sequence
//                     last_sequence = current_fasta;
//                     prev_start = current_start;
//                     prev_end = current_end;

//                     it_features++;
//                 }
//             }
//             else{
//                 // variant is being read from previous line: loc_on_line:0, pos: start location of current line

//                 if(pos + read_line_len - 1 - loc_on_line < current_end){
//                     // variant sequence is to be continued onto the next line
//                     current_fasta += (*it_assembly_fasta_lines);
//                     pos += (read_line_len - loc_on_line); // start location for next line
//                     being_read = true;
//                     it_assembly_fasta_lines++;
//                     loc_on_line = 0;
//                     continue;
//                 }
//                 else{
//                     // variant end on the current line
//                     current_fasta += (*it_assembly_fasta_lines).substr(0, current_end - pos + 1);
//                     loc_on_line += current_end - pos + 1; // next location to read on line
//                     pos = current_end + 1; // next location to read on assembly
//                     being_read = false;

//                     if(get<4>(*it_features)){
//                         // metablock
//                         metablock_idx = get<3>(*it_features);
//                         metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];
//                         // Add to buffer
//                         buffer_metablock.push_back( generate_formatted_fasta(metablock_tuple_string + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
//                         metablock_count++;
//                     }
//                     else{
//                         block_idx = get<3>(*it_features);
//                         buffer_block.push_back( generate_formatted_fasta( to_string(block_idx) + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
//                         block_count++;
//                     }
//                     // Set the variant fasta sequence to last_sequence
//                     last_sequence = current_fasta;
//                     prev_start = current_start;
//                     prev_end = current_end;

//                     it_features++;
//                 }
//             }

//             // pos++;
//         }

//         if(metablock_count > 0)
//             save_fasta_from_buffers(buffer_metablock, true, fasta_outdir, batch_no_str);
//         if(block_count > 0)
//             save_fasta_from_buffers(buffer_block, false, fasta_outdir, batch_no_str);

//         cout<<"\tFASTA sequences extracted and saved: "<<assembly_idx<<"\n\n";

//         it_assembly_fasta++;
//         assembly_idx++;
//     }
//     cout<<"extract_and_save_variant_fasta ended: "<<start_assembly_idx<<" "<<fasta_filepaths.size()<<" "<<batch_no_str<<"\n";
// }


// start_assembly_idx: batch_no; and assembly_idx is incremented by batch_no
//      (consecutinve entries in fasta_filepaths correspond to assemblies with their assembly_idx separated by batch_no intermediate entries)

void extract_and_save_variant_fasta2(short batch_no, list<string> fasta_filepaths, string metablocks_dir,
                                    list<ulli> filtered_block_idx_list, string assembly_blocks_dir, unsigned long seq2write_batch,
                                    string fasta_outdir, short batch_count, string batch_no_str){
    cout<<"extract_and_save_variant_fasta2 started: "<<batch_no<<" "<<fasta_filepaths.size()<<" "<<batch_no_str;
    cout<<" "<<filtered_block_idx_list.size()<<"\n";

    string last_sequence, current_fasta, filename, metablock_tuple_string;
    list<string> buffer_metablock, buffer_block, assembly_fasta_lines;
    ulli prev_start, prev_end, current_start, current_end, pos, block_count, block_idx, metablock_count, metablock_idx;
    list<string>::iterator it_assembly_fasta = fasta_filepaths.begin();
    list<string>::iterator it_assembly_fasta_lines;
    unsigned long start_assembly_idx = batch_no;
    unsigned long assembly_idx = start_assembly_idx; //, end_assembly_idx = (start_assembly_idx + fasta_filepaths.size() - 1);
    short count, metablock_tuple_string_length, read_line_len, loc_on_line;
    tup_uub coords;
    list<tup_uubub> feature_coords_list; // <start, end, strand, feature_idx, is_metablock>
    // in case of metablock, the feature_idx is a placeholder mapping to the corresponding metablock tuple
    map<ulli, string> metablock_tuple_idx_map;
    list<ulli>::iterator it_filtered_block_idx_list;
    list<tup_uubub>::iterator it_features;
    bool iterate, being_read, continue_on_current_line;

    ifstream inFile;

    while(it_assembly_fasta != fasta_filepaths.end()){

        buffer_metablock.clear();
        buffer_block.clear();
        feature_coords_list.clear();
        metablock_tuple_idx_map.clear();
        prev_start = 0;
        prev_end = 0;
        pos = 1;
        loc_on_line = 0;

        // Load all long blocks' coords with first occurrences in current assembly
        // The merged block idx list as well as first occurrences list are sorted on block idx by design

        filename = assembly_blocks_dir + "fo_" + to_string(assembly_idx);
        // cout<<"\t"<<filename;
        inFile.open( filename.c_str(), ios::binary);
        if(!inFile.fail()){

            inFile.read((char*) (&block_count), sizeof(block_count));
            // cout<<" "<<block_count;
            it_filtered_block_idx_list = filtered_block_idx_list.begin();
            iterate = (it_filtered_block_idx_list != filtered_block_idx_list.end());

            while(iterate && block_count--){
                inFile.read((char*) (&block_idx), sizeof(block_idx));
                inFile.read((char*) (&coords), sizeof(coords));
                while(it_filtered_block_idx_list != filtered_block_idx_list.end() && (*it_filtered_block_idx_list) < block_idx)
                    it_filtered_block_idx_list++;
                if(it_filtered_block_idx_list == filtered_block_idx_list.end())
                    iterate = false;
                else{
                    if((*it_filtered_block_idx_list)==block_idx){
                        feature_coords_list.push_back( tup_uubub(get<0>(coords), get<1>(coords), get<2>(coords), block_idx, false) );
                        it_filtered_block_idx_list++;
                    }
                }
            }
        }

        inFile.close();
        // remove(filename.c_str());

        cout<<"\tFiltered blocks located: "<<assembly_idx<<" "<<feature_coords_list.size()<<"\n";

        // count = 0;
        // for(it_features=feature_coords_list.begin(); it_features!=feature_coords_list.end() && count<10; it_features++, count++){
        //     cout<<"("<< get<0>(*it_features) << ","<< get<1>(*it_features) << ","<< get<2>(*it_features) << ",";
        //     cout<< get<3>(*it_features) << ","<< get<4>(*it_features) << ") ";
        // }
        // cout<<"\n";


        // Load all metablocks' coords with first occurrences in current assembly

        filename = metablocks_dir + "fo_" + to_string(assembly_idx);
        // cout<<"\t"<<filename;
        inFile.open( filename.c_str(), ios::binary);
        if(!inFile.fail()){
            
            inFile.read((char*) (&metablock_count), sizeof(metablock_count));
            // cout<<" "<<metablock_count;
            metablock_idx = 0;

            while(metablock_count--){

                inFile.read((char*)(&metablock_tuple_string_length), sizeof(metablock_tuple_string_length));
                metablock_tuple_string.resize(metablock_tuple_string_length);
                inFile.read( &metablock_tuple_string[0], metablock_tuple_string_length);
                inFile.read((char*) (&coords), sizeof(coords));

                feature_coords_list.push_back( tup_uubub(get<0>(coords), get<1>(coords), get<2>(coords), metablock_idx, true) );
                metablock_tuple_idx_map[metablock_idx++] = metablock_tuple_string;
            }
        }

        inFile.close();
        remove(filename.c_str());

        cout<<"\tMetablocks located and appended: "<<assembly_idx<<" "<<feature_coords_list.size()<<"\n";

        // Sort the current assembly first occurrence features by their start positions
        feature_coords_list.sort();

        // count = 0;
        // for(it_features=feature_coords_list.begin(); it_features!=feature_coords_list.end() && count<10; it_features++, count++){
        //     cout<<"("<< get<0>(*it_features) << ","<< get<1>(*it_features) << ","<< get<2>(*it_features) << ",";
        //     cout<< get<3>(*it_features) << ","<< get<4>(*it_features) << ") ";
        // }
        // cout<<"\n";


        // Load assembly fasta sequence
        // Note that coords start from location index 1 (and not 0)
        assembly_fasta_lines = get_lines( (*it_assembly_fasta).c_str() );
        it_features = feature_coords_list.begin();
        iterate = (it_features != feature_coords_list.end());
        being_read = false;
        block_count = 0;
        metablock_count = 0;

        for(it_assembly_fasta_lines = assembly_fasta_lines.begin();
                it_features != feature_coords_list.end() && it_assembly_fasta_lines != assembly_fasta_lines.end();
                /*it_features++, it_assembly_fasta_lines++*/){

            // Write buffers
            if(metablock_count == seq2write_batch){
                save_fasta_from_buffers(buffer_metablock, true, fasta_outdir, batch_no_str);
                metablock_count = 0;
            }
            if(block_count == seq2write_batch){
                save_fasta_from_buffers(buffer_block, false, fasta_outdir, batch_no_str);
                block_count = 0;
            }

            current_start = get<0>(*it_features);
            current_end = get<1>(*it_features);

            if((*it_assembly_fasta_lines)[0] == '>'){
                last_sequence = "";
                prev_end = pos-1; // pos: next location to read on assembly (1 indexed)
                prev_start = prev_end;
                loc_on_line = 0; // points to the current location on the line (0 indexed) to be read
                if(being_read){
                    cout<<"ALERT!!! ERROR IN PARSING AND LOCATING VARIANTS: OVERFLOWING CONTIG BOUNDARIES!!! "<< (*it_assembly_fasta_lines);
                    cout<<" "<<current_fasta<<"\n";
                }
                it_assembly_fasta_lines++;
                continue;
            }

            // being_read = current_start<pos;
            read_line_len = (*it_assembly_fasta_lines).size();

            if( (pos + read_line_len - 1 - loc_on_line) < current_start){
                pos += (read_line_len - loc_on_line);
                it_assembly_fasta_lines++;
                loc_on_line = 0;
                continue;
            }

            // cout<< *it_assembly_fasta_lines <<" "<< current_start<<" "<<current_end<<" "<<loc_on_line<<" "<<pos<<" ";
            // cout<<being_read<<" "<<prev_start<<" "<<prev_end<<" "<<buffer_metablock.size()<<" "<<buffer_block.size()<<"\n";

            if(!being_read){

                if(current_start < pos){
                    // current variant overlaps previous variant
                    // (can be contained inside previous one as well since a block - that is a part of a metablock - is a variant to be located)
                    current_fasta = last_sequence.substr(current_start - prev_start);

                    if(current_end <= prev_end){
                        // variant is contained inside the previous variant located
                        // ** DO NOT UPDATE THE SAVED COORDS AND SEQUENCE FOR PREVIOUS VARIANT after saving the current one
                        current_fasta = current_fasta.substr(0, current_end - current_start + 1);
                        // Add Fasta to the corresponding buffer
                        if(get<4>(*it_features)){
                            // metablock
                            metablock_idx = get<3>(*it_features);
                            metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];
                            buffer_metablock.push_back( generate_formatted_fasta(metablock_tuple_string + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
                            metablock_count++;
                        }
                        else{
                            block_idx = get<3>(*it_features);
                            buffer_block.push_back( generate_formatted_fasta( to_string(block_idx) + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
                            block_count++;
                        }
                        it_features++;
                        continue;
                    }
                    else{
                        // being_read = true;
                        current_fasta += (*it_assembly_fasta_lines).substr(loc_on_line);
                    }
                }
                else{
                    current_fasta = (*it_assembly_fasta_lines).substr(loc_on_line + current_start - pos);
                    loc_on_line += current_start - pos;
                    pos = current_start;
                }

                if(pos + read_line_len - 1 - loc_on_line < current_end){
                    // variant sequence is to be continued onto the next line
                    pos += (read_line_len - loc_on_line); // start location for next line
                    being_read = true;
                    it_assembly_fasta_lines++;
                    loc_on_line = 0;
                    continue;
                }
                else{
                    // variant end on the current line
                    current_fasta = current_fasta.substr(0, current_end - current_start + 1);
                    loc_on_line += current_end - pos + 1; // next location to read on line
                    pos = current_end + 1; // next location to read on assembly
                    being_read = false;

                    if(get<4>(*it_features)){
                        // metablock
                        metablock_idx = get<3>(*it_features);
                        metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];
                        // Add to buffer
                        buffer_metablock.push_back( generate_formatted_fasta(metablock_tuple_string + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
                        metablock_count++;
                    }
                    else{
                        block_idx = get<3>(*it_features);
                        buffer_block.push_back( generate_formatted_fasta( to_string(block_idx) + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
                        block_count++;
                    }
                    // Set the variant fasta sequence to last_sequence
                    last_sequence = current_fasta;
                    prev_start = current_start;
                    prev_end = current_end;

                    it_features++;
                }
            }
            else{
                // variant is being read from previous line: loc_on_line:0, pos: start location of current line

                if(pos + read_line_len - 1 - loc_on_line < current_end){
                    // variant sequence is to be continued onto the next line
                    current_fasta += (*it_assembly_fasta_lines);
                    pos += (read_line_len - loc_on_line); // start location for next line
                    being_read = true;
                    it_assembly_fasta_lines++;
                    loc_on_line = 0;
                    continue;
                }
                else{
                    // variant end on the current line
                    current_fasta += (*it_assembly_fasta_lines).substr(0, current_end - pos + 1);
                    loc_on_line += current_end - pos + 1; // next location to read on line
                    pos = current_end + 1; // next location to read on assembly
                    being_read = false;

                    if(get<4>(*it_features)){
                        // metablock
                        metablock_idx = get<3>(*it_features);
                        metablock_tuple_string = metablock_tuple_idx_map[metablock_idx];
                        // Add to buffer
                        buffer_metablock.push_back( generate_formatted_fasta(metablock_tuple_string + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
                        metablock_count++;
                    }
                    else{
                        block_idx = get<3>(*it_features);
                        buffer_block.push_back( generate_formatted_fasta( to_string(block_idx) + ":" + to_string(assembly_idx) + "(" + to_string(current_start) + "," + to_string(current_end) + ")", current_fasta) );
                        block_count++;
                    }
                    // Set the variant fasta sequence to last_sequence
                    last_sequence = current_fasta;
                    prev_start = current_start;
                    prev_end = current_end;

                    it_features++;
                }
            }

            // pos++;
        }

        if(metablock_count > 0)
            save_fasta_from_buffers(buffer_metablock, true, fasta_outdir, batch_no_str);
        if(block_count > 0)
            save_fasta_from_buffers(buffer_block, false, fasta_outdir, batch_no_str);

        cout<<"\tFASTA sequences extracted and saved: "<<assembly_idx<<"\n\n";

        it_assembly_fasta++;
        assembly_idx += batch_count;
    }
    cout<<"extract_and_save_variant_fasta2 ended: "<<start_assembly_idx<<" "<<fasta_filepaths.size()<<" "<<batch_no_str<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_feature_fasta_generator.cpp "<<argc<<"\n";
    // string assembly_fasta_list_file = argv[1];
    // unsigned long start_assembly_idx = stoul(argv[2]);
    // unsigned long seq2write_batch = stoul(argv[3]);
    // string metablocks_dir = argv[4];
    // string filtered_blocks_file = argv[5];
    // string assembly_blocks_dir = argv[6];
    // string fasta_outdir = argv[7];
    // string batch_no_str = argv[8];

    string all_assembly_fasta_list_file = argv[1];
    unsigned long seq2write_batch = stoul(argv[2]);
    string metablocks_dir = argv[3];
    string filtered_blocks_file = argv[4];
    string assembly_blocks_dir = argv[5];
    string fasta_outdir = argv[6];
    string batch_no_str = argv[7];
    short batch_no = static_cast<short>(stoi(batch_no_str));
    short batch_count = static_cast<short>(stoi(argv[8]));

    // list<string> fasta_filepaths = get_lines(assembly_fasta_list_file.c_str());

    list<string> fasta_filepaths = get_lines(all_assembly_fasta_list_file.c_str());
    list<string>::iterator it_assembly_fasta = fasta_filepaths.begin();
    for(unsigned long assembly_idx=0; it_assembly_fasta != fasta_filepaths.end(); assembly_idx++){
        if(assembly_idx%batch_count == batch_no)
            it_assembly_fasta++;
        else
            it_assembly_fasta = fasta_filepaths.erase(it_assembly_fasta);
    }
    cout<<"\t"<<batch_no<<": "<<fasta_filepaths.size()<<" assemblies\n";
    
    // for(auto it=fasta_filepaths.begin(); it!=fasta_filepaths.end(); it++)
    //     cout<< *it <<"\n";

    list<ulli> filtered_block_idx_list = get_block_indices(filtered_blocks_file.c_str()); //sorted by design

    // cout<<filtered_block_idx_list.size()<<"\n";
    // short count=0;
    // for(auto it = filtered_block_idx_list.begin(); count<10; count++, it++)
    //     cout<< *it << " ";
    // cout<<"\n";

    // remove(assembly_fasta_list_file.c_str());

    // extract_and_save_variant_fasta( start_assembly_idx, fasta_filepaths, metablocks_dir, filtered_block_idx_list, assembly_blocks_dir,
    //                                 seq2write_batch, fasta_outdir, batch_no_str);

    extract_and_save_variant_fasta2(batch_no, fasta_filepaths, metablocks_dir, filtered_block_idx_list, assembly_blocks_dir,
                                    seq2write_batch, fasta_outdir, batch_count, batch_no_str);

    return 0;
}