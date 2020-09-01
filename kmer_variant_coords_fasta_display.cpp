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


void fetch_variant_fasta(string assembly_fasta_file, ulli variant_start, ulli variant_end){
    // cout<<"fetch_variant_fasta started: "<<assembly_fasta_file<<" "<<variant_start<<" "<<variant_end<<"\n";

    string header, contig_header, current_fasta = "";
    list<string>::iterator it_assembly_fasta_lines;
    ulli pos=1, local_start, local_end, local_coords=1;
    short read_line_len, loc_on_line=0;
    bool found=false, being_read=false, continue_on_current_line=false;

    ifstream inFile;

    // Load assembly fasta sequence
    // Note that coords start from location index 1 (and not 0)
    list<string> assembly_fasta_lines = get_lines( assembly_fasta_file.c_str() );

    for(it_assembly_fasta_lines = assembly_fasta_lines.begin(); !found && it_assembly_fasta_lines != assembly_fasta_lines.end();
            it_assembly_fasta_lines++){

        if((*it_assembly_fasta_lines)[0] == '>'){
            loc_on_line = 0; // points to the current location on the line (0 indexed) to be read
            if(being_read){
                cout<<"ALERT!!! ERROR IN PARSING AND LOCATING VARIANTS: OVERFLOWING CONTIG BOUNDARIES!!! "<< (*it_assembly_fasta_lines);
                cout<<" "<<current_fasta<<"\n";
            }
            local_coords = 1; // reset local_coords to 1 at the start of the contig
            contig_header = *it_assembly_fasta_lines;
            continue;
        }

        read_line_len = (*it_assembly_fasta_lines).size();

        if( (pos + read_line_len - 1) < variant_start){
            // loc_on_line = 0
            pos += (read_line_len);
            local_coords += (read_line_len);
            continue;
        }

        if(!being_read){
            current_fasta = (*it_assembly_fasta_lines).substr(loc_on_line + variant_start - pos);
            loc_on_line += variant_start - pos;
            local_coords += variant_start - pos;
            pos = variant_start;
            local_start = local_coords;

            if(pos + read_line_len - 1 - loc_on_line < variant_end){
                // variant sequence is to be continued onto the next line
                pos += (read_line_len - loc_on_line); // start location for next line
                local_coords += (read_line_len - loc_on_line);
                being_read = true;
                loc_on_line = 0;
                continue;
            }
            else{
                // variant end on the current line
                current_fasta = current_fasta.substr(0, variant_end - variant_start + 1);
                being_read = false;
                found = true;
                local_end = local_start + variant_end - variant_start;
                continue;
            }
        }
        else{
            // variant is being read from previous line: loc_on_line:0, pos: start location of current line

            if(pos + read_line_len - 1 < variant_end){
                // variant sequence is to be continued onto the next line
                current_fasta += (*it_assembly_fasta_lines);
                pos += (read_line_len); // start location for next line
                local_coords += read_line_len;
                being_read = true;
                loc_on_line = 0;
                continue;
            }
            else{
                // variant end on the current line
                current_fasta += (*it_assembly_fasta_lines).substr(0, variant_end - pos + 1);
                being_read = false;
                found = true;
                local_end = local_start + variant_end - variant_start;
                continue;
            }
        }
    }

    header = to_string(variant_start) + "-" + to_string(variant_end) + "::" + contig_header + ":" + to_string(local_start) + "-" + to_string(local_end);
    cout<< generate_formatted_fasta(header, current_fasta) <<"\n\n";
    // cout<<"fetch_variant_fasta ended: "<<assembly_fasta_file<<" "<<variant_start<<" "<<variant_end<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_variant_coords_fasta_display.cpp "<<argc<<"\n";

    string all_assembly_fasta_list_file = argv[1];
    string assembly_idx_str = argv[2];
    unsigned long assembly_idx = stoul(assembly_idx_str);
    ulli variant_start = strtoulli(argv[3]);
    ulli variant_end = strtoulli(argv[4]);

    // list<string> fasta_filepaths = get_lines(assembly_fasta_list_file.c_str());

    string assembly_fasta_file = "";
    unsigned long idx;

    list<string> fasta_filepaths = get_lines(all_assembly_fasta_list_file.c_str());
    list<string>::iterator it_assembly_fasta = fasta_filepaths.begin();

    for(idx=0; it_assembly_fasta != fasta_filepaths.end() && idx!=assembly_idx; idx++){
        it_assembly_fasta++;
    }
    if(idx==assembly_idx)
        assembly_fasta_file = *it_assembly_fasta;
    else{
        cout<<"ERROR!!! ASSEMBLY FILE NOT FOUND!!! "<<assembly_idx<<" "<<idx<<"\n";
        return -1;
    }

    cout<<assembly_fasta_file<<"\n";

    fetch_variant_fasta(assembly_fasta_file, variant_start, variant_end);

    return 0;
}