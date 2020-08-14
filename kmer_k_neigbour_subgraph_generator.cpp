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


unsigned long max (const unsigned long& a, const unsigned long& b) {
  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}


auto neighbour_pair_compare = [](ppuuul a, ppuuul b) { return a.second > b.second; };


void generate_k_neighbour_subgraph( string neighbour_pair_partition_dir, unsigned long assembly_count, short partition_idx, short k_neighbours,
                                    unsigned long min_presence_count, string neighbour_subgraph_outdir){
    cout<<"generate_k_neighbour_subgraph started: "<<neighbour_pair_partition_dir<<" "<<assembly_count<<" "<<partition_idx<<" ";
    cout<<k_neighbours<<" "<<neighbour_subgraph_outdir<<"\n";

    unsigned long assembly_idx;
    string filename, outstring="", file_prefix = neighbour_pair_partition_dir + to_string(partition_idx) + "_";
    ofstream outFile;
    list< ppuuul > pair_counts_list;
    list< ppuuul >::iterator it_pair_counts;
    list<string> lines, lsplit;
    list<string>::iterator it_coords;
    ulli block_1_idx, block_2_idx;

    // Load and merge neighbour pairs from all assemblies

    for(assembly_idx=0; assembly_idx<assembly_count; assembly_idx++){
        it_pair_counts = pair_counts_list.begin();
        filename = file_prefix + to_string(assembly_idx);
        lines = get_lines(filename.c_str());

        if(lines.size() == 0){
            cout << assembly_idx << " : NO_ENTRIES_FOUND!!!\n";
            continue;
        }

        lsplit = split(*(lines.begin()), ',');
        
        it_coords = lsplit.begin();
        while(it_coords!=lsplit.end()){
            block_1_idx = strtoulli(*it_coords);
            it_coords++;
            block_2_idx = strtoulli(*it_coords);
            it_coords++;

            while(it_pair_counts != pair_counts_list.end() && (it_pair_counts->first).first < block_1_idx)
                it_pair_counts++;

            while(it_pair_counts != pair_counts_list.end() && ((it_pair_counts->first).first == block_1_idx
                                                                && (it_pair_counts->first).second < block_2_idx))
                it_pair_counts++;

            if(it_pair_counts == pair_counts_list.end())
                pair_counts_list.push_back(make_pair( make_pair(block_1_idx, block_2_idx), 1));
            else if((it_pair_counts->first).first == block_1_idx && (it_pair_counts->first).second == block_2_idx)
                it_pair_counts->second += 1;
            else
                pair_counts_list.insert(it_pair_counts, make_pair( make_pair(block_1_idx, block_2_idx), 1) );
        }

        cout<<assembly_idx<<"\t"<<pair_counts_list.size()<<" : ";
        short count = 0;
        for(it_pair_counts = pair_counts_list.begin(); it_pair_counts != pair_counts_list.end() && count<5; it_pair_counts++, count++)
            cout<<"(("<<(it_pair_counts->first).first<<","<<(it_pair_counts->first).second<<"),"<<it_pair_counts->second<<")\t";
        cout<<"\n";
    }

    // Select k neighbours for each block which have the maximum paired occurence with the corresponding block
    priority_queue< ppuuul, vector< ppuuul >, decltype(neighbour_pair_compare) > current_block_neighbours(neighbour_pair_compare);
    bool rand_select_start = false;
    short rand_start, pos, excess;
    unsigned long current_max = 0, current_min_required = 0;

    for(it_pair_counts = pair_counts_list.begin(); it_pair_counts != pair_counts_list.end(); it_pair_counts++){
        if(!current_block_neighbours.empty()){
            if( ((current_block_neighbours.top()).first).first != (it_pair_counts->first).first){

                if(outstring != ""){
                    outstring.pop_back();
                    outstring += "\n";
                }
                outstring += to_string(((current_block_neighbours.top()).first).first) + ",";

                // // current_min_required = 0.95*current_max;
                // if(current_max>min_presence_count)
                //     current_min_required = (0.95*current_max > current_max-10)?0.95*current_max:current_max-10;
                // else
                //     current_min_required = min_presence_count;

                // current_min_required = 0.98*current_max;
                if(current_max>min_presence_count){
                    current_min_required = max( max(0.95*current_max, current_max-10), min_presence_count);
                    // current_min_required = max( max(0.90*current_max, current_max-10), min_presence_count);
                    // current_min_required = max( max(0.98*current_max, current_max-5), min_presence_count);
                }
                else
                    current_min_required = min_presence_count;

                // cout<<((current_block_neighbours.top()).first).first<<" "<<current_max<<" "<<current_min_required<<"\n";

                while(!current_block_neighbours.empty() && (current_block_neighbours.top()).second < current_min_required)
                    current_block_neighbours.pop();

                // Select random sample if more than k neighbours detected
                if(current_block_neighbours.size() > k_neighbours){
                    rand_select_start = true;
                    excess = current_block_neighbours.size() - k_neighbours + 1;
                    rand_start = rand() % excess;
                    pos = 0;
                }
                while(!current_block_neighbours.empty()){
                    // cout<<"(("<<((current_block_neighbours.top()).first).first<<","<<((current_block_neighbours.top()).first).second;
                    // cout<<"),"<<(current_block_neighbours.top()).second<<")";
                    // current_block_neighbours.pop();
                    // if(rand_select_start)
                    //     if(rand_start == pos || pos >= excess)
                    //         cout<<"**";
                    // cout<<"\t";
                    if(rand_select_start){
                        if(rand_start == pos){
                            if((current_block_neighbours.top()).second >= min_presence_count){
                                outstring += to_string(((current_block_neighbours.top()).first).second) + ",";
                                outstring += to_string(assembly_count - (current_block_neighbours.top()).second + 1) + ",";
                            }
                        }
                        pos++;
                        if(pos >= excess)
                            rand_select_start = false;
                    }
                    else
                        if((current_block_neighbours.top()).second >= min_presence_count){
                            outstring += to_string(((current_block_neighbours.top()).first).second) + ",";
                            outstring += to_string(assembly_count - (current_block_neighbours.top()).second + 1) + ",";
                        }
                    current_block_neighbours.pop();
                }
                // cout<<"\n";
                rand_select_start = false;
                current_max = 0;
            }
        }
        if(current_block_neighbours.size() < k_neighbours){
            current_block_neighbours.push( *it_pair_counts );
            if(it_pair_counts->second > current_max)
                current_max = it_pair_counts->second;
        }
        else{
            while((current_block_neighbours.top()).second < it_pair_counts->second && current_block_neighbours.size() > k_neighbours){
                // cout << "\t\tRemoved neighbour pair: "<<"(("<<((current_block_neighbours.top()).first).first<<",";
                // cout<<((current_block_neighbours.top()).first).second <<"),"<<(current_block_neighbours.top()).second<<")\n";
                current_block_neighbours.pop();
            }
            if((current_block_neighbours.top()).second < it_pair_counts->second){
                // cout << "\t\tRemoved neighbour pair: "<<"(("<<((current_block_neighbours.top()).first).first<<",";
                // cout<<((current_block_neighbours.top()).first).second <<"),"<<(current_block_neighbours.top()).second<<")\n";
                current_block_neighbours.pop();
                current_block_neighbours.push( *it_pair_counts );
                if(it_pair_counts->second > current_max)
                    current_max = it_pair_counts->second;
            }
            else if ((current_block_neighbours.top()).second == it_pair_counts->second){
                current_block_neighbours.push( *it_pair_counts );
                if(it_pair_counts->second > current_max)
                    current_max = it_pair_counts->second;
            }
        }
    }

    if(outstring != ""){
        outstring.pop_back();
        outstring += "\n";
    }
    outstring += to_string(((current_block_neighbours.top()).first).first) + ",";
    // cout<<((current_block_neighbours.top()).first).first<<"\n";

    // current_min_required = 0.95*current_max;
    if(current_max>min_presence_count)
        current_min_required = (0.95*current_max > current_max-10)?0.95*current_max:current_max-10;
    else
        current_min_required = min_presence_count;

    // cout<<((current_block_neighbours.top()).first).first<<" "<<current_max<<" "<<current_min_required<<"\n";

    while(!current_block_neighbours.empty() && (current_block_neighbours.top()).second < current_min_required)
        current_block_neighbours.pop();
    
    // Select random sample if more than k neighbours detected
    if(current_block_neighbours.size() > k_neighbours){
        rand_select_start = true;
        excess = current_block_neighbours.size() - k_neighbours + 1;
        rand_start = rand() % excess;
        pos = 0;
    }

    // // Select random sample if more than k neighbours detected
    // if(current_block_neighbours.size() > k_neighbours){
    //     rand_select_start = true;
    //     excess = current_block_neighbours.size() - k_neighbours + 1;
    //     rand_start = rand() % excess;
    //     pos = 0;
    // }
    // if(outstring != ""){
    //     outstring.pop_back();
    //     outstring += "\n";
    // }
    // outstring += to_string(((current_block_neighbours.top()).first).first) + ",";

    while(!current_block_neighbours.empty()){
        // cout<<"(("<<((current_block_neighbours.top()).first).first<<","<<((current_block_neighbours.top()).first).second;
        // cout<<"),"<<(current_block_neighbours.top()).second<<")";
        // current_block_neighbours.pop();
        // if(rand_select_start)
        //     if(rand_start == pos || pos >= excess)
        //         cout<<"**";
        // cout<<"\t";
        if(rand_select_start){
            if(rand_start == pos){
                if((current_block_neighbours.top()).second >= min_presence_count){
                    outstring += to_string(((current_block_neighbours.top()).first).second) + ",";
                    outstring += to_string(assembly_count - (current_block_neighbours.top()).second + 1) + ",";
                }
            }
            pos++;
            if(pos >= excess)
                rand_select_start = false;
        }
        else
            if((current_block_neighbours.top()).second >= min_presence_count){
                outstring += to_string(((current_block_neighbours.top()).first).second) + ",";
                outstring += to_string(assembly_count - (current_block_neighbours.top()).second + 1) + ",";
            }
        current_block_neighbours.pop();
    }
    if(outstring != ""){
        outstring.pop_back();
        outstring += "\n";
    }
    // cout<<"\n";

    // Save the selected k-nearest-neighbour subgraph edges
    // Format: <start_block_idx>,<neighbour_1_block_idx>,...,<neighbour_k_block_idx>
    outFile.open(neighbour_subgraph_outdir + to_string(partition_idx), ios::out);
    outFile << outstring;
    outFile.close();

    cout<<"generate_k_neighbour_subgraph ended: "<<neighbour_subgraph_outdir<<" "<<partition_idx<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_k_neigbour_subgraph_generator.cpp "<<argc<<"\n";
    short k_neighbours = static_cast<short>(stoi(argv[1]));
    unsigned long assembly_count = stoul(argv[2]);
    short partition_idx = static_cast<short>(stoi(argv[3]));
    unsigned long min_presence_count = stoul(argv[4]);
    string neighbour_pair_partition_dir = argv[5];
    string neighbour_subgraph_outdir = argv[6];

    generate_k_neighbour_subgraph(neighbour_pair_partition_dir, assembly_count, partition_idx, k_neighbours, min_presence_count,
                                neighbour_subgraph_outdir);

    return 0;
}

