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


list< pulliulli > get_suffix_tuple_list(string filename){
    ifstream inFile;
    inFile.open(filename,ios::binary);//|ios::in);
    ulli binned_count;
    inFile.read((char*) (&binned_count), sizeof(binned_count));
    list< pulliulli > suffix_tuple_list;
    // cout<<filename<<"\t"<<binned_count<<"\n";
    for (ulli idx=0; idx<binned_count; idx++){
        ulli suffix_idx;
        inFile.read((char*) (&suffix_idx), sizeof(suffix_idx));
        pulliulli current_tuple = make_pair(suffix_idx, 1);
        suffix_tuple_list.push_back(current_tuple);
    }
    inFile.close();
    return suffix_tuple_list;
}


list< pullipullib > get_suffix_position_paired_list(string filename_prefix){
    ifstream inFile_suffix, inFile_pos;
    inFile_suffix.open(filename_prefix+"mer",ios::binary);//|ios::in);
    inFile_pos.open(filename_prefix+"pos",ios::binary);//|ios::in);

    ulli binned_count, check, suffix_idx;
    inFile_suffix.read((char*) (&binned_count), sizeof(binned_count));
    inFile_pos.read((char*) (&check), sizeof(check));
    assert(binned_count==check);

    list< pullipullib > suffix_position_paired_list;
    // cout<<filename_prefix<<"\t"<<binned_count<<"\n";

    pullib pos_strand_pair;

    for (ulli idx=0; idx<binned_count; idx++){
        inFile_suffix.read((char*) (&suffix_idx), sizeof(suffix_idx));
        inFile_pos.read((char*) (&pos_strand_pair), sizeof(pos_strand_pair));
        suffix_position_paired_list.push_back( pullipullib(suffix_idx, pos_strand_pair) );
    }

    inFile_suffix.close();
    inFile_pos.close();
    return suffix_position_paired_list;
}


list<pulliulli> merge_tuple_lists(list<pulliulli>& tuple_list_1, list<pulliulli>& tuple_list_2){
    list< pulliulli > merged_tuple_list;
    list< pulliulli >::iterator list_it_1 = tuple_list_1.begin();
    list< pulliulli >::iterator list_it_2 = tuple_list_2.begin();

    while(list_it_1 != tuple_list_1.end() && list_it_2 != tuple_list_2.end()){
        if(list_it_1->first == list_it_2->first){
            pulliulli merged_tuple = make_pair(list_it_1->first, (list_it_1->second) + (list_it_2->second));
            merged_tuple_list.push_back(merged_tuple);
            list_it_1++;
            list_it_2++;
        }
        else if(list_it_1->first < list_it_2->first){
            merged_tuple_list.push_back(*list_it_1);
            list_it_1++;
        }
        else{
            merged_tuple_list.push_back(*list_it_2);
            list_it_2++;
        }
    }
    while(list_it_1 != tuple_list_1.end()){
        merged_tuple_list.push_back(*list_it_1);
        list_it_1++;
    }
    while(list_it_2 != tuple_list_2.end()){
        merged_tuple_list.push_back(*list_it_2);
        list_it_2++;
    }
    return merged_tuple_list;
}


list< pulliulli > merge_suffix_tuple_list_wrapper(queue< list< pulliulli > >& list_queue){
    list< pulliulli > first_tuple_list;
    while(!list_queue.empty()){
        first_tuple_list = list_queue.front();
        list_queue.pop();
        if(list_queue.empty()){
            return first_tuple_list;
        }
        list< pulliulli > second_tuple_list = list_queue.front();
        list_queue.pop();
        list_queue.push(merge_tuple_lists(first_tuple_list, second_tuple_list));
    }
    return first_tuple_list;
}


list<string> get_lines_list(const char* FILENAME){
    ifstream input_file(FILENAME);
    string line;
    list<string> lines;
    while (getline( input_file, line )){
        lines.push_back(line);
    }
    return lines;
}


list<ulli> filter_tuples_for_min_presence(const list <pulliulli>& tuple_list, unsigned long threshold=1){
    list<ulli> filtered_list;
    for(list< pulliulli >::const_iterator suf_tup_it = tuple_list.begin(); suf_tup_it != tuple_list.end(); suf_tup_it++){
        if(suf_tup_it->second >= threshold)
            filtered_list.push_back(suf_tup_it->first);
    }
    return filtered_list;
}


list<pullipullib> get_required_suffix_position_paired_list(list<ulli>& filtered_list, list< pullipullib >& suffix_position_paired_list){
    list< pullipullib > required_suffix_position_paired_list;
    list< ulli >::iterator list_it_1 = filtered_list.begin();
    list< pullipullib >::iterator list_it_2 = suffix_position_paired_list.begin();

    while(list_it_1!=filtered_list.end() && list_it_2!=suffix_position_paired_list.end()){

        if(*list_it_1 == list_it_2->first){
            required_suffix_position_paired_list.push_back(*list_it_2);
            list_it_1++;
            list_it_2++;
        }
        else if(*list_it_1 < list_it_2->first)
            list_it_1++;
        else
            list_it_2++;
    }
    return required_suffix_position_paired_list;
}


/*void save_assemblywise_filtered_kmers(list<ulli>& filtered_list, unsigned long assembly_count, string prefix_idx_str,
                                        string binned_kmer_dir, string filtered_kmer_dir){

    cout<<"save_assemblywise_filtered_kmers started: "<<prefix_idx_str<<"\n";

    for(unsigned long assembly_idx=0; assembly_idx<assembly_count; assembly_idx++){
        string filename_prefix = binned_kmer_dir + to_string(assembly_idx) + "_" + prefix_idx_str + "_";
        list< pullipullib > suffix_position_paired_list = get_suffix_position_paired_list(filename_prefix);
    
        list< pullipullib > required_suffix_position_paired_list = get_required_suffix_position_paired_list(filtered_list, suffix_position_paired_list);
    
        // Remove binned kmer files
        string remove_file_cmd = "rm -rf " + filename_prefix + "*";
        system(remove_file_cmd.c_str());
    
        ofstream outFile;
        string retained_kmer_filename = filtered_kmer_dir + to_string(assembly_idx) + "_" + prefix_idx_str;
        ulli retained_count = required_suffix_position_paired_list.size();
        cout<<"Retained kmers "<<retained_kmer_filename<<" "<<retained_count<<"\n";
        outFile.open(retained_kmer_filename,ios::binary);//|ios::out);
        outFile.write((char*) (&retained_count), sizeof(retained_count));

        for (list<pullipullib>::iterator it = required_suffix_position_paired_list.begin(); it != required_suffix_position_paired_list.end(); it++){
            pullipullib& suffix_pos_strand_pair = *it;
            outFile.write((char*) (&suffix_pos_strand_pair), sizeof(suffix_pos_strand_pair));
        }
        outFile.close();
    }
}*/


int main(int argc, char** argv){
    // ifstream inFile;
    cout<<"kmer_filtering.cpp "<<argc<<"\n";
    queue < list<pulliulli> > list_queue;
    queue < list<pulliulli> > intermediate_list_queue;
    short intermediate_merge_count = 32;
    // list<string> assemblies = get_assembly_names(argv[1]);
    unsigned long assembly_count = stoul(argv[1]);
    string prefix_idx_str = argv[2];
    string binned_kmer_dir = argv[3];
    if(binned_kmer_dir.length()>1 && binned_kmer_dir[binned_kmer_dir.length()-1]!='/')
        binned_kmer_dir += "/";
    unsigned long min_presence_count = stoul(argv[4]);
    string filtered_kmer_dir = argv[5];
    if(filtered_kmer_dir.length()>1 && filtered_kmer_dir[filtered_kmer_dir.length()-1]!='/')
        filtered_kmer_dir += "/";


    for(unsigned long assembly_idx=0; assembly_idx<assembly_count; assembly_idx++){
        string filename = binned_kmer_dir + to_string(assembly_idx) + "_" + prefix_idx_str + "_mer";
        list< pulliulli > suffix_tuple_list = get_suffix_tuple_list(filename);

        if(suffix_tuple_list.size() > 0)
            list_queue.push(suffix_tuple_list);
        if(list_queue.size() == intermediate_merge_count){
            // cout<<"Pre-merge queue size: "<<list_queue.size()<<"\n";
            intermediate_list_queue.push(merge_suffix_tuple_list_wrapper(list_queue));
            // cout<<"Post-merge queue size: "<<list_queue.size()<<"\n";
        }
    }
    // cout<<"Tuple lists added to queue\n";

    if(list_queue.size() > 0)
        intermediate_list_queue.push(merge_suffix_tuple_list_wrapper(list_queue));

    list< pulliulli > merged_tuple_list = merge_suffix_tuple_list_wrapper(intermediate_list_queue);
    // cout<<"Tuple lists merged\t"<<merged_tuple_list.size()<<"\n";


    list< ulli > filtered_list = filter_tuples_for_min_presence(merged_tuple_list, min_presence_count);
    // cout<<"Merged tuple list filtered\t"<<filtered_list.size()<<"\n";


    // save_assemblywise_filtered_kmers(filtered_list, assembly_count, prefix_idx_str, binned_kmer_dir, filtered_kmer_dir);


    ofstream outFile;
    string retained_kmer_filename = filtered_kmer_dir + prefix_idx_str;
    ulli filtered_kmer_count = filtered_list.size();
    cout<<retained_kmer_filename<<" "<<filtered_kmer_count<<"\n";
    outFile.open(retained_kmer_filename, ios::binary);
    outFile.write((char*) (&filtered_kmer_count), sizeof(filtered_kmer_count));
    for (list< ulli >::iterator it = filtered_list.begin(); it != filtered_list.end(); it++){
        ulli& suffix_idx = *it;
        outFile.write((char*) (&suffix_idx), sizeof(suffix_idx));
    }
    outFile.close();

    cout<<"kmer_filtering.cpp completed: "<<prefix_idx_str<<"\n";

    return 0;
}