#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <fstream>


using namespace std;


typedef unsigned long long int ulli;


string get_nucleotide_string(ulli value, short length){
    string result = string(length, 'A');
    int pos = length-1;
    while(value>0){
        char nuc;
        switch(value%4){
            case 0: nuc='A'; break;
            case 1: nuc='C'; break;
            case 2: nuc='G'; break;
            case 3: nuc='T'; break;
        }
        result[pos--] = nuc;
        value /= 4;
    }
    return result;
}


string get_kmer_string(int prefix_idx, ulli suffix_idx, short kmer_len=25, short prefix_len=5){
    return get_nucleotide_string(prefix_idx, prefix_len) + get_nucleotide_string(suffix_idx, kmer_len-prefix_len);
}


string get_reverse_kmer(string kmer){
    string rev_kmer = string(kmer.length(), 'A');
    unsigned short pos = kmer.length()-1;
    for(string::iterator it=kmer.begin(); it!=kmer.end(); it++){
        switch(*it){
            case 'A': rev_kmer[pos--] = 'T';
            break;
            case 'C': rev_kmer[pos--] = 'G';
            break;
            case 'G': rev_kmer[pos--] = 'C';
            break;
            case 'T': rev_kmer[pos--] = 'A';
            break;
        }
    }
    return rev_kmer;
}


void generate_kmer_fasta_from_idx(string assembly_name, string retained_binned_kmer_dir, string outdir, short kmer_len=25,
                                short prefix_count=1024, short prefix_len=5, char kmer_separator='Z'){
    ofstream outFile;
    string outfilename = outdir + assembly_name + ".fasta";
    outFile.open(outfilename, ios_base::app);
    outFile<<">"<<assembly_name<<"\n";

    ifstream inFile;

    string kmer;
    ulli total_kmer_count = 0;

    for(short prefix_idx=0; prefix_idx<prefix_count; prefix_idx++){
        string retained_kmer_filename = retained_binned_kmer_dir + assembly_name + "_" + to_string(prefix_idx);
        inFile.open(retained_kmer_filename,ios::binary);

        ulli retained_count;
        inFile.read((char*) (&retained_count), sizeof(retained_count));
        total_kmer_count += retained_count;

        for(ulli idx=0; idx<retained_count; idx++){
            ulli suffix_idx;
            inFile.read((char*) (&suffix_idx), sizeof(suffix_idx));
            kmer = get_kmer_string(prefix_idx, suffix_idx, kmer_len, prefix_len);
            outFile<<kmer<<kmer_separator;
        }

        inFile.close();
    }
    outFile<<get_reverse_kmer(kmer)<<kmer_separator;
    cout<<"generate_kmer_fasta_from_idx ended: "<<assembly_name<<" "<<total_kmer_count<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_fasta_generator.cpp "<<argv<<"\n";
    string assembly_name = argv[1];
    string retained_binned_kmer_dir = argv[2];
    string fasta_outdir = argv[3];
    short kmer_len = static_cast<short>(stoi(argv[4]));

    generate_kmer_fasta_from_idx(assembly_name, retained_binned_kmer_dir, fasta_outdir, kmer_len);

    return 0;
}