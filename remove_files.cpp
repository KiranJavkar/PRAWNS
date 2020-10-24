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

int main(int argc, char** argv){
    cout<<"remove_files.cpp "<<argc<<"\n";
    string listing_filename = argv[1];
    short remove_start = static_cast<short>(stoi(argv[2]));
    short remove_offset = static_cast<short>(stoi(argv[3]));

	unsigned long current_idx=0, to_remove_idx = remove_start;
	string text;
	ifstream input_file(listing_filename);
    while (getline( input_file, text )){
    	if(current_idx++ == to_remove_idx){
    		remove(text.c_str());
    		to_remove_idx += remove_offset;
    	}
    }

    return 0;
}