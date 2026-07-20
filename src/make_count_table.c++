// Copyright 2020 Nadia Davidson. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** 
 ** Take a table of trsnscript counts and convert to gene-level counts.
 ** Assumes reads are uniquely mapped to transcripts.
 **
 ** Author: Nadia Davidson
 ** Modified: 2020
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <regex>
#include <stdlib.h>
#include <algorithm>    

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << endl;
  cerr << "Usage: <read/transcripts table> > make_count_table <ref_table> <id regular expression> > <out table>" << endl;
  cerr << endl;
}

// Main does the I/O and call multi_gene for each read and its associated
// set of transcript alignments.
int main(int argc, char **argv){

  //wrong number of arguements. Print help.
  if(argc!=3){
    print_usage();
    exit(1);
  }

  //Get the gene names and positions.
  //This is taken from process_transcriptome_alignment_table
  //Should probably separate this out as a function to remove redundancy
  /** 
   ** Now read in the gene to transcript ID mapping 
   **/
  ifstream file; 
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  //assume the first line is the header
  //use this to work out which columns to get
  map<string,string > trans_gene_map;
  string line;
  getline (file,line);
  int n_name=0; int n_name2=0; //column index
  istringstream line_stream(line);
  for(int i=0; line_stream ; i++){
    string temp;
    line_stream >> temp;
    if(temp=="name") n_name=i;
    if(temp=="name2") n_name2=i;
  }
  //loop over transcripts and fill information into data structure.
  while ( getline (file,line) ){ 
    istringstream line_stream(line);
    string gene; string trans; 
    for(int i=0; line_stream ; i++){
      string temp;
      line_stream >> temp;
      if(i==n_name) trans=temp; 
      if(i==n_name2) gene=temp; 
    }
    trans_gene_map[trans]=gene;
  }
  file.close();
  cerr << "Done reading in transcript IDs." << endl;

  //Now read the countTable. Store everything as a map
  map<string,int> gene_counts_map;
  
  /**********  Read all the alignments ****************/
  cerr << "Reading the input transcript alignments." << endl;
  string current_read="";
  //  string line;
  set<string> matched_genes;

  // anno_prefix comes from argv and is constant for the whole run, so compile
  // the regex once here rather than reconstructing it on every read (regex
  // construction is far more expensive than a match).
  string anno_reg=argv[2];
  bool have_anno_re = !anno_reg.empty();
  regex anno_re;
  if(have_anno_re) anno_re = regex(anno_reg);

  while(getline(std::cin,line) ){
    istringstream line_stream(line);
    string read, transcript;
    line_stream >> read >> transcript;

    // If we've moved to a new read, flush counts for the previous one
    if(read != current_read && current_read!=""){
        for(const auto& gene : matched_genes){
            gene_counts_map[gene]++;
        }
        matched_genes.clear();
    }

    current_read = read;

    // Resolve `transcript` (the raw SAM subject ID) to a real transcript
    // ID that exists as a key in trans_gene_map, using the same validated
    // fallback chain as extract_trans_id() / make_simple_read_table.c++.
    // The old code did `trans_gene_map[m[1].str()]` unconditionally --
    // std::map::operator[] default-constructs an empty string for any
    // missing key, so an unresolved (or wrongly-extracted, e.g. via a
    // permissive anno_prefix like "(.*)") ID silently produced gene=""
    // instead of skipping the read, corrupting every count.
    string trans_id = transcript; // default: Strategy 3, unresolved passthrough
    bool resolved = false;

    // Strategy 1: configured anno_prefix regex, but only trust the match
    // if it resolves to a real transcript ID. A greedy pattern like "(.*)"
    // will always "match" without producing a usable ID.
    if(have_anno_re){
      smatch m;
      if(regex_search(transcript, m, anno_re) && m.size() > 1){
        string candidate = m[1].str();
        if(trans_gene_map.count(candidate)){
          trans_id = candidate;
          resolved = true;
        }
      }
    }

    // Strategy 2: legacy UCSC-style headers with embedded __range= coords.
    if(!resolved){
      size_t range_pos = transcript.find("__range=");
      if(range_pos != string::npos){
        string before_range = transcript.substr(0, range_pos);
        size_t last_underscore = before_range.rfind('_');
        if(last_underscore != string::npos){
          string candidate = before_range.substr(last_underscore+1);
          if(trans_gene_map.count(candidate)){
            trans_id = candidate;
            resolved = true;
          }
        }
        if(!resolved && trans_gene_map.count(before_range)){
          trans_id = before_range;
          resolved = true;
        }
      }
    }
    // Strategy 3 (default set above): plain transcript ID, no embedded
    // coordinates -- used as-is when Strategies 1 and 2 don't resolve.

    // Only count reads that actually resolve to a known transcript --
    // never insert an empty-string "gene" for unresolved IDs.
    if(trans_gene_map.count(trans_id)){
      matched_genes.insert(trans_gene_map[trans_id]);
    }
  }
  // Flush the final read
  for (const auto& gene : matched_genes) {
    gene_counts_map[gene]++;
  }
  file.close();

  //loop over count map and output
  map<string, int>::iterator it = gene_counts_map.begin();
  while(it!=gene_counts_map.end()){
    cout << it->first << "\t" << it->second << endl;
    it++;
  }

  cerr << "Finished." << endl;

  return(0);
}
