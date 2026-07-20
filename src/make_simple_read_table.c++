// Copyright 2019 Nadia Davidson. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** make_simple_read_table is a simple program to count the number of
 ** discordant pairs that support a fusions.
 **
 **
 ** Author: Nadia Davidson
 ** Modified: Feb 2026
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <regex>
#include <stdlib.h>

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << endl;
  cerr << "Usage: samtools view <bam file> | cut -f1-3 | make_simple_read_table <candidate_fusion_table> <trans_gene_map> <regex> >  <out table>" << endl;
  cerr << endl;
}

//struct for fusion information
struct fusion_candidate {
  string read ;
  int break_min ;
  int break_max ;
  pair<string,string> fusion ;
};

// the real stuff starts here.
int main(int argc, char **argv){

  if(argc!=4){
    print_usage();
    exit(1);
  }

  /**
   ** Read the list of candidate fusions
   **/
  ifstream file;
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  string line;
  vector<fusion_candidate>  candidate_table;
  vector< pair<string,string> > fusion_list;
  vector< string > genes_of_interest;
  // read in the table of candidate fusions
  const regex fusion_id_re("(.*):(.*)"); // constant pattern, compile once
  while ( getline (file,line) ){
    fusion_candidate new_cand;
    istringstream line_stream(line);
    line_stream >> new_cand.read;
    line_stream >> new_cand.break_min;
    line_stream >> new_cand.break_max;
    string temp;
    line_stream >> temp;
    smatch m; //extract the gene id
    regex_search(temp,m,fusion_id_re);
    //    if(m[1].str()<m[2].str())
      new_cand.fusion=make_pair(m[1].str(), m[2].str());
      //    else
      //     new_cand.fusion=make_pair(m[2].str(), m[1].str());
    //sort
    fusion_list.push_back(new_cand.fusion);
    candidate_table.push_back(new_cand);
    genes_of_interest.push_back(m[1].str());
    genes_of_interest.push_back(m[2].str());
  }
  file.close();
  // remove duplicates in the fusion list
  sort( fusion_list.begin(), fusion_list.end() );
  fusion_list.erase( unique( fusion_list.begin(), fusion_list.end() ), fusion_list.end() );
  sort( genes_of_interest.begin(), genes_of_interest.end() );
  genes_of_interest.erase( unique( genes_of_interest.begin(), genes_of_interest.end() ), genes_of_interest.end() );
  cerr << "Done reading in candidate fusions" << endl;

  /** 
   ** Now read in the gene to transcript ID mapping 
   **/
  file.open(argv[2]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[2] << endl;
    exit(1);
  }
  //assume the first line is the header
  getline (file,line);
  int n_name=0; int n_name2=0;
  istringstream line_stream(line);
  for(int i=0; line_stream ; i++){
    string temp;
    line_stream >> temp;
    if(temp=="name") n_name=i;
    if(temp=="name2") n_name2=i;
  }
  unordered_map<string, vector<string> > gene_trans_map ;
  while ( getline (file,line) ){
    istringstream line_stream(line);
    string gene;
    string trans;
    for(int i=0; line_stream ; i++){
      string temp;
      line_stream >> temp;
      if(i==n_name) trans=temp;
      if(i==n_name2) gene=temp;
    }
    //only include genes in candidate fusions
    if(find(genes_of_interest.begin(), genes_of_interest.end(), gene) != genes_of_interest.end())
      gene_trans_map[gene].push_back(trans);
  }
  file.close();
  cerr << "Done reading in transcript IDs" << endl;

  // Flat set of every transcript ID we actually care about (i.e. every
  // value in gene_trans_map), used below to validate extracted IDs before
  // trusting them -- the same principle as extract_trans_id() in
  // process_transcriptome_align_table.c++.
  unordered_set<string> valid_trans_ids;
  for(auto const& kv : gene_trans_map)
    for(auto const& t : kv.second)
      valid_trans_ids.insert(t);

  /** 
   ** Now read in the mapped reads
   **/
  unordered_map< string, vector < pair < string , bool > > > trans_read_map;
  vector<string> trans_names;
  int i=0;
  string sam_read;
  string sam_flag;
  string sam_trans;
  while(!cin.eof()){
    cin >> sam_read >> sam_flag >> sam_trans ;
    trans_read_map[sam_trans].push_back(make_pair(sam_read,sam_flag=="0"));
    if(i%1000000==0)
      cerr << i << " alignments read" << endl;
    i++;
  }
  cerr << "Done reading in bam file" << endl;
  // Loop over all the trans_read_map and fix the names:
  string anno_reg=argv[3];
  // anno_reg is constant for the whole run, so compile the regex once here
  // rather than reconstructing it for every entry in the loop below.
  bool have_anno_re = !anno_reg.empty();
  regex anno_re;
  if(have_anno_re) anno_re = regex(anno_reg);
  unordered_map< string, vector< pair < string, bool > > > trans_read_map_fixed;
  auto tr_itr=trans_read_map.begin();
  for(;tr_itr!=trans_read_map.end() ; tr_itr++){
    const string& subject_id = tr_itr->first;
    string key = subject_id; // default: Strategy 3, unresolved passthrough

    // Strategy 1: configured anno_prefix regex, but only trust the match
    // if it resolves to a real transcript ID. A permissive/greedy prefix
    // like "(.*)" will always "match" without producing a usable ID, so
    // unconditionally trusting it (the old behaviour) breaks any BAM
    // aligned against UCSC-style "__range=" headers.
    bool resolved = false;
    if(have_anno_re){
      smatch m;
      if(regex_search(subject_id, m, anno_re) && m.size() > 1){
        string candidate = m[1].str();
        if(valid_trans_ids.count(candidate)){
          key = candidate;
          resolved = true;
        }
      }
    }

    // Strategy 2: legacy UCSC-style headers with embedded __range= coords.
    if(!resolved){
      size_t range_pos = subject_id.find("__range=");
      if(range_pos != string::npos){
        string before_range = subject_id.substr(0, range_pos);
        size_t last_underscore = before_range.rfind('_');
        if(last_underscore != string::npos){
          string candidate = before_range.substr(last_underscore+1);
          if(valid_trans_ids.count(candidate)){
            key = candidate;
            resolved = true;
          }
        }
        if(!resolved && valid_trans_ids.count(before_range)){
          key = before_range;
          resolved = true;
        }
      }
    }
    // Strategy 3 (default set above): plain transcript ID, no embedded
    // coordinates -- used as-is when Strategies 1 and 2 don't resolve.

    trans_read_map_fixed[key]=tr_itr->second;
  }
  trans_read_map.clear();
  cerr << "Done getting trans ids" << endl;

  /**
   ** now loop over each gene and get all the reads aligning to it
   **/
  unordered_map<string, pair< unordered_set < string> , unordered_set< string > > > gene_reads_sense;
  unordered_map<string, pair< unordered_set < string> , unordered_set< string > > > gene_reads_antisense;
  auto gt_itr=gene_trans_map.begin();
  //loop over genes and get transcripts id
  for(;gt_itr!=gene_trans_map.end(); gt_itr++){
    vector<pair<string,bool>> reads;
    //for all transcripts look up the reads
    for(int t=0; t < gt_itr->second.size(); t++){
      vector<pair<string,bool> > reads_to_add = trans_read_map_fixed[gt_itr->second.at(t)];
      reads.insert(reads.end(), reads_to_add.begin(), reads_to_add.end());
    }
    // remove redundancy
    sort( reads.begin(), reads.end() );
    reads.erase( unique( reads.begin(), reads.end() ), reads.end() );   
    //now loop over the reads and separate the read start and ends
    for(int r=0; r<reads.size(); r++){
      //separate the read id and pair end
      string read_id=reads.at(r).first.substr(0,reads.at(r).first.find("/"));//reads.at(r).size()-1); 
      char read_end=reads.at(r).first.back(); //split by anti-sense to fusion and R1 or R2
      if( (read_end=='1') & reads.at(r).second ) 
	gene_reads_sense[gt_itr->first].first.insert(read_id);
      if( (read_end=='1') & !reads.at(r).second ) 
	gene_reads_antisense[gt_itr->first].first.insert(read_id);
      if( (read_end=='2') & !reads.at(r).second ) 
	gene_reads_sense[gt_itr->first].second.insert(read_id);
      if( (read_end=='2') & reads.at(r).second ) 
	gene_reads_antisense[gt_itr->first].second.insert(read_id); 

    }
  }
  trans_read_map_fixed.clear();
  gene_trans_map.clear();
  cerr << "Done assigning reads to genes" << endl;

  /**
   ** Calculate the counts for each fusion
   **/
  //loop over the fusion list
  map<pair<string,string >, int> spanning_reads;
  for(int f=0; f<fusion_list.size(); f++){
    //check from intersection of read ids.
    unordered_set<string> g1_r1=gene_reads_sense[fusion_list.at(f).first].first;
    unordered_set<string> g1_r2=gene_reads_antisense[fusion_list.at(f).first].second;
    unordered_set<string> g2_r1=gene_reads_antisense[fusion_list.at(f).second].first;
    unordered_set<string> g2_r2=gene_reads_sense[fusion_list.at(f).second].second;
    int total=0;
    unordered_set<string> temp_set1(g1_r1.begin(),g1_r1.end()); //sense supporting
    temp_set1.insert(g2_r2.begin(), g2_r2.end());
    total+=g1_r1.size() + g2_r2.size() - temp_set1.size();

    unordered_set<string> temp_set2(g2_r1.begin(),g2_r1.end()); //anti-sense supporting
    temp_set2.insert(g1_r2.begin(), g1_r2.end());
    total+=g1_r2.size() + g2_r1.size() - temp_set2.size();
    spanning_reads[fusion_list.at(f)]=total;
  }
  cerr << "Done calculating spanning pairs" << endl;

  cout << "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" << endl;
  vector<fusion_candidate>::iterator cand_itr=candidate_table.begin();
  for(;cand_itr!=candidate_table.end(); cand_itr++){
    cout << cand_itr->read << "\t" ;
    cout << cand_itr->break_min << "\t" ;
    cout << cand_itr->break_max << "\t" ;
    cout << cand_itr->fusion.first << ":" << cand_itr->fusion.second << "\t" ;
    cout << spanning_reads[cand_itr->fusion] << "\t" ;
    cout << "1" << endl;
  }

  return(0);
}
