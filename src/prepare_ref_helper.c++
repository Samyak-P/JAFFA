// prepare_ref_helper.c++
//
// C++ replacement for the Python glue logic previously embedded as
// heredocs inside prepare_jaffa_reference_minimal.sh. Provides the same
// four operations as standalone, compiled subcommands, consistent with
// the rest of JAFFA's tools/bin/ toolchain (no runtime Python dependency).
//
// Subcommands:
//   classify   <genome.fa>
//       Prints YES/NO -- whether sequence names contain unsafe characters
//       (|, :, /) or duplicates.
//
//   validate   <transcript.fa>
//       Checks transcript FASTA headers for duplicate IDs (hard error,
//       exit 1) and generic/auto-numbered-looking IDs (warning only).
//
//   sanitize   <genome.fa> <out.fa> <log.tsv>
//       Minimal-touch renaming: strips ENA pipe-format names to the last
//       field, replaces ':' and '/' with '_', de-duplicates collisions.
//       Writes a mapping log of only the names that actually changed.
//
//   rename-anno <anno.gtf> <log.tsv> <out.gtf>
//       Applies the sanitize log's mapping to an annotation's chromosome
//       column (column 1).
//
// Author: written for the JAFFA reference-prep pipeline, 2026.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <regex>
#include <cstring>

using namespace std;

static void print_usage(){
  cerr << "Usage:\n"
       << "  prepare_ref_helper classify <genome.fa>\n"
       << "  prepare_ref_helper validate <transcript.fa>\n"
       << "  prepare_ref_helper sanitize <genome.fa> <out.fa> <log.tsv>\n"
       << "  prepare_ref_helper rename-anno <anno.gtf> <log.tsv> <out.gtf>\n";
}

// Extract the first whitespace-delimited token from a FASTA header line,
// with the leading '>' already stripped by the caller.
static string first_token(const string& header_no_gt){
  size_t i = 0;
  while(i < header_no_gt.size() && !isspace((unsigned char)header_no_gt[i])) i++;
  return header_no_gt.substr(0, i);
}

static vector<string> read_fasta_names(const string& path){
  vector<string> names;
  ifstream in(path);
  if(!in){
    cerr << "Error: could not open " << path << endl;
    exit(1);
  }
  string line;
  while(getline(in, line)){
    if(!line.empty() && line[0] == '>'){
      names.push_back(first_token(line.substr(1)));
    }
  }
  return names;
}

static bool has_unsafe_char(const string& n){
  return n.find('|') != string::npos ||
         n.find(':') != string::npos ||
         n.find('/') != string::npos;
}

int cmd_classify(int argc, char** argv){
  if(argc != 3){ print_usage(); return 1; }
  vector<string> names = read_fasta_names(argv[2]);

  bool unsafe = false;
  for(const auto& n : names){
    if(has_unsafe_char(n)){ unsafe = true; break; }
  }

  unordered_set<string> seen;
  bool dup = false;
  for(const auto& n : names){
    if(!seen.insert(n).second){ dup = true; break; }
  }

  cout << ((unsafe || dup) ? "YES" : "NO") << endl;
  return 0;
}

int cmd_validate(int argc, char** argv){
  if(argc != 3){ print_usage(); return 1; }
  vector<string> names = read_fasta_names(argv[2]);

  unordered_map<string,int> counts;
  for(const auto& n : names) counts[n]++;

  vector<string> dupes;
  for(const auto& kv : counts) if(kv.second > 1) dupes.push_back(kv.first);

  // Case-insensitive match against ^(exon|CDS|mRNA|gene|transcript)[-_]?\d*$
  regex generic_pattern("^(exon|CDS|mRNA|gene|transcript)[-_]?[0-9]*$",
                         regex::icase);
  unordered_set<string> unique_names(names.begin(), names.end());
  vector<string> generic;
  for(const auto& n : unique_names){
    if(regex_match(n, generic_pattern)) generic.push_back(n);
  }

  if(!dupes.empty()){
    cerr << "ERROR: " << dupes.size()
         << " duplicate transcript IDs found in transcript FASTA." << endl;
    cerr << "Examples: ";
    for(size_t i=0; i<dupes.size() && i<5; ++i) cerr << dupes[i] << " ";
    cerr << endl;
    cerr << "This will cause makeblastdb to fail, or silently drop records in lookups." << endl;
    cerr << "Check your GTF/GFF3 for non-unique transcript_id / ID attributes" << endl;
    return 1;
  }

  if(!generic.empty()){
    cerr << "WARNING: " << generic.size()
         << " transcript IDs look generic/auto-numbered." << endl;
    cerr << "Examples: ";
    for(size_t i=0; i<generic.size() && i<5; ++i) cerr << generic[i] << " ";
    cerr << endl;
    cerr << "This usually means the source GFF3 lacks real unique transcript IDs." << endl;
    cerr << "Continuing, but check output carefully." << endl;
  }

  cerr << "Transcript ID validation passed: " << names.size()
       << " sequences, " << unique_names.size() << " unique." << endl;
  return 0;
}

static string sanitize_name(const string& name){
  string result = name;
  size_t pipe_pos = result.rfind('|');
  if(pipe_pos != string::npos){
    result = result.substr(pipe_pos + 1);
  }
  for(auto& c : result){
    if(c == ':' || c == '/') c = '_';
  }
  return result;
}

int cmd_sanitize(int argc, char** argv){
  if(argc != 5){ print_usage(); return 1; }
  string fasta_path = argv[2];
  string out_fa_path = argv[3];
  string log_path = argv[4];

  ifstream in(fasta_path);
  if(!in){ cerr << "Error: could not open " << fasta_path << endl; return 1; }
  ofstream out(out_fa_path);
  if(!out){ cerr << "Error: could not write " << out_fa_path << endl; return 1; }

  unordered_set<string> seen;
  vector<pair<string,string>> mapping; // original -> new, only where changed

  string line;
  while(getline(in, line)){
    if(!line.empty() && line[0] == '>'){
      string rest = line.substr(1);
      size_t sp = rest.find(' ');
      string name = (sp == string::npos) ? rest : rest.substr(0, sp);
      string desc = (sp == string::npos) ? "" : rest.substr(sp); // includes leading space

      string new_name = sanitize_name(name);
      string base = new_name;
      int n = 1;
      while(seen.count(new_name)){
        n++;
        new_name = base + "_dup" + to_string(n);
      }
      seen.insert(new_name);

      if(new_name != name){
        mapping.push_back({name, new_name});
      }
      out << ">" << new_name << desc << "\n";
    } else {
      out << line << "\n";
    }
  }

  ofstream log(log_path);
  if(!log){ cerr << "Error: could not write " << log_path << endl; return 1; }
  log << "#original\tsanitized\n";
  for(const auto& p : mapping){
    log << p.first << "\t" << p.second << "\n";
  }

  cerr << "Sanitized " << mapping.size()
       << " of the sequence names (rest left unchanged)" << endl;
  return 0;
}

int cmd_rename_anno(int argc, char** argv){
  if(argc != 5){ print_usage(); return 1; }
  string anno_path = argv[2];
  string log_path = argv[3];
  string out_anno_path = argv[4];

  unordered_map<string,string> mapping;
  {
    ifstream log(log_path);
    if(!log){ cerr << "Error: could not open " << log_path << endl; return 1; }
    string line;
    bool first = true;
    while(getline(log, line)){
      if(first){ first = false; continue; } // skip header
      if(line.empty()) continue;
      size_t tab = line.find('\t');
      if(tab == string::npos) continue;
      string orig = line.substr(0, tab);
      string neu = line.substr(tab + 1);
      // strip trailing \r if present (CRLF safety)
      if(!neu.empty() && neu.back() == '\r') neu.pop_back();
      mapping[orig] = neu;
    }
  }

  ifstream anno(anno_path);
  if(!anno){ cerr << "Error: could not open " << anno_path << endl; return 1; }
  ofstream out(out_anno_path);
  if(!out){ cerr << "Error: could not write " << out_anno_path << endl; return 1; }

  string line;
  while(getline(anno, line)){
    if(!line.empty() && line[0] == '#'){
      out << line << "\n";
      continue;
    }
    size_t tab = line.find('\t');
    if(tab == string::npos){
      out << line << "\n";
      continue;
    }
    string chrom = line.substr(0, tab);
    string rest = line.substr(tab);
    auto it = mapping.find(chrom);
    if(it != mapping.end()){
      out << it->second << rest << "\n";
    } else {
      out << line << "\n";
    }
  }
  return 0;
}

int main(int argc, char** argv){
  if(argc < 2){
    print_usage();
    return 1;
  }
  string subcmd = argv[1];

  if(subcmd == "classify")     return cmd_classify(argc, argv);
  if(subcmd == "validate")     return cmd_validate(argc, argv);
  if(subcmd == "sanitize")     return cmd_sanitize(argc, argv);
  if(subcmd == "rename-anno")  return cmd_rename_anno(argc, argv);

  cerr << "Unknown subcommand: " << subcmd << endl;
  print_usage();
  return 1;
}
