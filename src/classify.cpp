/*
 * Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"
#include "readcounts.hpp"
#include "taxdb.hpp"
#include "gzstream.h"
#include "uid_mapping.hpp"
#include <sstream>

const size_t DEF_WORK_UNIT_SIZE = 500000;
int New_taxid_start = 1000000000;

using namespace std;
using namespace kraken;

#define USE_KHSET_FOR_EXACT_COUNTING

#ifdef EXACT_COUNTING
  #ifdef USE_KHSET_FOR_EXACT_COUNTING
    #include "khset.h"
    using READCOUNTS = ReadCounts< khset64_t >;
  #else
    #include <unordered_set>
    using READCOUNTS = ReadCounts< unordered_set<uint64_t> >;
  #endif
#else
  using READCOUNTS = ReadCounts<HyperLogLogPlusMinus<uint64_t> >;
#endif



void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename);
inline void print_sequence(ostringstream* oss_ptr, const DNASequence& dna);
string hitlist_string(const vector<uint32_t> &taxa, const vector<char>& ambig_list);


set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);
double get_seconds(struct timeval time1, struct timeval time2);

tuple<
  unordered_map<uint32_t, uint32_t>,vector<uint32_t>, vector<uint8_t>
> get_hit_count_map(DNASequence &dna,
		    unordered_map<uint32_t, READCOUNTS>& my_taxon_counts);
uint32_t classify_hit_count_map(DNASequence &dna,
                       unordered_map<uint32_t, READCOUNTS>& my_taxon_counts,
                       unordered_map<uint32_t, uint32_t> read_hit_counts,
                       unordered_map<uint32_t, uint32_t> bc_hit_counts);
bool handle_call(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       vector<uint32_t> taxa, vector<uint8_t> ambig_list, uint32_t call);
unordered_map<uint32_t, READCOUNTS> taxon_counts; // stats per taxon

int Num_threads = 1;
vector<string> DB_filenames;
vector<string> Index_filenames;
bool Quick_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Print_kraken_report = false;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
bool Print_sequence = false;
bool Print_Progress = true;
bool full_report = false;

bool Map_UIDs = false;
string UID_to_TaxID_map_filename;
map<uint32_t, vector<uint32_t> > UID_to_taxids_map;
QuickFile UID_to_TaxID_map_file;

uint32_t Minimum_hit_count = 1;
unordered_map<uint32_t, uint32_t> Parent_map;
unordered_map<uint32_t, vector<uint32_t> > Uid_dict;
string Classified_output_file, Unclassified_output_file, Kraken_output_file, Report_output_file, TaxDB_file;

ostream *Kraken_output;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Report_output;
vector<ofstream*> Open_fstreams;
vector<ogzstream*> Open_gzstreams;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;
TaxonomyDB<uint32_t> taxdb;
static vector<KrakenDB*> KrakenDatabases (DB_filenames.size());

struct db_status {
  db_status() : current_bin_key(0), current_min_pos(1), current_max_pos(0) {}
  uint64_t current_bin_key;
  int64_t current_min_pos;
  int64_t current_max_pos;
};

uint64_t total_classified = 0;
uint64_t total_sequences = 0;
uint64_t total_bases = 0;
uint32_t ambig_taxon = -1;

inline bool ends_with(std::string const & value, std::string const & ending)
{
        if (ending.size() > value.size()) return false;
            return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

ostream* cout_or_file(string file, bool append = false) {
    if (file == "-")
      return &cout;

    if (ends_with(file, ".gz")) {
      ogzstream* ogzs = new ogzstream(file.c_str());
      ogzs->exceptions( ifstream::failbit | ifstream::badbit );
      Open_gzstreams.push_back(ogzs);
      return ogzs;
    } else {
      ofstream* ofs = append? new ofstream(file.c_str(), std::ofstream::app) : new ofstream(file.c_str());
      ofs->exceptions( ifstream::failbit | ifstream::badbit );
      Open_fstreams.push_back(ofs);
      return ofs;
    }
}

void loadKrakenDB(KrakenDB& database, string DB_filename, string Index_filename) {
  QuickFile db_file;
  db_file.open_file(DB_filename);
  if (Populate_memory) {
    db_file.load_file();
  }
  database = KrakenDB(db_file.ptr());
  QuickFile idx_file;
  idx_file.open_file(Index_filename);
  if (Populate_memory)
    idx_file.load_file();

  KrakenDBIndex db_index(idx_file.ptr());
  database.set_index(&db_index);
}

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);
  
  if (Map_UIDs) {
    if (DB_filenames.size() > 1) {
      cerr << "Cannot use more than one database with UID mapping!" << endl;
      return 1;
    }

    cerr << "Reading UID mapping file " << UID_to_TaxID_map_filename << endl;
    UID_to_TaxID_map_file.open_file(UID_to_TaxID_map_filename);

    // Always Populate memory
    //if (Populate_memory) {
    UID_to_TaxID_map_file.load_file();
    //}
  }

  if (Populate_memory)
    cerr << "Loading database(s)... " << endl;

  static vector<QuickFile> idx_files (DB_filenames.size());
  static vector<QuickFile> db_files (DB_filenames.size());
  static vector<KrakenDBIndex> db_indices (DB_filenames.size());


  // TODO: Check DB_filenames and Index_filesnames have the same length
  for (size_t i=0; i < DB_filenames.size(); ++i) {
    cerr << " Database " << DB_filenames[i] << endl;
    db_files[i].open_file(DB_filenames[i]);
    if (Populate_memory)
      db_files[i].load_file();

    KrakenDatabases.push_back(new KrakenDB(db_files[i].ptr()));
    idx_files[i].open_file(Index_filenames[i]);
    if (Populate_memory)
      idx_files[i].load_file();
    db_indices[i] = KrakenDBIndex(idx_files[i].ptr());
    KrakenDatabases[i]->set_index(&db_indices[i]);
  }

  // TODO: Check all databases have the same k
  uint8_t kmer_size = KrakenDatabases[0]->get_k();
  for (size_t i = 1; i < KrakenDatabases.size(); ++i) {
    uint8_t kmer_size_i = KrakenDatabases[i]->get_k();
    if (kmer_size_i != kmer_size) {
      fprintf(stderr, "Different k-mer sizes in databases 1 and %lu: %i vs %i!\n", i+1, (int)kmer_size, (int)kmer_size_i);
      exit(1);
    }
  };
  KmerScanner::set_k(kmer_size);

  if (Populate_memory)
    cerr << "\ncomplete." << endl;


  if (!TaxDB_file.empty()) {
      taxdb = TaxonomyDB<uint32_t>(TaxDB_file, false);
      Parent_map = taxdb.getParentMap();
  } else {
      cerr << "TaxDB argument is required!" << endl;
      return 1;
  }

  if (Print_classified) {
    Classified_output = cout_or_file(Classified_output_file);
  }

  if (Print_unclassified) {
    Unclassified_output = cout_or_file(Unclassified_output_file);
  }

  if (! Kraken_output_file.empty()) {
    if (Kraken_output_file == "off" || Kraken_output_file == "-") {
      Print_kraken = false;
    //else if (Kraken_output_file == "-") {
    //  Kraken_output = &cout;
    } else {
      cerr << "Writing Kraken output to " << Kraken_output_file << endl;
      Kraken_output = cout_or_file(Kraken_output_file);
    }
  } else {
    Kraken_output = &cout;
  }

  //cerr << "Print_kraken: " << Print_kraken << "; Print_kraken_report: " << Print_kraken_report << "; k: " << uint32_t(KrakenDatabases[0]->get_k()) << endl;

  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (int i = optind; i < argc; i++)
    process_file(argv[i]);
  gettimeofday(&tv2, NULL);

  report_stats(tv1, tv2);

  if (!Report_output_file.empty() && Report_output_file != "off") {
    gettimeofday(&tv1, NULL);
    std::cerr << "Writing report file to " << Report_output_file <<"  ..\n";
    for (size_t i = 0; i < DB_filenames.size(); ++i) {
      const auto fname = DB_filenames[i] + ".counts";
      ifstream ifs(fname);
      bool counts_file_gd = false;
      if (ifs.good()) {
        if (ifs.peek() == std::ifstream::traits_type::eof()) {
          cerr << "Kmer counts file is empty - trying to regenerate ..." << endl;
        } else {
          ifs.close();
          counts_file_gd = true;
        }
      } 
      if (!counts_file_gd) {
        ofstream ofs(fname);
        cerr << "Writing kmer counts to " << fname << "... [only once for this database, may take a while] " << endl;
        auto counts = KrakenDatabases[i]->count_taxons();
        for (auto it = counts.begin(); it != counts.end(); ++it) {
          ofs << it->first << '\t' << it->second << '\n';
        }
        ofs.close();
      }
      taxdb.readGenomeSizes(fname);
    }
     Report_output = cout_or_file(Report_output_file, true);
  
    TaxReport<uint32_t,READCOUNTS> rep = TaxReport<uint32_t, READCOUNTS>(*Report_output, taxdb, taxon_counts, false);
    if (HLL_PRECISION > 0) {
      if (full_report) {
        rep.setReportCols(vector<string> { 
          "%",
          "reads", 
          "taxReads",
          "kmers",
          "taxKmers",
          "kmersDB",
          "taxKmersDB",
          "dup",
          "cov", 
          "taxID", 
          "rank", 
          "taxName"});
      } else {
        rep.setReportCols(vector<string> { 
          "%",
          "reads", 
          "taxReads",
          "kmers",
          "dup",
          "cov", 
          "taxID", 
          "rank", 
          "taxName"});
      }
    } else {
      rep.setReportCols(vector<string> { 
        "%",
        "reads", 
        "taxReads",
        "taxID", 
        "rank", 
        "taxName"});
    }
    rep.printReport("kraken");
    gettimeofday(&tv2, NULL);
    fprintf(stderr, "Report finished in %.3f seconds.\n", get_seconds(tv1,tv2));
  }
  cerr << "Finishing up ...";

  for (size_t i = 0; i < Open_fstreams.size(); ++i) {
    ofstream* ofs = Open_fstreams[i];
    ofs->close();
  }

  for (size_t i = 0; i < Open_gzstreams.size(); ++i) {
    ogzstream* ogzs = Open_gzstreams[i];
    ogzs->close();
  }

  return 0;
}

double get_seconds(struct timeval time1, struct timeval time2) {
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;
  return(seconds);
}

void report_stats(struct timeval time1, struct timeval time2) {
  double seconds = get_seconds(time1, time2);

  cerr << "\r";
  fprintf(stderr, 
          "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
          (unsigned long long) total_sequences, total_bases / 1.0e6, seconds,
          total_sequences / 1.0e3 / (seconds / 60),
          total_bases / 1.0e6 / (seconds / 60) );
  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) (total_sequences - total_classified),
          (total_sequences - total_classified) * 100.0 / total_sequences);
}

void process_file(char *filename) {

  uint64_t total_classified = 0;
  uint64_t total_sequences = 0;
  uint64_t total_bases = 0;
  string file_str(filename);
  BCReader *reader;
  reader = new BCReader(file_str);

  #pragma omp parallel
  {
    vector<DNASequence> cur_bc;
    ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;

    while (reader->is_valid()) { // Loop over barcodes
      cur_bc.clear();
      size_t total_nt = 0;

      cur_bc = reader->next_bc();
      
      unordered_map<uint32_t, READCOUNTS> my_taxon_counts;
      uint64_t my_total_classified = 0;
      kraken_output_ss.str("");
      classified_output_ss.str("");
      unclassified_output_ss.str("");

      unordered_map<uint32_t, uint32_t> bc_hit_counts;
      vector<unordered_map<uint32_t, uint32_t>> all_read_hit_counts;
      vector<vector<uint32_t>> all_read_taxa;
      vector<vector<uint8_t>> all_read_ambig;

      // Kmer Assignment Step: Loop over reads in the barcode keeping the following for each read
      // 
      // - map from taxa_id -> the number of kmers mapping to that taxa
      // - a vector of taxa for each read (0 == ambig)
      // - a vector with only 0/1 indicating if a kmer is (1 == ambiguous)
      // Also, count the number of kmers per taxa for all reads in the BC
      for (size_t i = 0; i < cur_bc.size(); i++) {
        tuple<
	  unordered_map<uint32_t, uint32_t>, vector<uint32_t>, vector<uint8_t>
	> read_hit_counts_tuple = get_hit_count_map(cur_bc[i], my_taxon_counts);
        all_read_taxa.push_back(get<1>(read_hit_counts_tuple));
        all_read_ambig.push_back(get<2>(read_hit_counts_tuple));

        unordered_map<uint32_t, uint32_t> read_hit_counts = get<0>(read_hit_counts_tuple);
        all_read_hit_counts.push_back(read_hit_counts);
        for (auto it = read_hit_counts.begin(); it != read_hit_counts.end(); ++it) {
          uint32_t taxon = it->first;
          bc_hit_counts[taxon]++;
        }
      }

      unordered_map<uint32_t, uint32_t> pruned_hit_counts;
      pruned_hit_counts = prune_tree(2, bc_hit_counts, Parent_map);


      // Read Assignment Step: Loop over reads in the barcode assigning each read to a taxa
      // 
      // Uses the maps/vectors stored from the kmer assignment step to
      // assign each read to a taxa. Use info from the whole BC to:
      // 1) improve the taxonomic rank of calls
      // 2) trim low abundance paths from the taxa tree
      for (size_t i = 0; i < cur_bc.size(); i++) {
        uint32_t call = classify_hit_count_map(cur_bc[i], my_taxon_counts, all_read_hit_counts[i], pruned_hit_counts);

        my_total_classified += handle_call(
          cur_bc[i],
          kraken_output_ss,
          classified_output_ss,
          unclassified_output_ss,
          all_read_taxa[i],
          all_read_ambig[i],
          call
        );
      }

      // Write output
      #pragma omp critical(write_output)
      {
        total_classified += my_total_classified;
        for (auto it = my_taxon_counts.begin(); it != my_taxon_counts.end(); ++it) {
          taxon_counts[it->first] += std::move(it->second);
        }

        if (Print_kraken)
          (*Kraken_output) << kraken_output_ss.str();
        if (Print_classified)
          (*Classified_output) << classified_output_ss.str();
        if (Print_unclassified)
          (*Unclassified_output) << unclassified_output_ss.str();
        
        total_sequences += cur_bc.size();
        total_bases += total_nt;
        if (Print_Progress) {  
          fprintf(stderr, 
            "\r Processed %lu sequences (%.2f%% classified)",
            total_sequences, total_classified * 100.0 / total_sequences
          );
        }
      }
    }
  }  // end parallel section

  delete reader;
}



inline void print_sequence(ostringstream* oss_ptr, const DNASequence& dna) {
      if (Fastq_input) {
        (*oss_ptr) << "@" << dna.header_line << endl
            << dna.seq << endl
            << "+" << endl
            << dna.quals << endl;
      }
      else {
        (*oss_ptr) << ">" << dna.header_line << endl
            << dna.seq << endl;
      }
}


string hitlist_string(const vector<uint32_t> &taxa, const vector<uint8_t> &ambig)
{
  int64_t last_code;
  int code_count = 1;
  ostringstream hitlist;

  if (ambig[0])   { last_code = -1; }
  else            { last_code = taxa[0]; }

  for (size_t i = 1; i < taxa.size(); i++) {
    int64_t code;
    if (ambig[i]) { code = -1; }
    else          { code = taxa[i]; }

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code >= 0) {
    hitlist << last_code << ":" << code_count;
  }
  else {
    hitlist << "A:" << code_count;
  }
  return hitlist.str();
}







set<uint32_t> get_ancestry(uint32_t taxon) {
  set<uint32_t> path;

  while (taxon > 0) {
    path.insert(taxon);
    taxon = Parent_map[taxon];
  }
  return path;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qfcC:U:Ma:r:sI:p:")) != -1) {
    switch (opt) {
      case 'd' :
        DB_filenames.push_back(optarg);
        break;
      case 'i' :
        Index_filenames.push_back(optarg);
        break;
      case 't' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive thread count");
        #ifdef _OPENMP
        if (sig > omp_get_num_procs())
          errx(EX_USAGE, "thread count exceeds number of processors");
        Num_threads = sig;
        omp_set_num_threads(Num_threads);
        #endif
        break;
      case 'p' :
        HLL_PRECISION = stoi(optarg);
        break;
      case 'q' :
        Quick_mode = true;
        break;
      case 'm' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive minimum hit count");
        Minimum_hit_count = sig;
        break;
      case 'f' :
        Fastq_input = true;
        break;
      case 'c' :
        Only_classified_kraken_output = true;
        break;
      case 'C' :
        Print_classified = true;
        Classified_output_file = optarg;
        break;
      case 'U' :
        Print_unclassified = true;
        Unclassified_output_file = optarg;
        break;
      case 'o' :
        Kraken_output_file = optarg;
        break;
      case 'r' :
        Report_output_file = optarg;
        break;
      case 's' :
        Print_sequence = true;
        break;
      case 'a' :
        TaxDB_file = optarg;
        break;
      case 'u' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive work unit size");
        Work_unit_size = sig;
        break;
      case 'M' :
        Populate_memory = true;
        break;
      case 'I' :
        UID_to_TaxID_map_filename = optarg;
        Map_UIDs = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filenames.empty()) {
    cerr << "Missing mandatory option -d" << endl;
    usage();
  }
  if (Index_filenames.empty()) {
    cerr << "Missing mandatory option -i" << endl;
    usage();
  }
  if (optind == argc && !Populate_memory) {
    cerr << "No sequence data files specified" << endl;
  }
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -r filename      Output file for Kraken report output" << endl
       << "  -a filename      TaxDB" << endl
       << "  -I filename      UID to TaxId map" << endl
       << "  -p #             Precision for unique k-mer counting, between 10 and 18" << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -f               Input is in FASTQ format" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -s               Print read sequence in Kraken output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "At least one FASTA or FASTQ file must be specified." << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}


tuple<unordered_map<uint32_t, uint32_t>, vector<uint32_t>, vector<uint8_t>>
get_hit_count_map(DNASequence &dna, unordered_map<uint32_t, READCOUNTS>& my_taxon_counts) {
  // Take a single sequence and return a triple
  // - map from taxa_id -> the number of kmers mapping to that taxa
  // - a vector of taxa for each read (0 == ambig)
  // - a vector with only 0/1 indicating if a kmer is (1 == ambiguous)
  // As a side effect keep a list of all 

  vector<uint32_t> taxa;
  vector<uint8_t> ambig_list;
  unordered_map<uint32_t, uint32_t> hit_counts;
  uint64_t *kmer_ptr;
  uint32_t taxon = 0;

  vector<db_status> db_statuses(KrakenDatabases.size());

  if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {  // check that the seq is longer than K
    size_t n_kmers = 1 + dna.seq.size() - KrakenDatabases[0]->get_k();
    taxa.reserve(n_kmers);
    ambig_list.reserve(n_kmers);
    KmerScanner scanner(dna.seq);

    while ((kmer_ptr = scanner.next_kmer()) != NULL) { // Iterate over kmers
      taxon = 0;
      if (scanner.ambig_kmer()) {
        ambig_list.push_back(1);
      } else {
        uint64_t cannonical_kmer = KrakenDatabases[0]->canonical_representation(*kmer_ptr);
        ambig_list.push_back(0);
        // go through multiple databases to map k-mer
        for (size_t i=0; i<KrakenDatabases.size(); ++i) {
          uint32_t* val_ptr = KrakenDatabases[i]->kmer_query(
            cannonical_kmer,
            &db_statuses[i].current_bin_key,
            &db_statuses[i].current_min_pos,
            &db_statuses[i].current_max_pos
          );
          if (val_ptr) {
            taxon = *val_ptr;
            break;
          }
        }

        my_taxon_counts[taxon].add_kmer(cannonical_kmer);

        if (taxon) {
          hit_counts[taxon]++;
        }
      }
      taxa.push_back(taxon);
    }
  }
  return make_tuple(hit_counts, taxa, ambig_list);
}


uint32_t classify_hit_count_map(DNASequence &dna, 
                       unordered_map<uint32_t, READCOUNTS>& my_taxon_counts,
                       unordered_map<uint32_t, uint32_t> read_hit_counts,
                       unordered_map<uint32_t, uint32_t> bc_hit_counts) {
  uint32_t call = 0;
  call = resolve_tree(read_hit_counts, Parent_map, bc_hit_counts);
  // TODO USE BC COUNTS TO PROMOTE CALL
  my_taxon_counts[call].incrementReadCount();
  return call;
}


bool handle_call(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       vector<uint32_t> taxa, vector<uint8_t> ambig_list, uint32_t call) {

  if (Print_unclassified && !call) 
    print_sequence(&uoss, dna);

  if (Print_classified && call)
    print_sequence(&coss, dna);


  if (! Print_kraken)
    return call;

  if (call) {
    koss << "C\t";
  }
  else {
    if (Only_classified_kraken_output)
      return false;
    koss << "U\t";
  }
  koss << dna.id << '\t' << call << '\t' << dna.seq.size() << '\t';

  if (taxa.empty())
    koss << "0:0";
  else
    koss << hitlist_string(taxa, ambig_list);

  if (Print_sequence)
      koss << "\t" << dna.seq;

  koss << "\n";
  return call;
}
