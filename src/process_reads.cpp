#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "seqreader.hpp"
#include "readcounts.hpp"
#include "taxdb.hpp"
#include "classify_io.cpp"
#include <sstream>


tuple<uint64_t, uint64_t, uint64_t> process_file(
  char *filename, unordered_map<uint32_t, READCOUNTS> taxon_counts,
  bool Print_kraken, bool Print_classified, bool Print_unclassified, bool Print_Progress) {

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
        tuple<unordered_map<uint32_t, uint32_t>, vector<uint32_t>, vector<uint8_t>>
        read_hit_counts_tuple = get_hit_count_map(cur_bc[i], my_taxon_counts);
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
        uint32_t call = classify_hit_count_map(
          cur_bc[i],
          my_taxon_counts,
          all_read_hit_counts[i],
          pruned_hit_counts
        );

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

      #pragma omp critical(write_output)
      {
        print_progress(
          my_total_classified, my_taxon_counts, taxon_counts,
          total_classified, total_sequences, total_bases,
          Print_kraken, Print_classified, Print_unclassified, Print_Progress
        );
      }
    }
  }  // end parallel section

  delete reader;
  return make_tuple(total_classified, total_sequences, total_bases)
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