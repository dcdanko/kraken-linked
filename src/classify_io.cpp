

void print_progress(
    uint64_t my_total_classified, unordered_map<uint32_t, READCOUNTS> my_taxon_counts,
    unordered_map<uint32_t, READCOUNTS> taxon_counts,
    uint64_t total_classified, uint64_t total_sequences, uint64_t total_bases,
    bool Print_kraken, bool Print_classified, bool Print_unclassified, bool Print_Progress
    ){
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