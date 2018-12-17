/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
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
#include "seqreader.hpp"

using namespace std;

namespace kraken {

  BCReader::BCReader(string filename) {
    BCFastqReader::BCFastqReader(filename);
  }

  vector<DNASequence> BCReader::next_bc() {
    vector<DNASequence> bc;
    DNASequence next_seq;
    while (!cur_seq){
      cur_seq = bc_fastq_reader.next_sequence();
    }
    bc.push_back(cur_seq);
    next_seq = bc_fastq_reader.next_sequence();
    while(next_seq.bc == cur_seq.bc){
      bc.push_back(next_seq)
      next_seq = bc_fastq_reader.next_sequence();
    }
    cur_seq = next_seq;
    return bc;
  }

  bool BCReader::is_valid() {
    return BCFastqReader.is_valid()
  }

  BCFastqReader::BCFastqReader(string filename) {
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }

  DNASequence BCFastqReader::next_sequence() {
    DNASequence dna;

    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }

    string line;
    getline(file, line);
    if (line.empty()) {
      valid = false;  // Sometimes FASTQ files have empty last lines
      return dna;
    }
    if (line[0] != '@') {
      if (line[0] != '\r')
        warnx("malformed fastq file - sequence header (%s)", line.c_str());
      valid = false;
      return dna;
    }

    size_t bc_start = line.find("BX:");
    if (bc_start == -1){
      return dna;
    }
    size_t bc_end = line.find(' ', pos=bc_start);
    dna.bc = line.substr(bc_start, bc_end - bc_start);

    dna.header_line = line.substr(1);
    istringstream line_ss(dna.header_line);
    line_ss >> dna.id;


    getline(file, dna.seq);

    getline(file, line);
    if (line.empty() || line[0] != '+') {
      if (line[0] != '\r')
        warnx("malformed fastq file - quality header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    getline(file, dna.quals);

    return dna;
  }

  bool BCFastqReader::is_valid() {
    return valid;
  }


} // namespace
