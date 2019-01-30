/*
 * Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
 * Portions (c) 2019-2020, David Danko <dcd3001@med.cornell.edu> as part of KrakenLinked
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

#include "assert_helpers.h"
#include "kraken_headers.hpp"
#include "krakenutil.hpp"
#include <unordered_set>
#include <algorithm>

using namespace std;

namespace kraken {

  // Build a node->parent unordered_map from NCBI Taxonomy nodes.dmp file
  unordered_map<uint32_t, uint32_t> build_parent_map(string filename) {
    unordered_map<uint32_t, uint32_t> pmap;
    uint32_t node_id, parent_id;
    string line;
    ifstream ifs(filename.c_str());
    if (ifs.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "error opening %s", filename.c_str());
    }

    while (ifs.good()) {
      getline(ifs, line);
      if (line.empty())
        break;
      sscanf(line.c_str(), "%d\t|\t%d", &node_id, &parent_id);
      pmap[node_id] = parent_id;
    }
    pmap[1] = 0;
    return pmap;
  }

  // Return lowest common ancestor of a and b
  // LCA(0,x) = LCA(x,0) = x
  // Default ancestor is 1 (root of tree)
  uint32_t lca(const unordered_map<uint32_t, uint32_t> &parent_map,
    uint32_t a, uint32_t b)
  {
    if (a == 0 || b == 0)
      return a ? a : b;

    unordered_set<uint32_t> a_path;
    while (a > 1) {
      a_path.insert(a);
      auto a_it = parent_map.find(a);
      if (a_it == parent_map.end()) {
        cerr << "No parent for " << a << "!\n";
        break;
      } 
      a = a_it->second;
    }
    while (b > 1) {
      if (a_path.count(b) > 0)
        return b;

      auto b_it = parent_map.find(b);
      if (b_it == parent_map.end()) {
        cerr << "No parent for " << b << "!\n";
        break;
      } 
      b = b_it->second;
    }
    return 1;
  }

  uint32_t lca_vec(const unordered_map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b) {
    if (a == 0 || b == 0)
      return a ? a : b;

    // create a path from a to the root
    std::vector<uint32_t> a_path;
    do {
      if (a == b) 
        return a; 
      a_path.push_back(a);
      a = parent_map.at(a);
    } while (a != a_path.back());

    // search for b in the path from a to the root
    uint32_t last_b = 0;
    do {
      if (std::find(a_path.begin(), a_path.end(), b) != a_path.end())
        return b;

      last_b = b;
      b = parent_map.at(b);
    } while (last_b != b);
    return 1;
  }

  // Remove low abundance paths from the tree demoting their kmer
  // counts to their parents.
  unordered_map<uint32_t, uint32_t> prune_tree(
    uint32_t min_abundance,
    const unordered_map<uint32_t, uint32_t> &inp_hit_counts,
    const unordered_map<uint32_t, uint32_t> &parent_map) {
    uint32_t taxon, parent, abund;
    unordered_map<uint32_t, uint32_t> pruned_counts;
    unordered_map<uint32_t, uint32_t> hit_counts = inp_hit_counts;
    unordered_set<uint32_t> internal_nodes;
    bool modified;
    
    // Run loop until the pruned tree is the same as the unpruned tree
    modified = true;
    while(modified){
      modified = false;

      /*
       * Find the ancestors of all nodes in the unpruned tree and
       * keep track of them as internal (not leaf) nodes.
       *
       * Also add the hits from non-leaf nodes to the pruned tree
       * since we only remove leaf nodes each iteration. Note that
       * not all internal nodes have hits.
       */
      for(auto it=hit_counts.begin(); it!= hit_counts.end(); ++it){
        taxon = it->first;
      	auto parent_node = parent_map.find(taxon);
      	while(parent_node != parent_map.end()){
      	  internal_nodes.insert(parent_node->second);
      	  if(hit_counts.count(parent_node->second) > 0){
      	    pruned_counts[parent_node->second] = hit_counts.at(parent_node->second);
      	  }
      	  parent_node = parent_map.find(parent_node->second);
      	}
      }

      /*
       * Loop through leaf nodes in the unpruned tree. If they are
       * above threshold add the leaf node to the pruned tree. If
       * they are below threshold add the leaf node's abundance to
       * the parent (assume a count of zero if the parent is not
       * already present).
       */
      for(auto it=hit_counts.begin(); it!= hit_counts.end(); ++it){
        taxon = it->first;
        abund = it->second;
        if(internal_nodes.count(taxon) == 0){
	  auto parent_node = parent_map.find(taxon);
          if(abund >= min_abundance || parent_node == parent_map.end()){
            pruned_counts[taxon] = abund;
          } else {
	    parent = parent_node->second;
	    if(pruned_counts.count(parent) > 0){
	      pruned_counts[parent] += abund;
	    } else {
	      pruned_counts[parent] = abund;
	    }
	    modified = true;
          }
        }
      }

      // Prepare for the next iteration of the loop
      hit_counts = pruned_counts;
      pruned_counts = unordered_map<uint32_t, uint32_t>();
      internal_nodes.clear();
    }

    return hit_counts;
  }

  // Tree resolution: take all hit taxa (plus ancestors), then
  // return leaf of highest weighted leaf-to-root path.
  uint32_t resolve_tree(const unordered_map<uint32_t, uint32_t> &read_hit_counts,
                        const unordered_map<uint32_t, uint32_t> &parent_map,
                        const unordered_map<uint32_t, uint32_t> &bc_hit_counts) {
  
    set<uint32_t> max_taxa;
    uint32_t max_taxon = 0, max_score = 0;
    unordered_map<uint32_t, uint32_t> pruned_read_hit_counts;
    uint32_t taxon, parent;

    /*
     * Iterate through (taxon, kmer) hits in the read.
     * 
     * If taxon is in the pruned-barcode-tree add the(taxon, kmer)
     * to pruned-read-tree.
     * If taxon is not in the pruned-barcode-tree find the first
     * ancestor that is and add (ancestor, read-kmer + ancestor-kmer || 0)
     * to the pruned-read-tree.
     */
    for (auto it = read_hit_counts.begin(); it != read_hit_counts.end(); ++it) {
      taxon = it->first;
      if (bc_hit_counts.count(taxon) > 0) { // in pruned-bc-tree
        if(pruned_read_hit_counts.count(taxon) > 0){ // this node may be an ancestor of a pruned node 
          pruned_read_hit_counts[taxon] += it->second;
        } else {
          pruned_read_hit_counts[taxon] = it->second;
        }
      } else { // not in pruned-bc-tree
        parent = parent_map.at(taxon);
        while (bc_hit_counts.count(parent) == 0){ // find first ancestor in pruned-bc-tree
          parent = parent_map.at(parent);
        }
        if(pruned_read_hit_counts.count(parent) > 0){
          pruned_read_hit_counts[parent] += it->second;
        } else {
          pruned_read_hit_counts[parent] = it->second;
        }
      }
    }

    // Sum each taxon's LTR path in the pruned-read-tree
    for (auto it = pruned_read_hit_counts.begin(); it != pruned_read_hit_counts.end(); ++it) {
      uint32_t taxon = it->first;
      uint32_t node = taxon;
      uint32_t score = 0;
      while (node > 0) {
        auto it2 = pruned_read_hit_counts.find(node);
        if (it2 != pruned_read_hit_counts.end()) {
          score += it2->second;
        }
        auto node_it = parent_map.find(node);
        if (node_it == parent_map.end()) {
          break;
        } else if (node_it->second == node) {
          break;
        } else {
          node = node_it->second;
        }
      }

      if (score > max_score) {
        max_taxa.clear();
        max_score = score;
        max_taxon = taxon;
      }
      else if (score == max_score) {
        if (max_taxa.empty())
          max_taxa.insert(max_taxon);
        max_taxa.insert(taxon);
      }
    }

    // If two LTR paths are tied for max, return LCA of all
    if (max_taxa.size() > 1) {
      set<uint32_t>::iterator sit = max_taxa.begin();
      max_taxon = *sit;
      for (sit++; sit != max_taxa.end(); sit++)
        max_taxon = lca(parent_map, max_taxon, *sit);
    } else if(max_taxa.size() == 1){
      set<uint32_t>::iterator sit = max_taxa.begin();
      max_taxon = *sit;
    }
    return max_taxon;
  }


  /*
   * Return the best unambigous promotion present in the barcode.
   */
  uint32_t promote_call(uint32_t call, uint32_t max_hops,
                        const unordered_map<uint32_t, uint32_t> &parent_map,
                        const unordered_map<uint32_t, uint32_t> &bc_hit_counts){
    uint32_t taxon, parent, hops_taken;
    unordered_map<uint32_t, vector<uint32_t>> bc_child_map;
    vector<uint32_t> child_vec;

    for(auto it=bc_hit_counts.begin(); it!= bc_hit_counts.end(); ++it){
      taxon = it->first;
      auto parent_node = parent_map.find(taxon);
      if(parent_node != parent_map.end()){
	      parent = parent_node->second;
	      bc_child_map[parent].push_back(taxon);
      }
    }
    hops_taken = 0;
    while (hops_taken < max_hops) {
      if(bc_child_map.count(call) == 0) // call already is a leaf
        return call;
      child_vec  = bc_child_map[call];
      if(child_vec.size() >= 2) // ambiguous promotion
        return call;
      call = child_vec.front();
      hops_taken++;
    }
  }


  uint8_t KmerScanner::k = 0;
  uint64_t KmerScanner::kmer_mask = 0;
  uint32_t KmerScanner::mini_kmer_mask = 0;

  // Create a scanner for the string over the interval [start, finish)
  KmerScanner::KmerScanner(string &seq, size_t start, size_t finish) {
    if (! k)
      errx(EX_SOFTWARE, "KmerScanner created w/o setting k");
    if (finish > seq.size())
      finish = seq.size();

    kmer = 0;
    ambig = 0;
    str = &seq;
    curr_pos = start;
    pos1 = start;
    pos2 = finish;
    loaded_nt = 0;
    if (pos2 - pos1 + 1 < k)
      curr_pos = pos2;
  }

  uint8_t KmerScanner::get_k() { return k; }

  void KmerScanner::set_k(uint8_t n) {
    if (k)  // Only allow one setting per execution
      return;
    k = n;
    kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (k * 2);
    mini_kmer_mask = ~0;
    mini_kmer_mask >>= sizeof(mini_kmer_mask) * 8 - k;
  }

  uint64_t *KmerScanner::next_kmer() {
    bool skip_pos = false;
    if (curr_pos >= pos2)
      return NULL;
    if (loaded_nt)  
      loaded_nt--;
    while (loaded_nt < k) {
      if (skip_pos) {
	skip_pos = false;
      } else {
        loaded_nt++;
        kmer <<= 2;
        ambig <<= 1;
      }
      switch ((*str)[curr_pos++]) {
        case 'A': case 'a':
          break;
        case 'C': case 'c':
          kmer |= 1;
          break;
        case 'G': case 'g':
          kmer |= 2;
          break;
        case 'T': case 't':
          kmer |= 3;
          break;
	case '\n': case '\r':
	  --loaded_nt;
	  skip_pos = true;
	  continue;
	  break;
        default:
          ambig |= 1;
          break;
      }
      kmer &= kmer_mask;
      ambig &= mini_kmer_mask;
    }
    return &kmer;
  }

  bool KmerScanner::ambig_kmer() {
    return !! ambig;
  }
}
