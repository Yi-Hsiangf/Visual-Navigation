/**
BSD 3-Clause License

Copyright (c) 2018, Vladyslav Usenko and Nikolaus Demmel.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <fstream>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>

#include <visnav/common_types.h>
#include <visnav/serialization.h>

using namespace std;
namespace visnav {

bool comp(const pair<FrameCamId, double>& a,
          const pair<FrameCamId, double>& b) {
  return a.second < b.second;
}

class BowDatabase {
 public:
  BowDatabase() {}

  inline void insert(const FrameCamId& fcid, const BowVector& bow_vector) {
    // TODO SHEET 3: add a bow_vector that corresponds to frame fcid to the
    // inverted index. You can assume the image hasn't been added before.
    UNUSED(fcid);
    UNUSED(bow_vector);

    for (size_t i = 0; i < bow_vector.size(); i++) {
      WordId wordid = bow_vector[i].first;
      WordValue wordvalue = bow_vector[i].second;
      std::pair<FrameCamId, WordValue> inverted(fcid, wordvalue);
      inverted_index[wordid].push_back(inverted);
    }
  }

  inline void query(const BowVector& bow_vector, size_t num_results,
                    BowQueryResult& results) const {
    // TODO SHEET 3: find num_results closest matches to the bow_vector in the
    // inverted index. Hint: for good query performance use std::unordered_map
    // to accumulate scores and std::partial_sort for getting the closest
    // results. You should use L1 difference as the distance measure. You can
    // assume that BoW descripors are L1 normalized.
    UNUSED(bow_vector);   // known : img bow_vector
    UNUSED(num_results);  // known
    UNUSED(results);      // unknown BowQueryResult =
                          // std::vector<std::pair<FrameCamId, double>>;
    double score = 0;
    std::unordered_map<FrameCamId, double> frame_score_map;
    for (size_t j = 0; j < bow_vector.size(); j++) {
      WordId wordid = bow_vector[j].first;
      WordValue weight = bow_vector[j].second;
      // tbb::concurrent_unordered_map <WordId,
      // tbb::concurrent_vector<std::pair<FrameCamId, WordValue>>>;
      for (size_t frame_idx = 0; frame_idx < inverted_index.at(wordid).size();
           frame_idx++) {
        FrameCamId tree_frame_id = inverted_index.at(wordid)[frame_idx].first;
        WordValue tree_frame_weight =
            inverted_index.at(wordid)[frame_idx].second;
        score = abs(weight - tree_frame_weight) - abs(weight) -
                abs(tree_frame_weight);

        if (frame_score_map.find(tree_frame_id) ==
            frame_score_map.end())  // not found: initial
          frame_score_map[tree_frame_id] = 2 + score;
        else  // found
          frame_score_map[tree_frame_id] += score;
      }
    }

    vector<pair<FrameCamId, double>> frame_score_vec{frame_score_map.begin(),
                                                     frame_score_map.end()};
    int s = frame_score_vec.size();
    // cout << "num result" << num_results << endl;
    if (s < num_results) num_results = s;
    std::partial_sort(frame_score_vec.begin(),
                      frame_score_vec.begin() + num_results,
                      frame_score_vec.end(), comp);

    // std::cout << "score: " << std::endl;
    for (size_t i = 0; i < num_results; i++) {
      results.push_back(frame_score_vec[i]);
      // std::cout << frame_score_vec[i].second << " " << std::endl;
    }
  }

  void clear() { inverted_index.clear(); }

  void save(const std::string& out_path) {
    BowDBInverseIndex state;
    for (const auto& kv : inverted_index) {
      for (const auto& a : kv.second) {
        state[kv.first].emplace_back(a);
      }
    }
    std::ofstream os;
    os.open(out_path, std::ios::binary);
    cereal::JSONOutputArchive archive(os);
    archive(state);
  }

  void load(const std::string& in_path) {
    BowDBInverseIndex inverseIndexLoaded;
    {
      std::ifstream os(in_path, std::ios::binary);
      cereal::JSONInputArchive archive(os);
      archive(inverseIndexLoaded);
    }
    for (const auto& kv : inverseIndexLoaded) {
      for (const auto& a : kv.second) {
        inverted_index[kv.first].emplace_back(a);
      }
    }
  }

  const BowDBInverseIndexConcurrent& getInvertedIndex() {
    return inverted_index;
  }

 protected:
  BowDBInverseIndexConcurrent inverted_index;
};  // namespace visnav

}  // namespace visnav
