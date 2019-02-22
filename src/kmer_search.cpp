/*
Copyright (c) 2006-2018 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.
SINA is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

SINA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with SINA.  If not, see <http://www.gnu.org/licenses/>.

Additional permission under GNU GPL version 3 section 7

If you modify SINA, or any covered work, by linking or combining it
with components of ARB (or a modified version of that software),
containing parts covered by the terms of the
ARB-public-library-license, the licensors of SINA grant you additional
permission to convey the resulting work. Corresponding Source for a
non-source form of such a combination shall include the source code
for the parts of ARB used as well as that of the covered work.
*/

#include "kmer_search.h"
#include "kmer.h"
#include "idset.h"
#include "query_arb.h"
#include "helpers.h"
#include "timer.h"
#include "log.h"
#include "progress.h"

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <iostream>

#include <unordered_map>
#include <unordered_set>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/thread/mutex.hpp>

#include <tbb/tbb.h>
#include <tbb/cache_aligned_allocator.h>

#include <cstdio>
#include <sys/stat.h>

#include "zlib.h"

using namespace sina;

static const char* module_name = "Search (internal)";
static auto logger = Log::create_logger(module_name);


class kmer_search::impl {
public:
    unsigned int k;
    unsigned int n_kmers;
    unsigned int n_sequences{0};

    std::vector<std::string> sequence_names;
    std::vector<vlimap*> kmer_idx;

    query_arb* arbdb;
    timer_mt timeit;

    impl(query_arb* arbdb_, int k_);
    ~impl() {
        logger->info("Timings for Kmer Search: {}", timeit);
    }
    void find(const cseq& query, std::vector<cseq>& results, int max);
    void build();
    void store(const fs::path& filename);
    bool try_load(const fs::path& filename);
};


using kmer_search_key_t = std::pair<fs::path, int>;
static std::map<kmer_search_key_t,
                std::shared_ptr<kmer_search::impl>> indices;
static boost::mutex indices_access;

kmer_search*
kmer_search::get_kmer_search(const fs::path& filename, int k) {
    kmer_search_key_t key(filename, k);
    if (indices.count(key) == 0u) {
        std::shared_ptr<kmer_search::impl> pimpl(new impl(query_arb::getARBDB(filename), k));
        boost::mutex::scoped_lock lock(indices_access);
        indices[key] = pimpl;
    }
    return new kmer_search(indices[key]);
}

void
kmer_search::release_kmer_search(const fs::path& filename, int k) {
    kmer_search_key_t key(filename, k);
    boost::mutex::scoped_lock lock(indices_access);
    indices.erase(key);
}

kmer_search::kmer_search(std::shared_ptr<kmer_search::impl> pimpl_)
    : pimpl(pimpl_) {}
kmer_search::~kmer_search() = default;


class IndexBuilder {
    kmer_search::impl *idx;
    Progress *p;
public:
    std::vector<vlimap*> kmer_idx;

    void operator()(const tbb::blocked_range<size_t>& r) {
        size_t end = r.end();
        std::unordered_set<unsigned int> seen;
        for (size_t i = r.begin(); i < end; ++i) {
            const cseq& c = idx->arbdb->getCseqUncached(idx->sequence_names[i]);
            const auto& bases = c.const_getAlignedBases();
            for (unsigned int kmer: unique_prefix_kmers(bases, seen, idx->k, 1, BASE_A)) {
                //for (const auto& kmer: unique_kmers(bases, seen, idx->k)) {
                if (kmer_idx[kmer] == nullptr) {
                    kmer_idx[kmer] = new vlimap(idx->n_sequences);
                }
                kmer_idx[kmer]->push_back(i);
            }
            p->update();
        }
    }

    void join(IndexBuilder& other) {
        for (unsigned int i = 0; i < idx->n_kmers; ++i) {
            if (other.kmer_idx[i] != nullptr) {
                if (kmer_idx[i] != nullptr) {
                    kmer_idx[i]->append(*other.kmer_idx[i]);
                } else {
                    kmer_idx[i] = other.kmer_idx[i];
                    other.kmer_idx[i] = nullptr;
                }
            }
        }
    }

    ~IndexBuilder() {
        for (unsigned int i = 0; i < idx->n_kmers; ++i) {
            delete kmer_idx[i];
        }
    }

    IndexBuilder(IndexBuilder& x, tbb::split)
        : idx(x.idx), p(x.p), kmer_idx(x.idx->n_kmers, nullptr)
    {
    }

    IndexBuilder(kmer_search::impl *idx_, Progress *p_)
        : idx(idx_), p(p_), kmer_idx(idx->n_kmers, nullptr)
    {
    }
};

kmer_search::impl::impl(query_arb* arbdb_, int k_)
    : k(k_),
      n_kmers(1<<(k_*2)),
      kmer_idx(1<<(k_*2), nullptr),
      arbdb(arbdb_)
{
    fs::path dbpath = arbdb->getFileName();
    fs::path idxpath = fs::path(dbpath).replace_extension("sidx");
    if (fs::exists(idxpath) && fs::exists(dbpath)) {
        if (fs::last_write_time(idxpath) >= fs::last_write_time(dbpath)) {
            if (try_load(idxpath)) {
                return;
            } else {
                logger->error("Failed to load {} - rebuilding", idxpath);
            }
            logger->error("Reference {} newer than {} - rebuilding index.",
                          dbpath, idxpath);
        }
    }
    build();
    store(idxpath);
}

void
kmer_search::impl::build() {
    timestamp start;

    sequence_names = arbdb->getSequenceNames();
    n_sequences = sequence_names.size();

    Progress p("Building Index", n_sequences);

    IndexBuilder bi(this, &p);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n_sequences), bi);
    p += n_sequences - p.count();

    kmer_idx.clear();
    kmer_idx.reserve(n_kmers);
    int total = 0;
    p.restart("Compressing", n_kmers);
    for (unsigned int i=0; i < n_kmers; i++) {
        ++p;
        if (bi.kmer_idx[i] != nullptr) {
            if (bi.kmer_idx[i]->size() > n_sequences / 2) {
                bi.kmer_idx[i]->invert();
            }
            total += bi.kmer_idx[i]->size();
            kmer_idx.push_back(new vlimap(*bi.kmer_idx[i]));
        } else {
            kmer_idx.push_back(nullptr);
        }
    }

    logger->info("Built index from {} sequences ({} refs) in {}",
                 n_sequences, total, timestamp()-start);
}

const uint64_t idx_header_magic = 0x5844494b414e4953; // SINAKIDX
const uint16_t idx_header_vers  = 0;
struct idx_header {
    uint64_t magic{idx_header_magic};
    uint16_t vers{idx_header_vers};
    uint16_t k;
    uint32_t n_sequences;
};

void
kmer_search::impl::store(const fs::path& filename) {
    std::ofstream out(filename.native(), std::ofstream::binary);
    idx_header header;
    header.k = k;
    header.n_sequences = n_sequences;
    out.write((char*)&header, sizeof(idx_header));
    for (auto& name : sequence_names) {
        out << name << std::endl;
    }
    vlimap emptymap(n_sequences);
    for (unsigned int i = 0; i < n_kmers; i++) {
        if (kmer_idx[i] != nullptr && kmer_idx[i]->size() > 0) {
            emptymap.push_back(i);
        }
    }
    emptymap.write(out);

    size_t idxno = 0;
    for (auto inc : emptymap) {
        idxno += inc;
        kmer_idx[idxno]->write(out);
    }
}

bool
kmer_search::impl::try_load(const fs::path& filename) {
    std::ifstream in(filename.native(), std::ifstream::binary);
    idx_header header;
    in.read((char*)&header, sizeof(idx_header));
    if (k != header.k) {
        return false;
    }
    n_sequences = header.n_sequences;
    for (unsigned int i = 0; i < n_sequences; i++) {
        string name;
        getline(in, name);
        sequence_names.push_back(name);
    }
    vlimap emptymap(n_sequences);
    emptymap.read(in);
    size_t idxno = 0;
    int total = 0;
    for (auto inc : emptymap) {
        idxno += inc;
        vlimap *idx = new vlimap(n_sequences);
        idx->read(in);
        kmer_idx[idxno] = idx;
        total += idx->size();
    }
    logger->info("Index contains {} sequences ({} refs)", n_sequences, total);

    return true;
}

double
kmer_search::match(std::vector<cseq>& results,
                   const cseq& query,
                   int  /*min_match*/,
                   int max_match,
                   float  /*min_score*/,
                   float  /*max_score*/,
                   query_arb*  /*arb*/,
                   bool  /*noid*/,
                   int  /*min_len*/,
                   int  /*num_full*/,
                   int  /*full_min_len*/,
                   int  /*range_cover*/,
                   bool  /*leave_query_out*/) {
    find(query, results, max_match);
    if (results.empty()) {
        return 0;
    }
    return results[0].getScore();
}

void
kmer_search::find(const cseq& query, std::vector<cseq>& results, int max) {
    pimpl->find(query, results, max);
}

void
kmer_search::impl::find(const cseq& query, std::vector<cseq>& results, int max) {
    if (n_sequences == 0) {
        return;
    }
    timer& timing = timeit.get_timer();
    timing.start();
    const vector<aligned_base>& bases = query.const_getAlignedBases();
    idset::inc_t scores(n_sequences, 0);
    timing.stop("load query");

    std::unordered_set<unsigned int> seen(query.size()*2-1);

    int offset = 0;
    // for (unsigned int kmer: all_kmers(bases, k)) {
    // for (unsigned int kmer: unique_kmers(bases, seen, k)) {
    // for (unsigned int kmer: unique_prefix_kmers(bases, seen, k, 1, BASE_A)) {
    for (unsigned int kmer: prefix_kmers(bases, k, 1, BASE_A)) {
        if (kmer_idx[kmer] != nullptr) {
            offset += kmer_idx[kmer]->increment(scores);
        }
    }
    timing.stop("count kmers");

    using pair = std::pair<idset::inc_t::value_type, int>;
    std::vector<pair> ranks;
    ranks.reserve(n_sequences);
    int n = 0;
    for (auto score: scores) {
        ranks.emplace_back(score + offset, n++);
    }
    std::partial_sort(ranks.begin(), ranks.begin()+max, ranks.end(),
                      std::greater<pair>());
    timing.stop("sort result");

    results.clear();
    results.reserve(max);
    for (int i=0; i<max; i++) {
        results.emplace_back(arbdb->getCseq(sequence_names[ranks[i].second]));
        results.back().setScore(ranks[i].first);
    }
    timing.stop("load result");
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
