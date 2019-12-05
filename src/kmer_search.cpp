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
#include "cache.h"

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <iostream>

#include <unordered_map>
#include <unordered_set>
#include <mutex>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <tbb/tbb.h>
#include <tbb/cache_aligned_allocator.h>
#include <cstdio>
#include <sys/stat.h>

#include "zlib.h"

using namespace sina;

static const char* module_name = "Search (internal)";
static auto logger = Log::create_logger(module_name);

const uint64_t idx_header_magic = 0x5844494b414e4953; // SINAKIDX
const uint16_t idx_header_vers  = 0;
union idx_flags {
    uint16_t flags{0};
    struct {
        uint16_t k :8;
        uint16_t nofast :1;
        uint16_t reserved :7;
    };
    // must be sortable for use as map key below
    bool operator<(const idx_flags& r) const {
        if (k < r.k) return true;
        if (k > r.k) return false;
        return nofast < r.nofast;
    }
};

struct idx_header {
    uint64_t magic{idx_header_magic};
    uint16_t vers{idx_header_vers};
    uint32_t n_sequences{0};
    idx_flags flags{0};
};

class kmer_search::impl {
public:
    unsigned int k;
    unsigned int n_kmers;
    unsigned int n_sequences{0};

    bool nofast;

    std::vector<std::string> sequence_names;
    std::vector<vlimap*> kmer_idx;

    query_arb* arbdb;
    timer_mt timeit;

    using rank_result_type = std::vector<std::pair<idset::inc_t::value_type, int>>;
    fifo_cache<std::string, rank_result_type> cache{32};

    impl(query_arb* arbdb_, int k_, bool nofast_);
    ~impl() {
        logger->info("Timings for Kmer Search: {}", timeit);
    }
    void find(const cseq& query, result_vector& results, unsigned int max);
    void build();
    void store(const fs::path& filename);
    bool try_load(const fs::path& filename);
};


using kmer_search_key_t = std::pair<fs::path, idx_flags>;
static std::map<kmer_search_key_t, std::shared_ptr<kmer_search::impl>> indices;
static std::mutex indices_access;

kmer_search*
kmer_search::get_kmer_search(const fs::path& filename, int k, bool nofast) {
    idx_flags flags;
    flags.k = k;
    flags.nofast = nofast;
    kmer_search_key_t key(filename, flags);
    if (indices.count(key) == 0u) {
        std::shared_ptr<kmer_search::impl> pimpl(new impl(query_arb::getARBDB(filename), k, nofast));
        std::lock_guard<std::mutex> lock(indices_access);
        indices[key] = pimpl;
    }
    return new kmer_search(indices[key]);
}

void
kmer_search::release_kmer_search(const fs::path& filename, int k, bool nofast) {
    idx_flags flags;
    flags.k = k;
    flags.nofast = nofast;
    kmer_search_key_t key(filename, flags);
    std::lock_guard<std::mutex> lock(indices_access);
    indices.erase(key);
}

kmer_search::kmer_search(std::shared_ptr<kmer_search::impl> pimpl_)
    : pimpl(pimpl_) {}
kmer_search::~kmer_search() = default;
unsigned int kmer_search::size() const { return pimpl->n_sequences; }


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
            const auto& bases = c.getAlignedBases();
            if (idx->nofast) {
                for (const auto& kmer: unique_kmers(bases, seen, idx->k)) {
                    if (kmer_idx[kmer] == nullptr) {
                        kmer_idx[kmer] = new vlimap(idx->n_sequences);
                    }
                    kmer_idx[kmer]->push_back(i);
                }
            } else {
                for (unsigned int kmer: unique_prefix_kmers(bases, seen, (int)idx->k, 1, BASE_A)) {
                    if (kmer_idx[kmer] == nullptr) {
                        kmer_idx[kmer] = new vlimap(idx->n_sequences);
                    }
                    kmer_idx[kmer]->push_back(i);
                }
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

kmer_search::impl::impl(query_arb* arbdb_, int k_, bool nofast_)
    : k(k_),
      n_kmers(1<<(k_*2)),
      kmer_idx(1<<(k_*2), nullptr),
      nofast(nofast_),
      arbdb(arbdb_)
{
    fs::path dbpath = arbdb->getFileName();
    if (dbpath.compare(":") == 0) {
        logger->warn("Remote database found. Building in memory index.");
        build();
        return;
    }

    fs::path idxpath = fs::path(dbpath).replace_extension("sidx");
    if (fs::exists(idxpath) && fs::exists(dbpath)) {
        if (fs::last_write_time(idxpath) >= fs::last_write_time(dbpath)) {
            if (try_load(idxpath)) {
                return;
            }
        } else {
            logger->warn("Reference {} newer than {}", dbpath, idxpath);
        }
        logger->warn("Failed to load {} - rebuilding", idxpath);
    } else {
        logger->warn("No cached index found.");
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

void
kmer_search::impl::store(const fs::path& filename) {
    std::string native = filename.native();
    std::ofstream out(native, std::ofstream::binary);
    idx_header header;
    header.flags.k = k;
    header.flags.nofast = nofast;
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
    std::string native = filename.native();
    std::ifstream in(native, std::ifstream::binary);
    idx_header header;
    in.read((char*)&header, sizeof(idx_header));
    if (idx_header_magic != header.magic) {
        logger->error("Index file {} has wrong magic. Aborting.",
                      filename);
        exit(1);
    }
    if (idx_header_vers != header.vers) {
        logger->error("Index file {} created by different version.",
                      filename);
        return false;
    }
    if (k != header.flags.k) {
        logger->error("Index file {} build for k={} not k={}",
                      filename, header.flags.k, k);
        return false;
    }
    if (nofast != header.flags.nofast) {
        logger->error("Index file {} build for {} (want {})",
                      nofast?"no fast":"fast", nofast?"fast":"nofast");
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
kmer_search::match(result_vector&, const cseq&, int, int, float, float, query_arb*,
                   bool, int, int, int, int, bool) {
    throw std::runtime_error("Legacy family composition not implemented for internal search");
    return 0;
}

void
kmer_search::find(const cseq& query, result_vector& results, unsigned int max) {
    pimpl->find(query, results, max);
}

void
kmer_search::impl::find(const cseq& query, result_vector& results, unsigned int max) {
    if (max > n_sequences) {
        max = n_sequences;
    }
    if (max == 0) {
        return;
    }
    idset::inc_t scores(n_sequences, 0);
    using pair = std::pair<idset::inc_t::value_type, int>;
    std::vector<pair> ranks;

    std::string bases = query.getBases();
    if (!cache.try_get(bases, ranks)) {
        timer& timing = timeit.get_timer();
        timing.start();
        const vector<aligned_base>& bases = query.getAlignedBases();

        timing.stop("load query");

        std::unordered_set<unsigned int> seen(query.size()*2-1);

        int offset = 0;
        if (nofast) {
            // for (unsigned int kmer: unique_kmers(bases, seen, k)) {
            for (unsigned int kmer: all_kmers(bases, k, 1)) {
                if (kmer_idx[kmer] != nullptr) {
                    offset += kmer_idx[kmer]->increment(scores);
                }
            }
        } else { // fast
            // for (unsigned int kmer: unique_prefix_kmers(bases, seen, k, 1, BASE_A)) {
            for (unsigned int kmer: prefix_kmers(bases, k, 1, BASE_A)) {
                if (kmer_idx[kmer] != nullptr) {
                    offset += kmer_idx[kmer]->increment(scores);
                }
            }
        }
        timing.stop("count kmers");

        ranks.reserve(n_sequences);
        int n = 0;
        for (auto score: scores) {
            ranks.emplace_back(score + offset, n++);
        }
        timing.stop("store");
    }
    std::partial_sort(ranks.begin(), ranks.begin()+max, ranks.end(),  std::greater<pair>());

    results.clear();
    results.reserve(max);
    for (unsigned int i=0; i<max; i++) {
        results.emplace_back(ranks[i].first, &arbdb->getCseq(sequence_names[ranks[i].second]));
    }
    cache.store(std::move(bases), std::move(ranks));
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
