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

#include "search_filter.h"
#include "config.h"
#include "log.h"
#include "progress.h"
#include "famfinder.h"
#include "kmer_search.h"

#include <string>
using std::string;

#include <utility>
#include <vector>
using std::vector;

#include <map>
using std::map;

#include <iostream>

#include <sstream>
using std::stringstream;

#include <exception>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/algorithm/string.hpp>

#include "alignment_stats.h"
#include "cseq_comparator.h"

namespace sina{

static auto logger = Log::create_logger("search");

struct search_filter::options {
    fs::path pt_database;
    ENGINE_TYPE engine;
    string pt_port;
    bool search_all;
    string posvar_filter;
    bool fs_no_fast;
    int fs_kmer_len;
    int fs_kmer_mm;
    bool fs_kmer_norel;
    int kmer_candidates;
    float min_sim;
    bool ignore_super;
    int max_result;
    string lca_fields;
    vector<string> v_lca_fields;
    float lca_quorum;
    string copy_fields;
    vector<string> v_copy_fields;
    cseq_comparator comparator;
};

struct search_filter::options *search_filter::opts;

void
search_filter::get_options_description(po::options_description& main,
                                       po::options_description& adv){
    opts = new struct search_filter::options();

    po::options_description mid("Search & Classify");
    mid.add_options()
        ("search-db", po::value<fs::path>(&opts->pt_database),
         "reference db if different from -r/--db")
        ("search-engine", po::value<ENGINE_TYPE>(&opts->engine),
         "engine if different from --fs-engine")
        ("search-min-sim", po::value<float>(&opts->min_sim)->default_value(.7, ""),
         "required sequence similarity (0.7)")
        ("search-max-result", po::value<int>(&opts->max_result)->default_value(10, ""),
         "desired number of search results (10)")
        ("lca-fields", po::value<string>(&opts->lca_fields)->default_value(""),
         "names of fields containing source taxonomy (colon separated list)")
        ("lca-quorum", po::value<float>(&opts->lca_quorum)->default_value(.7, ""),
         "fraction of search result that must share resulting classification (0.7)")
        ;
    main.add(mid);

    po::options_description od("Search & Classify");
    od.add_options()
        ("search-port", po::value<string>(&opts->pt_port)
         ->default_value(fmt::format(":/tmp/sina_pt2_{}", getpid()), ""),
         "PT server port (:/tmp/sina_pt2_<pid>)")
        ("search-all", po::bool_switch(&opts->search_all), "do not use k-mer heuristic")
        ("search-no-fast", po::bool_switch(&opts->fs_no_fast), "don't use fast family search")
        ("search-kmer-candidates", po::value<int>(&opts->kmer_candidates)->default_value(1000,""),         "number of most similar sequences to acquire via kmer-step (1000)")
        ("search-kmer-len", po::value<int>(&opts->fs_kmer_len)->default_value(10, ""),
         "length of k-mers (12)")
        ("search-kmer-mm", po::value<int>(&opts->fs_kmer_mm)->default_value(0,""),
         "allowed mismatches per k-mer (0)")
        ("search-kmer-norel", po::bool_switch(&opts->fs_kmer_norel),
         "don't score k-mer distance relative to target length")
        ("search-ignore-super", po::bool_switch(&opts->ignore_super),
         "ignore sequences containing query")
        ("search-copy-fields", po::value<string>(&opts->copy_fields)->default_value(""),
         "copy fields from result sequences to query sequence")
        ;
    od.add(cseq_comparator::get_options_description("search-"));

    adv.add(od);
}


void
search_filter::validate_vm(boost::program_options::variables_map& vm,
                           po::options_description&  /*desc*/) {
    if (opts == nullptr) {
        throw std::logic_error("search options not parsed?!");
    }
 
    // we need a DB to search
    if (opts->pt_database.empty()) {
        if ((vm.count("db") != 0u) && !vm["db"].as<fs::path>().empty()) {
            opts->pt_database = vm["db"].as<fs::path>();
            if (vm["search-port"].defaulted()) {
                opts->pt_port = vm["ptport"].as<string>();
            }
        } else {
          throw std::logic_error("need search-db to search");
        }
    }

    if (vm.count("search-engine") == 0u) {
        opts->engine = famfinder::get_engine();
    }

    opts->comparator = cseq_comparator::make_from_variables_map(vm, "search-");

    boost::split(opts->v_lca_fields, opts->lca_fields,
                 boost::is_any_of(":,"));
    if (opts->v_lca_fields.back().empty()) {
        opts->v_lca_fields.pop_back();
    }
    boost::split(opts->v_copy_fields, opts->copy_fields,
                 boost::is_any_of(":,"));
    if (opts->v_copy_fields.back().empty()) {
        opts->v_copy_fields.pop_back();
    }
}

} // namespace sina

using namespace sina;

struct search_filter::priv_data {
    std::unique_ptr<search> index;
    query_arb *arb;
    vector<const cseq*> sequences;
};

search_filter::search_filter()
    :  data(new priv_data)
{
    data->arb = query_arb::getARBDB(opts->pt_database);

    if (opts->search_all) {
        Progress p("Caching Sequences", data->arb->getSeqCount());
        vector<string> names = data->arb->getSequenceNames();
        for (string& name: names) {
            data->sequences.push_back(&data->arb->getCseq(name));
            ++p;
        }
    } else if (opts->engine == ENGINE_ARB_PT) {
        data->index = std::unique_ptr<search>(
            query_pt_pool::get_pool(
                opts->pt_database,
                opts->fs_kmer_len,
                not opts->fs_no_fast,
                opts->fs_kmer_norel,
                opts->fs_kmer_mm,
                opts->pt_port
                ));
    } else if (opts->engine == ENGINE_SINA_KMER) {
        data->index = std::unique_ptr<search>(
            kmer_search::get_kmer_search(
                opts->pt_database,
                opts->fs_kmer_len,
                opts->fs_no_fast
                ));
    } else {
        throw std::runtime_error("Unknown engine");
    }
}

search_filter::search_filter(const search_filter& o) = default;
search_filter& search_filter::operator=(const search_filter& o) = default;
search_filter::~search_filter() = default;

template<typename F>
struct dereference {
    explicit dereference(F f) : _f(std::move(f)){}
    dereference() : _f(){}
    using result_type = typename F::result_type;

    template<typename A, typename B>
    result_type operator()(const A& a,
                           const B& b) {
        return _f(*a, *b);
    }

    template<typename A>
    result_type operator()(const A& a) {
        return _f(*a);
    }
    F _f;
};

tray
search_filter::operator()(tray t) {
    cseq *c;
    if (t.aligned_sequence != nullptr) {
        c = t.aligned_sequence;
    } else {
        t.log << "search: no sequence?!;";
        return t;
    }

    if (c->size() < 20) {
        t.log << "search:sequence too short (<20 bases);";
        return t;
    }

    t.search_result = new search::result_vector();
    auto& vc = *t.search_result;


    auto iupac_compare = [](const aligned_base& a, const aligned_base& b) {
        return a.comp(b);
    };
    auto contains_query = [&](search::result_item& item) {
        return boost::algorithm::contains(item.sequence->getAlignedBases(),
                                          c->getAlignedBases(),
                                          iupac_compare);
    };

    string bases = c->getBases();
    if (opts->search_all) {
        search::result_vector result;
        result.reserve(data->sequences.size());
        for (const cseq *r: data->sequences) {
            result.emplace_back(opts->comparator(*c, *r), r);
        }
        auto it = result.begin();
        auto middle = it + opts->max_result;
        auto end = result.end();
        
        do {
            partial_sort(it, middle, end, std::greater<search::result_item>());
            if (opts->ignore_super) {
                // sort sequences containing query to beginning
                // reset "start" to beginnin of non-identical section
                middle = it + opts->max_result;
                it = partition(it, middle, contains_query);
            }
            // repeat if we removed some just above
        } while (middle != end && it + opts->max_result > middle); 

        while (it != middle && it->score > opts->min_sim) {
            vc.push_back(*it);
            ++it;
        }
    } else {
        data->index->find(*c, vc, opts->kmer_candidates);
        int i=0;
        for (auto &r: vc) { //FIXME: remove bug tracing here
            ++i;
            if (r.sequence->getAlignedBases().data() == nullptr) {
                logger->error("BUG {} {}", *c, i);
                logger->error("  {}", *r.sequence);
            }
        }

        if (opts->ignore_super) {
            auto it = partition(vc.begin(), vc.end(), contains_query);
            vc.erase(it, vc.end());
        }

        for (auto& r: vc) {
            r.score = opts->comparator(*c, *r.sequence);
        }

        auto it = vc.begin();
        auto middle = vc.begin() + opts->max_result;
        auto end=vc.end();

        if (middle > end) {
            middle = end;
        }

        partial_sort(it, middle, end, std::greater<search::result_item>());

        while (it != middle && it->score > opts->min_sim) {
            ++it;
        }

        vc.erase(it, vc.end());
    }

    fmt::memory_buffer nearest;
    map<string, vector<vector<string> > > group_names_map;
    for (auto &i: vc) {
        const cseq& r = *i.sequence;
        data->arb->loadKey(r, "acc");
        data->arb->loadKey(r, "version");
        data->arb->loadKey(r, "start");
        data->arb->loadKey(r, "stop");

        for (string& s: opts->v_lca_fields) {
            data->arb->loadKey(r, s);
            string tax_path = r.get_attr<string>(s);
            if (tax_path == "Unclassified;") {
                continue;
            }
            vector<string> group_names;
            boost::split(group_names, tax_path, boost::is_any_of(";"));
            if (group_names.back().empty() || group_names.back()  == " ") {
                group_names.pop_back();
            }
            group_names_map[s].push_back(group_names);
        }

        fmt::format_to(nearest,
                       "{}.{}.{}.{}~{:.3f} ",
                       r.get_attr<string>("acc"),
                       r.get_attr<string>("version"),
                       r.get_attr<string>("start"),
                       r.get_attr<string>("stop"),
                       i.score);

        string acc = r.get_attr<string>("acc");
        for (string& s: opts->v_copy_fields) {
            data->arb->loadKey(r,s);
            string v = r.get_attr<string>(s);
            c->set_attr<string>(string("copy_")+acc+string("_")+s, v);
        }
    }
    c->set_attr<string>(query_arb::fn_nearest, fmt::to_string(nearest));

    for (string& s: opts->v_lca_fields) {
        stringstream result;
        vector<vector<string> > group_names = group_names_map[s];
        for (vector<string>& vs: group_names) {
            reverse(vs.begin(), vs.end());
        }
        int outliers = vc.size() * (1 - opts->lca_quorum);
        while (outliers >= 0 && !group_names.empty()) {
            auto it = group_names.begin();
            if (it->empty()) {
                group_names.erase(it);
                outliers--;
                continue;
            }
            string name = it->back();
            ++it;
            for (; it != group_names.end(); ++it) {
                if (it->empty() || it->back() != name) {
                    break;
                }
            }
            if (it != group_names.end()) {
                group_names.erase(it);
                outliers--;
                continue;
            }
            for (vector<string>& vs: group_names) {
                vs.pop_back();
            }
            result << name << ";";
        }
        string res = result.str();
        if (res.size()>1 && res.substr(res.size()-2) == ";;") {
            res = res.substr(res.size()-1);
        }
        if (res.empty() || res == ";") {
            res = "Unclassified;";
        }
        c->set_attr<string>(string("lca_")+s, res);
    }

    return t;
}


/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
