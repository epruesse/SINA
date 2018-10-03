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

#include <string>
using std::string;

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
#include <boost/progress.hpp>
#include <boost/algorithm/string.hpp>
using boost::algorithm::iequals;

#include "alignment_stats.h"
#include "cseq_comparator.h"

namespace sina{

static auto logger = Log::create_logger("search");

struct search_filter::options {
    string pt_database;
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

    //main.add_options()

    po::options_description mid("Search & Classify");
    mid.add_options()
        ("search-db",
         po::value<string>(&opts->pt_database),
         "reference db")
        ("search-min-sim",
         po::value<float>(&opts->min_sim)->default_value(.7, ""),
         "required sequence similarity (0.7)")
        ("search-max-result",
         po::value<int>(&opts->max_result)->default_value(10, ""),
         "desired number of search results (10)")
        ("lca-fields",
         po::value<string>(&opts->lca_fields)->default_value(""),
         "names of fields containing source taxonomy (colon separated list)")
        ("lca-quorum",
         po::value<float>(&opts->lca_quorum)->default_value(.7, ""),
         "fraction of search result that must share resulting classification (0.7)")
        ;
    main.add(mid);

    po::options_description od("Search & Classify");
    od.add_options()
        ("search-port",
#ifdef HAVE_GETPID
         po::value<string>(&opts->pt_port)
         ->default_value(":/tmp/sina_pt2_"
                         + boost::lexical_cast<std::string>(getpid()),""),
#else
         po::value<string>(&opts->pt_port)->default_value("localhost:4041"),
#endif
         "PT server port")
        ("search-all",
         po::bool_switch(&opts->search_all),
         "do not use k-mer heuristic")
//        ("search-filter",
//         po::value<string>(&opts->posvar_filter)->default_value("none"),
//         "select posvar filter")
        ("search-no-fast",
         po::bool_switch(&opts->fs_no_fast),
         "don't use fast family search")
        ("search-kmer-candidates",
         po::value<int>(&opts->kmer_candidates)->default_value(1000,""),
         "number of most similar sequences to acquire via kmer-step (1000)")
        ("search-kmer-len",
         po::value<int>(&opts->fs_kmer_len)->default_value(12,""),
         "length of k-mers (12)")
        ("search-kmer-mm",
         po::value<int>(&opts->fs_kmer_mm)->default_value(0,""),
         "allowed mismatches per k-mer (0)")
        ("search-kmer-norel",
         po::bool_switch(&opts->fs_kmer_norel),
         "don't score k-mer distance relative to target length")
        ("search-ignore-super",
         po::bool_switch(&opts->ignore_super),
         "ignore sequences containing query")
        ("search-copy-fields",
         po::value<string>(&opts->copy_fields)->default_value(""),
         "copy fields from result sequences to query sequence")
        ;
    od.add(cseq_comparator::get_options_description("search-"));

    adv.add(od);
}


void
search_filter::validate_vm(boost::program_options::variables_map& vm,
                           po::options_description& desc) {
    if (!opts) {
        throw std::logic_error("search options not parsed?!");
    }
    // if we're disabled, don't check further
    if (vm["search"].as<bool>() == false) {
        return;
    }
 
    // we need a DB to search
    if (vm.count("search-db") == 0) {  
        // default to alignment db if no search db configured
        if (vm.count("db") && !vm["db"].as<string>().empty()) {
           std::vector<string> cmd(2);
           cmd[0]="--search-db";
           cmd[1]=vm["db"].as<string>();
           po::store(po::command_line_parser(cmd).options(desc).run(), vm);
        } else {
          throw std::logic_error("need search-db to search");
        }
    }

    // search-port defaults to ptport if search-db==db
    if (vm["search-port"].defaulted() && 
        vm["db"].as<string>() == vm["search-db"].as<string>()) {
        std::vector<string> cmd(2);
        cmd[0]="--search-port";
        cmd[1]=vm["ptport"].as<string>();
        po::store(po::command_line_parser(cmd).options(desc).run(), vm);
    }

    opts->comparator = cseq_comparator::make_from_variables_map(vm, "search-");

    boost::split(opts->v_lca_fields, opts->lca_fields,
                 boost::is_any_of(":"));
    if (opts->v_lca_fields.back().empty()) {
        opts->v_lca_fields.pop_back();
    }
    boost::split(opts->v_copy_fields, opts->copy_fields,
                 boost::is_any_of(":"));
    if (opts->v_copy_fields.back().empty()) {
        opts->v_copy_fields.pop_back();
    }

}

} // namespace sina;

using namespace sina;

struct search_filter::priv_data {
    search *index;
    query_arb *arb;
    vector<cseq*> sequences;
};

search_filter::search_filter()
    :  data(new priv_data)
{
    data->arb = query_arb::getARBDB(opts->pt_database.c_str());

    if (opts->search_all) {
        logger->info("Caching Sequencences...");
        vector<string> names = data->arb->getSequenceNames();
        boost::progress_display p(data->arb->getSeqCount(), std::cerr);
        for (string& name: names) {
            data->sequences.push_back(&data->arb->getCseq(name));
            ++p;
        }
    } else {
        data->index = new query_pt(opts->pt_port.c_str(), opts->pt_database.c_str(),
                                   not opts->fs_no_fast,
                                   opts->fs_kmer_len,
                                   opts->fs_kmer_mm,
                                   opts->fs_kmer_norel);
    }
}

search_filter::search_filter(const search_filter& o)
    : data(o.data)
{
}

search_filter& search_filter::operator=(const search_filter& o) {
    data = o.data;
    return *this;
}

search_filter::~search_filter()
{
}

template<typename F>
struct dereference {
    explicit dereference(F f) : _f(f){}
    dereference() : _f(){}
    typedef typename F::result_type  result_type;

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

struct icontains {
    typedef bool result_type;
    const string bases;
    explicit icontains(const string& _bases) : bases(_bases) {}
    bool operator()(const cseq& c) {
        return boost::algorithm::icontains(c.getBases(), bases);
    }
};

struct iupac_contains {
    struct iupac_compare {
        typedef bool result_type;
        bool operator()(const aligned_base& a, const aligned_base& b) const {
            return a.comp(b);
        }
    };
    typedef bool result_type;
    vector<aligned_base> ref;
    explicit iupac_contains(cseq& c) : ref(c.getAlignedBases()) {}
    bool operator()(cseq& c) {
        return boost::algorithm::contains(c.getAlignedBases(), ref, iupac_compare());
    }
};

tray
search_filter::operator()(tray t) {
    cseq *c;
    if (t.aligned_sequence) {
        c = t.aligned_sequence;
    } else {
        t.log() << "search: no sequence?!;";
        return t;
    }

    if (c->size() < 20) {
        t.log() << "search:sequence too short (<20 bases);";
        return t;
    }

    t.search_result = new vector<cseq>();
    vector<cseq>& vc = *t.search_result; 

    string bases = c->getBases();
    if (opts->search_all) {
        for (cseq *r: data->sequences) {
            r->setScore(opts->comparator(*c, *r));
        }
        vector<cseq*>::iterator it = data->sequences.begin();
        vector<cseq*>::iterator middle = it + opts->max_result;
        const vector<cseq*>::iterator end = data->sequences.end();

        
        do {
            partial_sort(it, middle, end, dereference<std::greater<cseq> >());
            if (opts->ignore_super) {
                // sort sequences containing query to beginning
                // reset "start" to beginnin of non-identical section
                middle = it + opts->max_result;
                it = partition(it, middle, dereference<iupac_contains>(iupac_contains(*c)));
            }
            // repeat if we removed some just above
        } while (middle != end && it + opts->max_result > middle); 

        while (it != middle && (*it)->getScore() > opts->min_sim) {
            vc.push_back(**it);
            ++it;
        }
    } else {
        data->index->match(vc, *c, 1, opts->kmer_candidates, 0.3, 2.0, data->arb,
                           false, 50, 0, 0, 0,false);
        if (opts->ignore_super) {
            vector<cseq>::iterator it;
            it = partition(vc.begin(), vc.end(), iupac_contains(*c));
            vc.erase(it, vc.end());
        }

        for (cseq &r: vc) {
            r.setScore(opts->comparator(*c,r));
        }

        vector<cseq>::iterator 
            it=vc.begin(), 
            middle = vc.begin() + opts->max_result,
            end=vc.end();

        if (middle>end) middle=end;

        partial_sort(it, middle, end, std::greater<cseq>());

        while (it != middle && it->getScore() > opts->min_sim) ++it;
        vc.erase(it, vc.end());
    }

    stringstream result;
    map<string, vector<vector<string> > > group_names_map;
    for (cseq &r: vc) {
        data->arb->loadKey(r, "acc");
        data->arb->loadKey(r, "version");
        data->arb->loadKey(r, "start");
        data->arb->loadKey(r, "stop");

        for (string& s: opts->v_lca_fields) {
            data->arb->loadKey(r, s.c_str());
            string tax_path = r.get_attr<string>(s.c_str());
            if (tax_path == "Unclassified;") continue;
            vector<string> group_names;
            boost::split(group_names, tax_path, boost::is_any_of(";"));
            if (group_names.back().empty() || group_names.back()  == " ") {
                group_names.pop_back();
            }
            group_names_map[s].push_back(group_names);
        }
        result << r.get_attr<string>("acc")
               << "." << r.get_attr<string>("version")
               << "." << r.get_attr<string>("start")
               << "." << r.get_attr<string>("stop")
               << "~" << r.getScore()
               << " ";

        string acc = r.get_attr<string>("acc");
        for (string& s: opts->v_copy_fields) {
            data->arb->loadKey(r,s);
            string v = r.get_attr<string>(s);
            c->set_attr<string>(string("copy_")+acc+string("_")+s, v);
        }
    }
    c->set_attr<string>("nearest_slv", result.str());

    for (string& s: opts->v_lca_fields) {
        stringstream result;
        vector<vector<string> > group_names = group_names_map[s];
        for (vector<string>& vs: group_names) {
            reverse(vs.begin(), vs.end());
        }
        int outliers = vc.size() * (1 - opts->lca_quorum);
        while (outliers >= 0 && group_names.size() > 0) {
            vector<vector<string> >::iterator it = group_names.begin();
            if (it->size() == 0) {
                group_names.erase(it);
                outliers--;
                continue;
            }
            string name = it->back();
            ++it;
            for (; it != group_names.end(); ++it) {
                if (it->size() == 0 || it->back() != name) break;
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
