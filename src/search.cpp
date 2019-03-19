#include "search.h"
#include <ostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::istarts_with;
using boost::algorithm::iequals;

namespace sina {

void validate(boost::any& v,
              const std::vector<std::string>& values,
              ENGINE_TYPE* /*tt*/, int /*unused*/) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
    const std::string& s = validators::get_single_string(values);
    if (iequals(s, "pt-server")) {
        v = ENGINE_ARB_PT;
    } else if (iequals(s, "internal")) {
        v = ENGINE_SINA_KMER;
    } else {
        throw po::invalid_option_value(s);
    }
}

std::ostream& operator<<(std::ostream& out, const ENGINE_TYPE& t) {
    switch(t) {
    case ENGINE_ARB_PT: out << "pt-server"; break;
    case ENGINE_SINA_KMER: out << "internal"; break;
    default: out << "[UNKNOWN!]";
    }
    return out;
}
  
} // namespace sina
