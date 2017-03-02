/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.

SINA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _QUERY_PT_H_
#define _QUERY_PT_H_

#include <string>
#include <vector>
#include <exception>

namespace sina {

class cseq;
class query_arb;

class query_pt {
  query_pt& operator=(const query_pt&);
  query_pt(const query_pt&);

  void init();
  void exit();
  void restart();
 public:
  query_pt(const char* portname, const char* dbname);
  ~query_pt();

  double match(std::vector<cseq> &family, const cseq& query,
	       int min_match, int max_match, float min_score, float max_score,
	       query_arb *arb, bool noid, int minlen,
	       int num_full, int minlen_full, int range_cover, bool leave_query_out);

  double match(std::vector<cseq> &family, const cseq& sequence,
	       int min_match, int max_match, float min_score) {
    return match(family, sequence, min_match, max_match, min_score, 2.0,
		 0, false, 0, 0, 0, 0, false);
  };

  int turn_check(const cseq& query, bool all);

  void set_find_type_fast(bool fast);
  void set_probe_len(int len);
  void set_mismatches(int len);
  void set_sort_type(bool absolute);
  void set_range(int startpos, int stoppos);
  void unset_range();

  class exception : public std::exception {
      std::string message;
  public: 
      exception(std::string _message) throw();
      ~exception() throw();
      virtual const char* what() const throw();
  };

 private:
  struct priv_data;
  priv_data &data;
  bool find_type_fast;

};

} // namespace sina

#endif // _QUERY_PT_H
