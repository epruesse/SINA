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

#include "config.h"
#include "kmer_search.h"
#include "query_arb.h"

#include "query_pt.h"
#include "timer.h"

#include <vector>
#include <string>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#ifdef HAVE_TBB
#include <tbb/task_scheduler_init.h>
#endif

using namespace sina;
using std::cerr;
using std::endl;


int main(int argc, const char** argv) {
  timer t;

  if (argc<2) {
      std::cerr << "need arb db as argument" << std::endl;
      exit(1);
  }

#ifdef HAVE_TBB
  int n_threads = tbb::task_scheduler_init::default_num_threads();
  tbb::task_scheduler_init init(n_threads);
#endif

  t.start();
  query_arb* arbdb = query_arb::getARBDB(argv[1]);
  t.stop();
  /*
  query_pt pt("localhost:4040", argv[1]);
  t.stop();
  */

  t.stop();
  cerr << "Preloading sequence cache..." << endl;
  std::vector<std::string> seqNames = arbdb->getSequenceNames();
  for (int i = 0; i < seqNames.size(); i++) {
    //cseq target = arbdb->getCseq(seqNames[i]);
  }
  t.stop();
  kmer_search* search_index=kmer_search::get_kmer_search(argv[1], 10);

  t.stop();
  std::vector<cseq> family;
  for (int i = 0; i < 1000; i++) {
    cseq target = arbdb->getCseq(seqNames[i]);
    search_index->find(target, family, 10);
    float max_score = family[0].getScore();
    std::vector<cseq>::iterator self;
    self = std::find_if(family.begin(), family.end(),
			[&](const cseq &c) {
			  return c.getName() == target.getName();}
			);
    if (self == family.end()) {
      std::cerr << seqNames[i] << ": " << std::endl;
      for (int i=0; i<family.size(); i++) {
	std::cerr << family[i].getNameScore() << std::endl;
      }
      std::cerr << std::endl;
    }
  }
  t.stop();
  
  /*
  t.stop();
  for (int i = 0; i < 100; i++) {
    cseq target = arbdb->getCseq(seqNames[i]);
    std::vector<cseq> family;
    pt.match(family, target, 40, 40, 0);
  }
  */
  std::cerr << t << std::endl;
  

  return 0;
}
