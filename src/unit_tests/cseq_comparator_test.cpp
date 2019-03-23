/*
Copyright (c) 2012-2013 Arne Boeckman
Copyright (c) 2012-2013 Elmar Pruesse

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

#define BOOST_TEST_MODULE CSEQ_COMPARATOR_module
#include <boost/test/unit_test.hpp>

#include "../cseq_comparator.h"
using namespace sina;

#include <string>
using std::string;

#define EQUAL BOOST_CHECK_EQUAL
#define CASE BOOST_AUTO_TEST_CASE
#define FCASE BOOST_FIXTURE_TEST_CASE

BOOST_AUTO_TEST_SUITE(CSEQ_COMPARATOR_test);


cseq c1 ("", "----AGCUAGCU-------");
cseq c2 ("", "--------AGCUAGCU---");
cseq c3 ("", "----AGCUAGCUAGCU---");
cseq c4 ("", "--------AGCU-------");
cseq c5 ("", "----AGAGAGAG-------");
cseq c6 ("", "--------AGAGAGAG---");
cseq c7 ("", "----agagAGAGagag---");
cseq c8 ("", "----AG-CU-A-G-CU---");
cseq c9 ("", "----MGCU-A-G-CU--");
cseq c10("", "----NNCUBADG-CU----");
cseq c11("", "---VVVCUAGCU-----");

CASE(cover_basic) {
    cseq_comparator comp(CMP_IUPAC_OPTIMISTIC,
                         CMP_DIST_NONE,
                         CMP_COVER_ABS, false);
    
    EQUAL(comp(c1,c1), 8);
    EQUAL(comp(c1,c2), 4);
    EQUAL(comp(c1,c3), 8);
    EQUAL(comp(c1,c4), 4);
    EQUAL(comp(c1,c5), 4);
    EQUAL(comp(c1,c6), 2);
    EQUAL(comp(c1,c7), 4);
    EQUAL(comp(c2,c2), 8);
    EQUAL(comp(c2,c3), 8);
    EQUAL(comp(c2,c4), 4);
    EQUAL(comp(c2,c5), 2);
    EQUAL(comp(c2,c6), 4);
    EQUAL(comp(c2,c7), 4);
    EQUAL(comp(c3,c3), 12);
    EQUAL(comp(c3,c4), 4);
    EQUAL(comp(c3,c5), 4);
    EQUAL(comp(c3,c6), 4);
    EQUAL(comp(c3,c7), 6);
    EQUAL(comp(c4,c4), 4);
    EQUAL(comp(c4,c5), 2);
    EQUAL(comp(c4,c6), 2);
    EQUAL(comp(c4,c7), 2);
    EQUAL(comp(c5,c5), 8);
    EQUAL(comp(c5,c6), 4);
    EQUAL(comp(c5,c7), 8);
    EQUAL(comp(c6,c6), 8);
    EQUAL(comp(c6,c7), 8);
    EQUAL(comp(c7,c7), 12);
}

CASE(cover_lowercase) {
  cseq_comparator comp(CMP_IUPAC_OPTIMISTIC,
                       CMP_DIST_NONE,
                       CMP_COVER_ABS, true);

  EQUAL(comp(c1,c7),2);
  EQUAL(comp(c2,c7),2);
  EQUAL(comp(c3,c7),2);
  EQUAL(comp(c4,c7),2);
  EQUAL(comp(c5,c7),4);
  EQUAL(comp(c6,c7),4);
  EQUAL(comp(c7,c7),4);

}

CASE(iupac_compare) {
    cseq_comparator comp(CMP_IUPAC_PESSIMISTIC,
                         CMP_DIST_NONE,
                         CMP_COVER_ABS, false);

    //pessimistic tests
    EQUAL(comp(c1,c9), 3);
    EQUAL(comp(c1,c2), 4);
    EQUAL(comp(c1,c10),2);
    EQUAL(comp(c1,c11),6);


    //optimistic tests
    cseq_comparator comp2(CMP_IUPAC_OPTIMISTIC,
                          CMP_DIST_NONE,
                          CMP_COVER_ABS, false);
    EQUAL(comp2(c1,c9), 4);
    EQUAL(comp2(c1,c2), 4);
    EQUAL(comp2(c1,c10),4);
    EQUAL(comp2(c1,c11),8);

}

CASE(cover_either) {

    // CMP_COVER_QUERY and CMP_COVER_TARGET
    // should be the same, just reversed parameters
    cseq_comparator comp(CMP_IUPAC_PESSIMISTIC,
                         CMP_DIST_NONE,
                         CMP_COVER_QUERY, false);

    cseq_comparator comp2(CMP_IUPAC_PESSIMISTIC,
                         CMP_DIST_NONE,
                         CMP_COVER_TARGET, false);

    EQUAL(comp(c1,c1), comp2(c1,c1));
    EQUAL(comp(c1,c2), comp2(c2,c1));
    EQUAL(comp(c1,c3), comp2(c3,c1));
    EQUAL(comp(c1,c4), comp2(c4,c1));
    EQUAL(comp(c1,c5), comp2(c5,c1));
    EQUAL(comp(c1,c6), comp2(c6,c1));
    EQUAL(comp(c1,c7), comp2(c7,c1));
    EQUAL(comp(c2,c2), comp2(c2,c2));
    EQUAL(comp(c2,c3), comp2(c3,c2));
    EQUAL(comp(c2,c4), comp2(c4,c2));
    EQUAL(comp(c2,c5), comp2(c5,c2));
    EQUAL(comp(c2,c6), comp2(c6,c2));
    EQUAL(comp(c2,c7), comp2(c7,c2));
    EQUAL(comp(c3,c3), comp2(c3,c3));
    EQUAL(comp(c3,c4), comp2(c4,c3));
    EQUAL(comp(c3,c5), comp2(c5,c3));
    EQUAL(comp(c3,c6), comp2(c6,c3));
    EQUAL(comp(c3,c7), comp2(c7,c3));
    EQUAL(comp(c4,c4), comp2(c4,c4));
    EQUAL(comp(c4,c5), comp2(c5,c4));
    EQUAL(comp(c4,c6), comp2(c6,c4));
    EQUAL(comp(c4,c7), comp2(c7,c4));
    EQUAL(comp(c5,c5), comp2(c5,c5));
    EQUAL(comp(c5,c6), comp2(c6,c5));
    EQUAL(comp(c5,c7), comp2(c7,c5));
    EQUAL(comp(c6,c6), comp2(c6,c6));
    EQUAL(comp(c6,c7), comp2(c7,c6));
    EQUAL(comp(c7,c7), comp2(c7,c7));


}

CASE(cover_overlap) {
    cseq_comparator comp(CMP_IUPAC_PESSIMISTIC,
                         CMP_DIST_NONE,
                         CMP_COVER_OVERLAP, false);

    //OVERLAP: relative to length of overlapping part
    EQUAL(comp(c1,c1), 1);
    EQUAL(comp(c1,c2), 1);
    EQUAL(comp(c1,c3), 1);
    EQUAL(comp(c1,c4), 1);
    EQUAL(comp(c1,c5), 4.0f/8);
    EQUAL(comp(c1,c6), 2.0f/4);
    EQUAL(comp(c1,c8), 2.0f/8); 
    EQUAL(comp(c2,c2), 1);
    EQUAL(comp(c2,c3), 1);
    EQUAL(comp(c2,c4), 1);
    EQUAL(comp(c2,c5), 2.0f/4);
    EQUAL(comp(c2,c6), 2.0f/4);
    EQUAL(comp(c2,c7), 2.0f/4);
    EQUAL(comp(c2,c11), 1);
    EQUAL(comp(c10,c11), 2.0f/8);

}

CASE(cover_all) {

    cseq_comparator comp(CMP_IUPAC_PESSIMISTIC,
                         CMP_DIST_NONE,
                         CMP_COVER_ALL, false);


    EQUAL(comp(c1,c2), 4/12.0f);
    EQUAL(comp(c1,c3), 8/12.0f);
    EQUAL(comp(c1,c4), 4/8.0f);
    EQUAL(comp(c1,c5), 4/8.0f);
    EQUAL(comp(c1,c6), 2/12.0f);
    EQUAL(comp(c1,c7), 4/12.0f);
    EQUAL(comp(c2,c2), 8/8.0f);
    EQUAL(comp(c2,c3), 8/12.0f);
    EQUAL(comp(c2,c4), 4/8.0f);
    EQUAL(comp(c2,c5), 2/12.0f);
    EQUAL(comp(c2,c6), 4/8.0f);
    EQUAL(comp(c2,c7), 4/12.0f);
    EQUAL(comp(c3,c3), 12/12.0f);
    EQUAL(comp(c3,c4), 4/12.0f);
    EQUAL(comp(c3,c5), 4/12.0f);
    EQUAL(comp(c3,c6), 4/12.0f);
    EQUAL(comp(c3,c7), 6/12.0f);
    EQUAL(comp(c4,c4), 4/4.0f);
    EQUAL(comp(c4,c5), 2/8.0f);
    EQUAL(comp(c4,c6), 2/8.0f);
    EQUAL(comp(c4,c7), 2/12.0f);
    EQUAL(comp(c5,c5), 8/8.0f);
    EQUAL(comp(c5,c6), 4/12.0f);
    EQUAL(comp(c5,c7), 8/12.0f);
    EQUAL(comp(c6,c6), 8/8.0f);

}

CASE(cover_min){
    cseq_comparator comp(CMP_IUPAC_PESSIMISTIC,
                          CMP_DIST_NONE,
                          CMP_COVER_MIN, false);


      EQUAL(comp(c1,c2), 4/8.0f);
      EQUAL(comp(c1,c3), 1);
      EQUAL(comp(c1,c4), 1);
      EQUAL(comp(c1,c5), 4/8.0f);
      EQUAL(comp(c1,c6), 2/8.0f);
      EQUAL(comp(c1,c7), 4/8.0f);
      EQUAL(comp(c2,c2), 1);
      EQUAL(comp(c2,c3), 1);

}

CASE(cover_max){
    cseq_comparator comp(CMP_IUPAC_PESSIMISTIC,
                          CMP_DIST_NONE,
                          CMP_COVER_MAX, false);



//    cseq c1 ("",0, "----AGCUAGCU-------");
//    cseq c3 ("",0, "----AGCUAGCUAGCU---");

      EQUAL(comp(c1,c2), 4/8.0f);
      EQUAL(comp(c1,c3), 8/12.0f);

}


BOOST_AUTO_TEST_SUITE_END(); // CSEQ_COMPARATOR_test

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
