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

#include <vector>
#include <algorithm> 
#include <fstream>

#define BOOST_TEST_MODULE cseq_module
#include <boost/test/unit_test.hpp>

#include <boost/algorithm/string.hpp>

#include "../cseq.h"
#include "../aligned_base.h"
using sina::cseq;
using sina::aligned_base;

#include <string>
using std::string;

#include <iostream>

#define EQUAL BOOST_CHECK_EQUAL
#define EQUAL_COLLECTIONS BOOST_CHECK_EQUAL_COLLECTIONS
#define CASE BOOST_AUTO_TEST_CASE
#define FIXTURE_CASE BOOST_FIXTURE_TEST_CASE

BOOST_AUTO_TEST_SUITE(cseq_test);

const string rna = "AGCURYKMSWBDHVN";
const string rna_aligned = "--A-G---CUR-YKM-S---WBD-HVN---";
const string rna_aligned_complemented = "--T-C---GAR-YKM-S---WBD-HVN---";
const string rna_aligned_dots = "..A-G---CUR-YKM-S---WBD-HVN...";
const float score = 15.f;


/*

  TODO:
  - dot vs gap

*/

#define test_empty(cs)                      \
    EQUAL(cs.size(),0);                     \
    EQUAL(cs.getWidth(),0);                 \
    EQUAL(cs.getBases(), string());         \
    EQUAL(cs.getAlignedNoDots(), string()); \
    EQUAL(cs.end() - c.begin(), 0);         \
    EQUAL(cs.rend() - c.rbegin(), 0);       \
    EQUAL(cs.getName(), string());          \
    EQUAL(cs.getScore(), 0.f); 


#define test_data(cs, name, score, _aligned)                    \
    {                                                           \
        string aligned = (_aligned);                            \
        string unaligned = aligned;                             \
        boost::erase_all(unaligned, "-");                       \
        EQUAL(cs.size(), unaligned.size());                     \
        EQUAL(cs.getWidth(), aligned.size());                   \
        EQUAL(cs.getBases(), unaligned);                        \
        EQUAL(cs.getAlignedNoDots(), aligned);                  \
        EQUAL(cs.end() - cs.begin(), unaligned.size());         \
        EQUAL(cs.rend() - cs.rbegin(), unaligned.size());       \
        EQUAL(cs.getName(), name);                              \
        EQUAL(cs.getScore(), score);                            \
    }

CASE(test_constructor_empty) {
    cseq c;

    test_empty(c);
}

CASE(test_constructur_normal) {
    const string name("thename");
    const float score=15.f;
    cseq c(name.c_str(), score, rna.c_str());

    test_data(c, name, score, rna);
}

CASE(test_append,
     *boost::unit_test::expected_failures(1)) {
    cseq c;

    c.append(rna);
    test_data(c, "", 0, rna);
    c.append("");
    test_data(c, "", 0, rna);
    c.append(rna);
    test_data(c, "", 0, rna+rna);
    c.clearSequence();
    test_empty(c);
    c.append(rna_aligned);
    test_data(c, "", 0, rna_aligned);
    c.append(rna);
    test_data(c, "", 0, rna_aligned+rna);
    c.append(rna_aligned);
    test_data(c, "", 0, rna_aligned+rna+rna_aligned);


    aligned_base b(0,'A');
    std::cerr << "==================" << std::endl
              << "Triggering 1 message: "
              << "\"$ cseq::append(): wrong order! A(0<75)\""
              << std::endl;
    c.append(b);
    std::cerr << "Triggering 1 error: "
              << "c.getWidth() == 75 and aligned.size() == 76"
              << std::endl;
    test_data(c, "", 0, rna_aligned+rna+rna_aligned + "A");
    std::cerr << "==================" << std::endl << std::endl;

}


CASE(test_setWidth) {
    cseq c;
    string twentygaps="--------------------";

    c.setWidth(20);
    test_data(c, "", 0, twentygaps);
    c.setWidth(40);
    test_data(c, "", 0, twentygaps + twentygaps);
    c.setWidth(20);
    test_data(c, "", 0, twentygaps);

    c.setWidth(0);
    test_data(c, "", 0, "");
}


CASE(test_reverse) {
    string name("testtt");
    cseq c(name.c_str(), score, rna_aligned.c_str());
    string reversed = rna_aligned;
    std::reverse(reversed.begin(), reversed.end());

    c.reverse();
    test_data(c, name, score, reversed);
    c.reverse();
    test_data(c, name, score, rna_aligned);
}

CASE(test_lowercase) {
    string rna_aligned_lower = boost::to_lower_copy(rna_aligned);
    string rna_lower = boost::to_lower_copy(rna);
    cseq c("", 0, rna_lower.c_str());
    EQUAL(c.getAlignedNoDots(), rna_lower);
}

CASE(test_uppercase){

    string rna_aligned_lower = boost::to_lower_copy(rna_aligned);
    string rna_lower = boost::to_lower_copy(rna);
    cseq c("", 0, rna_lower.c_str());
    c.upperCaseAll();
    EQUAL(c.getAlignedNoDots(), rna);
    c.append(rna_aligned_lower);
    c.upperCaseAll();
    EQUAL(c.getAlignedNoDots(), rna + rna_aligned);
}

CASE(test_complement)
{
    //reference for complements:
    //http://www.animalgenome.org/edu/gene/genetic-code.html
    cseq c("",0,rna.c_str());
    c.complement();
    EQUAL(c.size(),rna.length());
    EQUAL(c.getBases(),"UCGAYRMKSWVHDBN");
}

class compression_data {
public:
    const cseq c_alig; 
    const cseq c_unalig;
    static const unsigned char c_alig_compr[];
    static const unsigned char c_unalig_compr[];

    compression_data() 
        : c_alig("", 0, rna_aligned.c_str()),
          c_unalig("", 0, rna.c_str())
    {
    }
};

const unsigned char
compression_data::c_alig_compr[] = {
    0x23, 0x00, 0x80, 0x00, 0x78, 0xda, 0x73, 0x74, 
    0x77, 0x0e, 0x0d, 0x8a, 0xf4, 0xf6, 0x0d, 0x0e,
    0x77, 0x72, 0xf1, 0x08, 0xf3, 0xd3, 0x63, 0x62, 
    0x62, 0x61, 0x64, 0x64, 0x02, 0x22, 0x08, 0xc5,
    0xc2, 0x40, 0x63, 0x00, 0x00, 0x40, 0xd5, 0x04,
    0xcc
        /* (used this to generate verification data)
        std::ofstream of("/tmp/thefile1");
        of << std::string((char*)&data.front(), data.size());
        */
};

const unsigned char
compression_data::c_unalig_compr[] = {
    0x23, 0x00, 0x80, 0x00, 0x78, 0xda, 0x73, 0x74, 
    0x77, 0x0e, 0x0d, 0x8a, 0xf4, 0xf6, 0x0d, 0x0e,
    0x77, 0x72, 0xf1, 0x08, 0xf3, 0xd3, 0x63, 0x60, 
    0x44, 0x05, 0x0c, 0x34, 0x06, 0x00, 0x3a, 0xad,
    0x04, 0xbd
};

FIXTURE_CASE(test_compress_unaligned, compression_data) {
    std::vector<unsigned char> data;
    cseq c = c_unalig;
    c.compressAligned(data);
    EQUAL_COLLECTIONS(data.begin(), data.end(), c_unalig_compr, 
                      c_unalig_compr + sizeof(c_unalig_compr));
}

FIXTURE_CASE(test_compress_aligned, compression_data) {
    std::vector<unsigned char> data;
    cseq c = c_alig;
    c.compressAligned(data);
    EQUAL_COLLECTIONS(data.begin(), data.end(), c_alig_compr, 
                      c_alig_compr + sizeof(c_alig_compr));
}

FIXTURE_CASE(test_decompress_unaligned, compression_data) {
    std::vector<unsigned char> data;
    cseq c;
    c.assignFromCompressed(c_unalig_compr, sizeof(c_unalig_compr));
    test_data(c, "", 0, rna.c_str());
}

FIXTURE_CASE(test_decompress_aligned, compression_data) {
    std::vector<unsigned char> data;
    cseq c;
    c.assignFromCompressed(c_alig_compr, sizeof(c_alig_compr));
    test_data(c, "", 0, rna_aligned.c_str());
}


CASE(test_dots){
    cseq c("",0,rna_aligned.c_str());
    EQUAL(c.getAligned(false),rna_aligned_dots);
}


BOOST_AUTO_TEST_SUITE_END(); // cseq_test

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

