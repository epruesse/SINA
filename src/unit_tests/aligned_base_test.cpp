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


#ifndef _aligned_base_TEST_HPP_
#define _aligned_base_TEST_HPP_

#define BOOST_TEST_MODULE aligned_base_module
#include <boost/test/unit_test.hpp>

#include "../aligned_base.h"
using namespace sina;
#define EQUAL BOOST_CHECK_EQUAL
#define TRUE(x) BOOST_CHECK_EQUAL(x,true)
#define FALSE(x) BOOST_CHECK_EQUAL(x,false)
#define CASE BOOST_AUTO_TEST_CASE
#define FCASE BOOST_FIXTURE_TEST_CASE

BOOST_AUTO_TEST_SUITE(aligned_base_test);

unsigned char A('A');
unsigned char G('G');
unsigned char C('C');
unsigned char U('U');
unsigned char T('T');
CASE(ctor_test){

    base_iupac baseA(A);
    base_iupac baseT(U);

    TRUE(baseA.has_A());
    FALSE(baseA.has_C());
    FALSE(baseA.has_G());
    FALSE(baseA.has_TU());

    FALSE(baseT.has_A());
    FALSE(baseT.has_C());
    FALSE(baseT.has_G());
    TRUE(baseT.has_TU());

    EQUAL((char)baseA,A);
    EQUAL((char)baseT,U);



}

CASE(ctor_illegal_param_test){
    //what happens if we use an illegal char

    unsigned char illegal('!');
    BOOST_CHECK_THROW(base_iupac baseIllegal(illegal),sina::aligned_base::bad_character_exception);


}

CASE(ctor_empty_param_test){
    base_iupac b;
    FALSE(b.has_A());
    FALSE(b.has_C());
    FALSE(b.has_G());
    FALSE(b.has_TU());
}

CASE(complement_test){
    base_iupac baseA(A);

    baseA.complement();
    EQUAL((char)baseA,U);

    base_iupac baseC(C);

    baseC.complement();
    EQUAL((char)baseC,G);

}

CASE(setLower_test){
    base_iupac baseA(A);
    base_iupac baseC(C);

    baseA.setLowerCase();
    baseC.setLowerCase();

    EQUAL((char)baseA,'a');
    EQUAL((char)baseC,'c');
}


CASE(setUpper_test){
    base_iupac baset('t');
    base_iupac baseg('g');

    baset.setUpperCase();
    baseg.setUpperCase();

    EQUAL((char)baset,'U');
    EQUAL((char)baseg,'G');
}

CASE(isLower_test){
    base_iupac baset('t');
    base_iupac baseG('G');

    TRUE(baset.isLowerCase());
    FALSE(baseG.isLowerCase());
}

CASE(comp_test){
    base_iupac baseT('T');
    base_iupac baseU('U');
    base_iupac baseu('u');
    base_iupac baset('t');
    base_iupac baseG('G');
    base_iupac baseC('C');

    TRUE(baseT.comp(baseU));
    TRUE(baseT.comp(baseu));
    TRUE(baseu.comp(baseT));
    TRUE(baseT.comp(baseT));
    TRUE(baset.comp(baseT));
    TRUE(baseG.comp(baseG));
    TRUE(baseC.comp(baseC));
    FALSE(baseT.comp(baseC));
    FALSE(baseG.comp(baseC));
    FALSE(baset.comp(baseC));

}

CASE(comp_pessimistic_test){
    base_iupac base('R');
    base_iupac base2('T');

    //R should not be equal to R because R is ambiguous
    FALSE(base.comp_pessimistic(base));
    base = 'A';
    TRUE(base.comp_pessimistic(base));
    base = 'G';
    TRUE(base.comp_pessimistic(base));
    FALSE(base.comp_pessimistic(base2));

}


CASE(pair_test){
    base_iupac baseT(T);
    base_iupac baseA(A);
    EQUAL(baseT.pair(baseA),1.0f);

}



CASE(is_ambig_test){
    base_iupac baseT('T');
    FALSE(baseT.is_ambig());

    base_iupac baseM('m');
    TRUE(baseM.is_ambig());

    base_iupac baseV('V');
    TRUE(baseV.is_ambig());

    base_iupac baseG('G');
    FALSE(baseG.is_ambig());
}

CASE(ambig_order_test){
    base_iupac baseT('T');
    EQUAL(baseT.ambig_order(),1);

    base_iupac baseM('M');
    EQUAL(baseM.ambig_order(),2);

    base_iupac baseD('D');
    EQUAL(baseD.ambig_order(),3);
}

CASE(cast_to_char_test){


    //iterate over small letter alphabet
    for(char i = 97; i <=122; i++ )
    {
        try
        {
            //skip t, it is identical to u :)
            if('t' == i)
            {
                i++;
            }
            base_iupac base(i);//it may either throw an exception
            EQUAL((char)base, i);//or it should convert correctly
        }
        catch(sina::aligned_base::bad_character_exception)
        {
            EQUAL(true,true);
        }

    }


}


BOOST_AUTO_TEST_SUITE_END(); // XXX_test

#endif /* _XXX_TEST_HPP_ */

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
