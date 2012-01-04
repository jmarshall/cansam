/*  bits/sign_traits.h -- Provides some missing type traits.

    Copyright (C) 2010-2012 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
    nor the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND ITS CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */

#ifndef BITS_SIGN_TRAITS_H
#define BITS_SIGN_TRAITS_H

// Arranges for traits::is_signed<IntType> and traits::make_unsigned<IntType>
// to be the corresponding C++11 or Boost.TypeTraits templates, if available.

#if __cplusplus >= 201103L
  #include <type_traits>
  namespace traits = std;
#else
  #include <boost/version.hpp>

  #if BOOST_VERSION >= 103500
    #include <boost/type_traits/is_signed.hpp>
    #include <boost/type_traits/make_unsigned.hpp>
    namespace traits = boost;
  #else
    #include <boost/type_traits/is_signed.hpp>

    // Debian 5.0 (lenny) ships a version of Boost without make_unsigned<>.
    // This defines it for integer types, which is all we need.
    namespace traits {
      using boost::false_type;
      using boost::true_type;
      using boost::is_signed;

      template <typename IntType>
      struct make_unsigned { typedef IntType type; };

      template<> struct make_unsigned<char>  { typedef unsigned char  type; };
      template<> struct make_unsigned<signed char>
					     { typedef unsigned char  type; };
      template<> struct make_unsigned<short> { typedef unsigned short type; };
      template<> struct make_unsigned<int>   { typedef unsigned int   type; };
      template<> struct make_unsigned<long>  { typedef unsigned long  type; };
      template<> struct make_unsigned<long long>
					 { typedef unsigned long long type; };
    }
  #endif
#endif

#endif
