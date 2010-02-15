#ifndef BITS_SIGN_TRAITS_H
#define BITS_SIGN_TRAITS_H

// Arranges for traits::is_signed<IntType> and traits::make_unsigned<IntType>
// to be the corresponding C++-0x or Boost.TypeTraits templates, if available.

#if __cplusplus >= 201001L  // Actual C++-0x value is yet to be determined.
  #include <type_traits>
  namespace traits = std;
#else
  #include <boost/version.hpp>

  #if BOOST_VERSION >= 104000  // Check this: 1.39.1 bad; 1.40 good?
    #include <boost/type_traits/is_signed.hpp>
    #include <boost/type_traits/make_unsigned.hpp>
    namespace traits = boost;
  #else
    #include <boost/type_traits/is_signed.hpp>

    // Debian ships a version of Boost that predates boost::make_unsigned<>.
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
