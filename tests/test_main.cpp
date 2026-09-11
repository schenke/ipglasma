// Generates doctest's own main(). Every other test_*.cpp in this directory
// just includes doctest.h without this macro defined.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
