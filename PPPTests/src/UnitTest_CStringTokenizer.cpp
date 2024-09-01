#include <PPPLib/MFC/CStringTokenizer.h>
#include "catch.hpp"
#include <cstring>

namespace novac
{
TEST_CASE("NextToken", "[CStringTokenizer]")
{
    SECTION("Single space as separators")
    {
        CStringTokenizer sut{ "the brown fox jumps over the lazy dog", " " };

        REQUIRE(0 == strcmp("the", sut.NextToken()));
        REQUIRE(0 == strcmp("brown", sut.NextToken()));
        REQUIRE(0 == strcmp("fox", sut.NextToken()));
        REQUIRE(0 == strcmp("jumps", sut.NextToken()));
        REQUIRE(0 == strcmp("over", sut.NextToken()));
        REQUIRE(0 == strcmp("the", sut.NextToken()));
        REQUIRE(0 == strcmp("lazy", sut.NextToken()));
        REQUIRE(0 == strcmp("dog", sut.NextToken()));
        REQUIRE(nullptr == sut.NextToken());
    }
}
}
