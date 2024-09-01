#include <PPPLib/MFC/CArray.h>
#include "catch.hpp"

namespace novac
{
TEST_CASE("CArray Construction behaves as expected", "[CArray]")
{
    SECTION("Default constructor")
    {
        CArray<int> defaultContructedCList;
        REQUIRE(defaultContructedCList.GetSize() == 0);
    }

    SECTION("Copy from source")
    {
        CArray<int> original;
        original.SetAtGrow(0, 5);
        original.SetAtGrow(1, 6);
        original.SetAtGrow(2, 7);
        original.SetAtGrow(3, 8);

        CArray<int> copy;
        copy.Copy(original);

        REQUIRE(copy.GetSize() == 4);

        // also make sure that the copy doesn't just reference the same memory...
        original.SetAtGrow(4, 9);

        REQUIRE(original.GetSize() == 5);
        REQUIRE(copy.GetSize() == 4);
    }
}

TEST_CASE("CArray retrieving elements behaves as expected", "[CArray]")
{
    SECTION("Bracket operator, get and set")
    {
        CArray<int> sut;
        sut.SetSize(5);

        REQUIRE(sut[3] == 0);

        sut[3] = 65;

        REQUIRE(sut[3] == 65);
    }

    SECTION("GetAt")
    {
        CArray<int> sut;
        sut.SetSize(5);
        sut[4] = 8765;

        REQUIRE(8765 == sut.GetAt(4));
    }

    SECTION("SetAt")
    {
        CArray<int> sut;
        sut.SetSize(5);
        sut.SetAt(4, 8765);

        REQUIRE(8765 == sut.GetAt(4));
    }

    SECTION("SetAtGrow")
    {
        CArray<int> sut;
        sut.SetSize(5);

        sut.SetAtGrow(76, 8765);

        REQUIRE(77 == sut.GetSize());
        REQUIRE(8765 == sut.GetAt(76));
    }
}

TEST_CASE("CArray operations behave as expected", "[CArray]")
{
    SECTION("SetSize")
    {
        CArray<int> sut;
        REQUIRE(sut.GetSize() == 0);

        sut.SetSize(1927);
        REQUIRE(sut.GetSize() == 1927);
    }

    SECTION("RemoveAll")
    {
        CArray<int> sut;
        sut.SetSize(5);

        REQUIRE(sut.GetSize() == 5);

        sut.RemoveAll();

        REQUIRE(sut.GetSize() == 0);
    }

    SECTION("RemoveAt - 1")
    {
        CArray<int> sut;
        sut.SetSize(5);
        sut[0] = 0;
        sut[1] = 1;
        sut[2] = 2;
        sut[3] = 3;
        sut[4] = 4;

        sut.RemoveAt(0);

        REQUIRE(sut.GetSize() == 4);
        REQUIRE(sut.GetAt(0) == 1);
        REQUIRE(sut.GetAt(1) == 2);
        REQUIRE(sut.GetAt(2) == 3);
        REQUIRE(sut.GetAt(3) == 4);
    }

    SECTION("RemoveAt - 2")
    {
        CArray<int> sut;
        sut.SetSize(5);
        sut[0] = 0;
        sut[1] = 1;
        sut[2] = 2;
        sut[3] = 3;
        sut[4] = 4;

        sut.RemoveAt(1, 3);

        REQUIRE(sut.GetSize() == 2);
        REQUIRE(sut.GetAt(0) == 0);
        REQUIRE(sut.GetAt(1) == 4);
    }
}
}