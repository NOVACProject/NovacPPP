#include <PPPLib/MFC/CList.h>
#include <PPPLib/MFC/CString.h>
#include "catch.hpp"

namespace novac
{
TEST_CASE("POSITION")
{
    SECTION("Construction from nullptr")
    {
        POSITION<int> p2 = nullptr;

        REQUIRE_FALSE(p2.HasNext());
    }

    SECTION("Assignment to nullptr, does not alter list")
    {
        CList<int> originalList;
        originalList.AddTail(1);
        originalList.AddTail(2);

        POSITION<int> sut = originalList.GetHeadPosition();

        REQUIRE(originalList.GetSize() == 2);

        sut = nullptr;

        REQUIRE(originalList.GetSize() == 2);
    }
}

TEST_CASE("REVERSE_POSITION")
{
    SECTION("Construction from nullptr")
    {
        REVERSE_POSITION<int> p2 = nullptr;

        REQUIRE_FALSE(p2.HasPrevious());
    }
}


TEST_CASE("CList Construction behaves as expected", "[CList]")
{
    SECTION("Default constructor")
    {
        CList<int> defaultContructedCList;
        REQUIRE(defaultContructedCList.GetSize() == 0);
        REQUIRE(defaultContructedCList.GetCount() == 0);
    }
}

TEST_CASE("CList append and remove behaves as expected", "[CList]")
{
    SECTION("AddTail - increases size")
    {
        CList<int> sut;
        REQUIRE(0 == sut.GetSize());
        sut.AddTail(2);
        REQUIRE(1 == sut.GetSize());
    }

    SECTION("AddTail - appends new element")
    {
        CList<int> sut;
        sut.AddTail(2);

        auto p = sut.GetTailPosition();
        REQUIRE(2 == sut.GetAt(p));
    }

    SECTION("AddHead - increases size")
    {
        CList<int> sut;
        REQUIRE(0 == sut.GetSize());
        sut.AddHead(2);
        REQUIRE(1 == sut.GetSize());
    }

    SECTION("AddHead - appends new element")
    {
        CList<int> sut;
        sut.AddHead(2);

        auto p = sut.GetHeadPosition();
        REQUIRE(2 == sut.GetAt(p));
    }

    SECTION("RemoveTail - removes last element")
    {
        CList<int> sut;
        sut.AddTail(2);
        sut.AddTail(4);
        sut.AddTail(8);
        REQUIRE(3 == sut.GetSize());

        sut.RemoveTail();

        REQUIRE(2 == sut.GetSize());

        auto p = sut.GetTailPosition();
        REQUIRE(4 == sut.GetAt(p));
    }

    SECTION("InsertBefore - inserts one element before")
    {
        CList<int> sut;
        sut.AddTail(2);

        auto p = sut.GetHeadPosition();

        sut.InsertBefore(p, 3);

        REQUIRE(2 == sut.GetSize());

        auto newHead = sut.GetHeadPosition();
        REQUIRE(3 == sut.GetAt(newHead));
    }
}

TEST_CASE("CList iteration behaves as expected", "[CList]")
{
    SECTION("Empty list")
    {
        CList<CString> sut;
        int nItems = 0;

        auto p = sut.GetHeadPosition();
        while (p.HasNext())
        {
            ++nItems;

            sut.GetNext(p);
        }

        REQUIRE(0 == nItems);
    }

    SECTION("Forward in 5 item list")
    {
        CList<CString> sut;
        sut.AddTail(CString("anders"));
        sut.AddTail(CString("berit"));
        sut.AddTail(CString("calle"));
        sut.AddTail(CString("dennis"));
        sut.AddTail(CString("eva"));

        int nItems = 0;

        auto p = sut.GetHeadPosition();
        while (p.HasNext())
        {
            ++nItems;

            auto element = sut.GetNext(p);
        }

        REQUIRE(5 == nItems);
    }

    SECTION("Forward in 5 item list - Comparison with nullptr")
    {
        CList<int> sut;
        sut.AddTail(1);
        sut.AddTail(2);
        sut.AddTail(3);
        sut.AddTail(4);
        sut.AddTail(5);

        int nItems = 0;

        auto p = sut.GetHeadPosition();
        while (p != nullptr)
        {
            ++nItems;

            sut.GetNext(p);
        }

        REQUIRE(5 == nItems);
    }

    SECTION("Backwards in 5 item list")
    {
        CList<int> sut;
        sut.AddTail(1);
        sut.AddTail(2);
        sut.AddTail(3);
        sut.AddTail(4);
        sut.AddTail(5);

        int nItems = 0;

        auto p = sut.GetTailPosition();
        while (p.HasPrevious())
        {
            ++nItems;

            sut.GetPrev(p);
        }

        REQUIRE(5 == nItems);
    }

    SECTION("Backwards in 5 item list - Comparison with nullptr")
    {
        CList<int> sut;
        sut.AddTail(1);
        sut.AddTail(2);
        sut.AddTail(3);
        sut.AddTail(4);
        sut.AddTail(5);

        int nItems = 0;

        auto p = sut.GetTailPosition();
        while (p != nullptr)
        {
            ++nItems;

            sut.GetPrev(p);
        }

        REQUIRE(5 == nItems);
    }
}
}