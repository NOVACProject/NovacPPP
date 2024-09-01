#include "catch.hpp"
#include <PPPLib/MFC/CFtpUtils.h>

namespace novac
{
TEST_CASE("SplitPathIntoServerAndDirectory behaves as expected", "[CFtpUtils]")
{
    SECTION("Splits server and directory")
    {
        CString path{ "ftp://127.0.0.1/some/directory/path/" };

        std::string server, directory;
        CFtpUtils utils;
        utils.SplitPathIntoServerAndDirectory(path, server, directory);

        REQUIRE(server == "127.0.0.1");
        REQUIRE(directory == "some/directory/path/");
    }
}

TEST_CASE("IsFilePermissions behaves as expected", "[CFtpUtils]")
{
    CFtpUtils sut;

    SECTION("Returns false for typical file names and dates")
    {
        REQUIRE_FALSE(sut.IsFilePermissions("2009"));
        REQUIRE_FALSE(sut.IsFilePermissions("wind"));
        REQUIRE_FALSE(sut.IsFilePermissions("drewcd"));
        REQUIRE_FALSE(sut.IsFilePermissions("dddddddddd"));
    }

    SECTION("Returns true for typical file permission strings")
    {
        REQUIRE(sut.IsFilePermissions("drwx------"));
        REQUIRE(sut.IsFilePermissions("drwxr-xr-x"));
        REQUIRE(sut.IsFilePermissions("-rwxr-xr-x"));
        REQUIRE(sut.IsFilePermissions("-rw-r-----"));
    }
}

TEST_CASE("ReadFtpDirectoryListing behaves as expected", "[CFtpUtils]")
{
    CFtpUtils sut;
    CFileInfo result;

    SECTION("No data, return false")
    {
        REQUIRE_FALSE(sut.ReadFtpDirectoryListing(" ", result));
    }

    SECTION("simple ftp listing line, return true and sets file info")
    {
        REQUIRE(sut.ReadFtpDirectoryListing("drwx------ 2    ftp    ftp    4096   Jun   12   2009    2007.08.29", result));
        REQUIRE(result.isDirectory == true);
        REQUIRE(result.fileName == "2007.08.29");
    }
}

}