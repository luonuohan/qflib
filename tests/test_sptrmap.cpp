#include <catch2/catch_test_macros.hpp>
#include <qflib/sptrmap.hpp>
#include <string>
#include <memory>

TEST_CASE("SPtrMap basic operations", "[sptrmap]")
{
    qf::SPtrMap<int> map;

    SECTION("empty map contains nothing") {
        REQUIRE_FALSE(map.contains("foo"));
        REQUIRE(map.list().empty());
    }

    SECTION("set and get") {
        auto sp = std::make_shared<int>(42);
        auto [name, ver] = map.set("mykey", sp);
        REQUIRE(name == "MYKEY");   // uppercased
        REQUIRE(ver == 1);
        REQUIRE(map.contains("mykey"));
        REQUIRE(*map.get("mykey") == 42);
    }

    SECTION("name is trimmed and uppercased") {
        auto sp = std::make_shared<int>(7);
        auto [name, ver] = map.set("  hello  ", sp);
        REQUIRE(name == "HELLO");
        REQUIRE(map.contains("hello"));
        REQUIRE(map.contains("  HELLO  "));
    }

    SECTION("version increments on update") {
        auto sp1 = std::make_shared<int>(1);
        auto sp2 = std::make_shared<int>(2);
        auto [n1, v1] = map.set("key", sp1);
        REQUIRE(v1 == 1);
        auto [n2, v2] = map.set("key", sp2);
        REQUIRE(v2 == 2);
        REQUIRE(*map.get("key") == 2);
        REQUIRE(map.version("key") == 2);
    }

    SECTION("list returns all keys") {
        map.set("alpha", std::make_shared<int>(1));
        map.set("beta", std::make_shared<int>(2));
        auto keys = map.list();
        REQUIRE(keys.size() == 2);
        // keys are uppercased and sorted (std::map)
        REQUIRE(keys[0] == "ALPHA");
        REQUIRE(keys[1] == "BETA");
    }

    SECTION("clear empties the map") {
        map.set("x", std::make_shared<int>(99));
        REQUIRE(map.contains("x"));
        map.clear();
        REQUIRE_FALSE(map.contains("x"));
        REQUIRE(map.list().empty());
    }

    SECTION("get returns nullptr for missing key") {
        auto ptr = map.get("nonexistent");
        REQUIRE(ptr == nullptr);
    }

    SECTION("names with internal blanks throw") {
        auto sp = std::make_shared<int>(1);
        REQUIRE_THROWS_AS(map.set("bad name", sp), qf::Exception);
    }

    SECTION("empty names throw") {
        auto sp = std::make_shared<int>(1);
        REQUIRE_THROWS_AS(map.set("", sp), qf::Exception);
        REQUIRE_THROWS_AS(map.set("   ", sp), qf::Exception);
    }
}

TEST_CASE("trim utility function", "[sptrmap]")
{
    REQUIRE(qf::trim("  hello  ") == "hello");
    REQUIRE(qf::trim("hello") == "hello");
    REQUIRE(qf::trim("   ") == "");
    REQUIRE(qf::trim("") == "");
    REQUIRE(qf::trim("\t\n hello \n\t") == "hello");
}
