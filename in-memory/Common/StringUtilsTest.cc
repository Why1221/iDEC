#include "Exception.h"
#include "StringUtils.hpp"

#include <iostream>
using namespace std;
using namespace npp;
using namespace StringUtils;

int main() {
  try {
    {
      auto res = partition("", "");
      NPP_ASSERT(res.size() == 3);

      for (const auto& p : res) {
        NPP_ASSERT(p.empty());
      }
    }

    {
      auto res = partition("", "abc");
      NPP_ASSERT(res.size() == 3);

      for (const auto& p : res) {
        NPP_ASSERT(p.empty());
      }
    }

    {
      auto res = partition("abc", "");
      NPP_ASSERT(res.size() == 3);
      int i = 0;
      for (const auto& p : res) {
        if (i == 0) {
          NPP_ASSERT(p == "abc");
        } else {
          NPP_ASSERT(p.empty());
        }
        ++i;
      }
    }

    {
      auto res = partition("ab.c", ".");
      NPP_ASSERT(res.size() == 3);
      NPP_ASSERT(res[0] == "ab");
      NPP_ASSERT(res[1] == ".");
      NPP_ASSERT(res[2] == "c");
    }

    {
      auto res = partition("ab.c.d", ".");
      NPP_ASSERT(res.size() == 3);
      NPP_ASSERT(res[0] == "ab");
      NPP_ASSERT(res[1] == ".");
      NPP_ASSERT(res[2] == "c.d");
    }

    {
      auto res = rpartition("ab.c", ".");
      NPP_ASSERT(res.size() == 3);
      NPP_ASSERT(res[0] == "ab");
      NPP_ASSERT(res[1] == ".");
      NPP_ASSERT(res[2] == "c");
    }

    {
      auto res = rpartition("ab.c.d", ".");
      NPP_ASSERT(res.size() == 3);
      NPP_ASSERT(res[0] == "ab.c");
      NPP_ASSERT(res[1] == ".");
      NPP_ASSERT(res[2] == "d");
    }

    {
      auto res = split("", "");
      NPP_ASSERT(res.size() == 1 && res.front().empty());
    }

    {
      auto res = split("abc", "");
      NPP_ASSERT(res.size() == 1 && res.front() == "abc");
    }

    {
      auto res = split("", "abc");
      NPP_ASSERT(res.size() == 1 && res.front().empty());
    }

    {
      auto res = split("iabcxabcyabc", "abc");
      NPP_ASSERT(res.size() == 4);
      NPP_ASSERT(res[0] == "i");
      NPP_ASSERT(res[1] == "x");
      NPP_ASSERT(res[2] == "y");
      NPP_ASSERT(res[3].empty());
    }

    {
      auto res = join({"", ""}, "");
      NPP_ASSERT(res.empty());
    }

    {
      auto res = join({""}, " ");
      NPP_ASSERT(res.empty());
    }
    {
      auto res = join({"x"}, "abc");
      NPP_ASSERT(res == "x");
    }
    {
      auto res = join({"x", "y", "z"}, "");
      NPP_ASSERT(res == "xyz");
    }
    {
      auto res = join({"x", "y", "z"}, " ");
      NPP_ASSERT(res == "x y z");
    }
    {
      auto res = endsWith("", "");
      NPP_ASSERT(res);
    }
    {
      auto res = endsWith("", "a");
      NPP_ASSERT(!res);
    }
    {
      auto res = endsWith("a", "");
      NPP_ASSERT(res);
    }
    {
      auto res = endsWith("abcd", "ed");
      NPP_ASSERT(!res);
    }
    {
      auto res = endsWith("abcd", "bcd");
      NPP_ASSERT(res);
    }
    {
      auto res = endsWith("abcd", "bbbcd");
      NPP_ASSERT(!res);
    }
    {
      auto res = startsWith("abcd", "bbbcd");
      NPP_ASSERT(!res);
    }
    {
      auto res = startsWith("abcd", "abe");
      NPP_ASSERT(!res);
    }
    {
      auto res = startsWith("abcd", "ab");
      NPP_ASSERT(res);
    }
    {
      auto res = rtrim("");
      NPP_ASSERT(res.empty());
    }
    {
      auto res = rtrim("       ");
      NPP_ASSERT(res.empty());
    }
    {
      auto res = rtrim("      \n \t \r ");
      NPP_ASSERT(res.empty());
    }
    {
      auto res = rtrim("     abc \n \t \r ");
      NPP_ASSERT(res == "     abc");
    }
    {
      auto res = ltrim("     abc \n \t \r ");
      NPP_ASSERT(res == "abc \n \t \r ");
    }
    {
      auto res = trim("     abc \n \t \r ");
      NPP_ASSERT(res == "abc");
    }
  } catch (const Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return EXIT_SUCCESS;
}