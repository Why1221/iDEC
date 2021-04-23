#include <iostream>
#include "Exception.h"
#include "FilenameUtils.hpp"

using namespace FilenameUtils;
using namespace std;
using namespace npp;

int main() {
  try {
    {
      auto res = pathJoin("a/b", "c");
      NPP_ASSERT(res == "a/b/c");
    }
    {
      auto res = pathJoin("a/b/", "c");
      NPP_ASSERT(res == "a/b/c");
    }

    {
      auto res = getBaseName("a/b/x.txt");
      NPP_ASSERT(res == "x");
    }

    {
      auto res = getBaseName("a/b/x");
      NPP_ASSERT(res == "x");
    }

    {
      auto res = getExtension("a/b/x.txt");
      NPP_ASSERT(res == ".txt");
    }
    {
      auto res = getExtension("a/b/x");
      NPP_ASSERT(res.empty());
    }
    {
      auto res = getExtension("a/b/x.tar.gz");
      NPP_ASSERT(res == ".gz");
    }
    {
      auto res = getName("a/b/x.tar.gz");
      NPP_ASSERT(res == "x.tar.gz");
    }
    {
      auto res = getName("x.tar.gz");
      NPP_ASSERT(res == "x.tar.gz");
    }
    {
      auto res = getPath("x.tar.gz");
      NPP_ASSERT(res.empty());
    }
    {
      auto res = getPath("a/b/x.tar.gz");
      NPP_ASSERT(res == "a/b");
    }
    {
      auto res = pathSplit("a/b/x.tar.gz");
      NPP_ASSERT(res.first == "a/b" && res.second == "x.tar.gz");
    }
    {
      auto res = pathSplit("a/");
      NPP_ASSERT(res.first == "a" && res.second.empty());
    }
    {
      auto res = pathextSplit("a/");
      NPP_ASSERT(res.first == "a/" && res.second.empty());
    }
    {
      auto res = pathextSplit("a/a.tz");
      NPP_ASSERT(res.first == "a/a" && res.second == ".tz");
    }
  } catch (const Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return EXIT_SUCCESS;
}