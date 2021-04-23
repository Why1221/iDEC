#ifndef _FILENAME_UTILS_HPP_
#define _FILENAME_UTILS_HPP_
#include "StringUtils.hpp"

// Names of APIs copied from
// https://commons.apache.org/proper/commons-io/javadocs/api-1.4/org/apache/commons/io/FilenameUtils.html
namespace FilenameUtils {
using namespace StringUtils;
using String = std::string;

// Joins a filename to a base path using normal command line style rules.
inline std::string pathJoin(std::string basePath,
                            std::string fullFilenameToAdd);
// Gets the base name, minus the full path and extension, from a full filename.
inline std::string getBaseName(std::string filename);
// Gets the extension of a filename.
inline std::string getExtension(std::string filename);
// Gets the name minus the path from a full filename.
inline std::string getName(std::string filename);
//  Gets the path from a full filename, which excludes the prefix.
inline String getPath(String filename);
//
inline std::pair<std::string /* path */, std::string /* name */> pathSplit(
    std::string filename);
//
inline std::pair<std::string /* root */, std::string /* ext */> pathextSplit(
    std::string filename);
}  // namespace FilenameUtils

namespace FilenameUtils {
std::string pathJoin(std::string basePath, std::string fullFilenameToAdd) {
  if (basePath.back() == '/') {
    return basePath + fullFilenameToAdd;
  }
  return basePath + "/" + fullFilenameToAdd;
}
std::string getBaseName(std::string filename) {
  auto name = getName(filename);
  auto res = rpartition(name, ".");
  return res.front();
}
std::string getExtension(std::string filename) {
  auto res = rpartition(filename, ".");
  return res[1] + res.back();
}
std::string getName(std::string filename) {
  auto res = rpartition(filename, "/");
  if (res[1] != "/") return filename;
  return res.back();
}
String getPath(String filename) {
  auto res = rpartition(filename, "/");
  if (res[1] != "/") return "";
  return res.front();
}
std::pair<std::string /* path */, std::string /* name */> pathSplit(
    std::string filename) {
  auto res = rpartition(filename, "/");
  if (res[1] != "/") {
    return {"", filename};
  }
  return {res.front(), res.back()};
}
std::pair<std::string /* root */, std::string /* ext */> pathextSplit(
    std::string filename) {
  auto res = rpartition(filename, ".");
  return {res.front(), res[1] + res.back()};
}

}  // namespace FilenameUtils

#endif // _FILENAME_UTILS_HPP_