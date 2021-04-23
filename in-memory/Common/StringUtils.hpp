#ifndef _STRING_UTILS_HPP_
#define _STRING_UTILS_HPP_
#include <locale>
#include <string>
#include <unordered_set>
#include <vector>

namespace StringUtils {
using String = std::string;
// mimick Python partition for str
inline std::vector<std::string> partition(const std::string& str,
                                          const std::string& delimiter);
// mimick Python rpartition for str
inline std::vector<std::string> rpartition(const std::string& str,
                                           const std::string& delimiter);
// mimick Python split for str
std::vector<String> split(const String& input, const String& delimiter);

// mimick Python join for str
inline String join(const std::vector<String>& container,
                   const String& delimiter);
template <typename Iterator>
inline String join(Iterator first, Iterator end, const String& delimiter);

// mimick Python endswith/startswith for str
inline bool endsWith(const String& s1, const String& s2);
inline bool startsWith(const String& s1, const String& s2);

// mimick Python trim

inline String ltrim(String s, const std::unordered_set<char>& toTrim = {});

inline String rtrim(String s, const std::unordered_set<char>& toTrim = {});

inline String trim(String s, const std::unordered_set<char>& toTrim = {});

inline String toLower(String s);
inline String toUpper(String s);
}  // end namespace StringUtils

namespace StringUtils {
std::vector<std::string> partition(const std::string& str,
                                   const std::string& delimiter) {
  std::vector<std::string> pat(3);
  if (delimiter.empty()) {
    pat[0] = str;
    return pat;
  }
  auto pos = str.find(delimiter);
  pat[0] = str.substr(0, pos);
  if (pos != std::string::npos) {
    pat[1] = delimiter;
    pat[2] = str.substr(pos + delimiter.size());
  }
  return pat;
}

std::vector<std::string> rpartition(const std::string& str,
                                    const std::string& delimiter) {
  std::vector<std::string> pat(3);
  auto pos = str.rfind(delimiter);
  pat[0] = str.substr(0, pos);
  if (pos != std::string::npos) {
    pat[1] = delimiter;
    pat[2] = str.substr(pos + delimiter.size());
  }
  return pat;
}

static void _split(const String& input, const String& delimiter,
                   std::vector<String>& out) {
  auto pos = input.find(delimiter);
  if (pos != std::string::npos) {
    out.push_back(input.substr(0, pos));
    _split(input.substr(pos + delimiter.size()), delimiter, out);
  } else {
    out.push_back(input);
  }
}

std::vector<String> split(const String& input, const String& delimiter) {
  if (delimiter.empty() || input.empty()) return {input};

  std::vector<String> out;
  _split(input, delimiter, out);
  return out;
}

template <typename Iterator>
String join(Iterator first, Iterator last, const String& delimiter) {
  size_t sz = 0;
  size_t cnt = 0;
  for (auto it = first; it != last; ++it) {
    sz += it->size();
    ++cnt;
  }

  if (cnt == 0u) {
    return "";
  } else if (cnt == 1u) {
    return *first;
  } else {
    sz += (cnt - 1) * delimiter.size();
    String out(sz, '\0');
    size_t gi = 0, i = 0;
    auto it = first;
    for (i = 0; i < it->size(); ++i, ++gi) {
      out[gi] = it->at(i);
    }
    ++it;
    for (; it != last; ++it) {
      for (i = 0; i < delimiter.size(); ++i, ++gi) {
        out[gi] = delimiter.at(i);
      }
      for (i = 0; i < it->size(); ++i, ++gi) {
        out[gi] = it->at(i);
      }
    }

    return out;
  }
}

String join(const std::vector<String>& container, const String& delimiter) {
  return join(container.cbegin(), container.cend(), delimiter);
}

bool endsWith(const String& s1, const String& s2) {
  return (s1.size() >= s2.size() &&
          s1.compare(s1.size() - s2.size(), String::npos, s2) == 0);
}

bool startsWith(const String& s1, const String& s2) {
  return (s1.size() >= s2.size() && s1.compare(0, s2.size(), s2) == 0);
}

static std::unordered_set<char>& getDefaultCTrim() {
  static std::unordered_set<char> DefaultCTrim({' ', '\n', '\r', '\t'});
  return DefaultCTrim;
}


String ltrim(String s, const std::unordered_set<char>& toTrim) {
  if (s.empty()) return "";
  if (toTrim.empty()) {
    return ltrim(s, getDefaultCTrim());
  }

  size_t first = 0u;
  while (first < s.size() && toTrim.count(s.at(first)) > 0) ++first;
  return s.substr(first);
}


String rtrim(String s, const std::unordered_set<char>& toTrim) {
  if (s.empty()) return "";
  if (toTrim.empty()) {
    return rtrim(s, getDefaultCTrim());
  }

  size_t last = s.size() - 1u;
  while (last != std::string::npos && toTrim.count(s.at(last)) > 0) --last;
  if (last == std::string::npos) return "";
  return s.substr(0, last + 1);
}


String trim(String s, const std::unordered_set<char>& toTrim) {
  auto out = ltrim(s, toTrim);
  return rtrim(out, toTrim);
}

String toLower(String s) {
  std::locale loc;
  for (auto& ch : s) {
    ch = std::tolower(ch, loc);
  }
  return s;
}

String toUpper(String s) {
  std::locale loc;
  for (auto& ch : s) {
    ch = std::toupper(ch, loc);
  }
  return s;
}

}  // end namespace StringUtils

#endif // _STRING_UTILS_HPP_