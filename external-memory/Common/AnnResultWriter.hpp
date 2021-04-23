#ifndef __ANNRESULT_HPP__
#define __ANNRESULT_HPP__

#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "Exception.h"

namespace AnnResults {

// default header & format for external memory implementations
// query id, k id (id for the k-th nearest neighbor), distance for the k-th nearest neighbor,
// ground truth distance, query time (microseconds), total IOs
const char * _DEFAULT_HEADER_E_ = "#qid,#kid,#rid,rdist,gdist,ratio,qtime(us),#io";
const char * _DEFAULT_FMT_E_ = "iiiiiffi";
const char * _DEFAULT_HEADER_I_ = "#qid,#kid,#rid,rdist,gdist,ratio,qtime(us)";
const char * _DEFAULT_FMT_I_ = "iiiiiff";
} // namespace AnnResults

class AnnResultWriter {
public:
  AnnResultWriter(const std::string &filename, bool allowOverwrite = false)
      : _fp(nullptr), _filename_cp(filename) {
#ifdef DEBUG
    fprintf(stdout, "Trying to open file %s\n", _filename_cp.c_str());
#endif
    if (!allowOverwrite && _exists()) {
      throw npp::Exception("AnnResultWriter::AnnResultWriter(): file " +
                               _filename_cp + " already exists", __FILE__, __LINE__);
    }
    if ((_fp = fopen(_filename_cp.c_str(), "w")) == nullptr)
      throw std::runtime_error(
          "AnnResultWriter::AnnResultWriter(): Failed to open file " +
          _filename_cp);
  }

  AnnResultWriter(const AnnResultWriter &) = delete;
  AnnResultWriter &operator=(const AnnResultWriter &) = delete;

  ~AnnResultWriter() {
#ifdef DEBUG
    fprintf(stdout, "Close file %s\n", _filename_cp.c_str());
#endif
    if (_fp)
      fclose(_fp);
  }

  bool writeRow(const char *fmt, ...) {
    assert(_fp != nullptr && "File MUST be opened");

    int count = 0;
    bool success = true;
    va_list args;
    va_start(args, fmt);

    while (*fmt != '\0') {
      switch (*fmt) {
      case 'i': {
        long long i = va_arg(args, long long);
        if (!_writeOne("%lld", i, (count == 0)))
          success = false;
        break;
      }
      case 'f':
      case 'd': {
        double d = va_arg(args, double);
        if (!_writeOne("%.6f", d, (count == 0)))
          success = false;
        break;
      }
      case 'c': {
        int c = va_arg(args, int);
        if (!_writeOne("%c", c, (count == 0)))
          success = false;
        break;
      }
      case 's': {
        char *s = va_arg(args, char *);
        if (!_writeOne("%s", s, (count == 0)))
          success = false;
        break;
      }
      default:
        throw std::runtime_error("Unsupported format \'" +
                                 std::to_string(*fmt) + "\'\n");
      }
      if (!success)
        break;
      ++fmt;
      ++count;
    }

    va_end(args);
    fprintf(_fp, "\n");

    return success;
  }

private:
  template <typename AnyPrintableType>
  bool _writeOne(const char *format, const AnyPrintableType &val,
                 bool isFirst = false) {

    int r1 = fprintf(_fp, "%s", (isFirst ? "" : ","));
    int r2 = fprintf(_fp, format, val);

    auto success = (r1 >= 0 && r2 >= 0);
#ifdef DEBUG
    if (!success) {
      perror("AnnResultWriter::_wrietOne() failed.\nError: ");
    }

#endif
    return success;
  }
  bool _exists() const {
    FILE *fp = fopen(_filename_cp.c_str(), "r");
    bool ex = (fp != nullptr);
    if (ex)
      fclose(fp);

    return ex;
  }
  FILE *_fp;
  std::string _filename_cp;
};

#endif
