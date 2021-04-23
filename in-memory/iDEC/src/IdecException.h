#ifndef _IDEC_EXCEPTION_H_
#define _IDEC_EXCEPTION_H_
#include <stdexcept>

class IdecException : public std::runtime_error {
 public:
  IdecException(const std::string &rFileName,  // filename
                unsigned int nLineNumber,      // line number
                const std::string &rMessage    // error message
                )
      : std::runtime_error(rFileName + ":" + std::to_string(nLineNumber) +
                           ": " + rMessage) {}
};

#define IDEC_REQUIRED(C)                                                   \
  do {                                                                     \
    if (!(C)) throw IdecException(__FILE__, __LINE__, #C " is required!"); \
  } while (false)

#define IDEC_REQUIRED_MSG(C, M)                                        \
  do {                                                                 \
    if (!(C))                                                          \
      throw IdecException(                                             \
          __FILE__, __LINE__,                                          \
          std::string(#C " is required! Message: ") + std::string(M)); \
  } while (false)

#endif  // _IDEC_EXCEPTION_H_