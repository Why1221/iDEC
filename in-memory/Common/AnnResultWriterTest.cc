#define DEBUG

#include <iostream>
#include <limits>
#include "AnnResultWriter.hpp"


const char* default_header = "#qid,#kid,#rid,rdist,gdist";
const char* default_format = "iiiff";

int main() {
  {
    AnnResultWriter writer("ann-result-writer-demo.txt", true);
    bool success = true;

    success = success && writer.writeRow("s", default_header);
    success = success && writer.writeRow(default_format, 1, 1, 100, 12.0, 10.0);
    success = success && writer.writeRow(default_format, 1, 2, 111, 21.0, 13.0);
    success = success && writer.writeRow(default_format, 2, 1, 201, 15.0, 11.0);
    success = success && writer.writeRow(default_format, 2, 2, 303, 19.0, 12.0);
    success = success && writer.writeRow(default_format, 3, 1, 404, 22.0, 17.0);
    success = success && writer.writeRow(default_format, 3, 2, 503, 30.0, 23.0);

    success = success && writer.writeRow(
                             default_format, 3, 2,
                             (long long)std::numeric_limits<int>::max() + 100ll,
                             30.0, 23.0);

    success = success &&
              writer.writeRow(default_format, 3, 2,
                              std::numeric_limits<unsigned>::max(), 30.0, 23.0);
    assert(success);
  }
  {
    try {
      AnnResultWriter writer("ann-result-writer-demo.txt");
    } catch (const npp::Exception& e) {
      std::cout << "You should see this output: " << e << std::endl;
    }
  }
  return 0;
}