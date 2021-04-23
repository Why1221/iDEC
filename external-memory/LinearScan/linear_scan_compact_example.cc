#include "Hdf5File.h"
#include "linear_scan_compact.h"
#include <string>

void ls_demo(const std::string &h5filename, int k, int dim) {

  Hdf5File h5(h5filename, Hdf5File::Mode::ReadOnly);

  std::vector<uint64_t> train, query;

  h5.read<uint64_t>(train, "train");
  h5.read<uint64_t>(query, "test");

  const int enc_dim = dim / 64;

  const int n = train.size() / enc_dim;
  const int qn = query.size() / enc_dim;

  for (int i = 0;i < 5;++ i) {
    for (int j = 0;j < enc_dim;++ j) fprintf(stdout, "%lu ", train[i*enc_dim +j]);
    fprintf(stdout, "\n");
  }

  return;

  FILE *fp = fopen("demo-results.txt", "w");
  if (!fp) {
    perror("fopen() failed");
  }

  for (int j = 0; j < qn; ++j) {
    std::vector<std::pair<size_t, float>> ans;
    linear_scan_compact(&train[0], n, &query[j * enc_dim], enc_dim, k, ans);
    fprintf(fp, "%d", j);
    for (int i = 0; i < k; i++) {
    //   std::cout << "ret_index[" << i << "]=" << ans[i].first
    //             << " dist=" << ans[i].second << std::endl;
      fprintf(fp, " %lu %.0f", ans[i].first, ans[i].second);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

int main() {

  const std::string sample_ds = "glove-hamming-128.h5";
  const int k = 5;
  const int dim = 128;
  ls_demo(sample_ds, k, dim);

  return 0;
}