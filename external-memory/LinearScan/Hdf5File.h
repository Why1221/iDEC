#ifndef __HDF5FILE_H_
#define __HDF5FILE_H_

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <memory>
#include <string>

using namespace HighFive;
class Hdf5File {
 public:
  enum class Mode { ReadOnly, ReadWrite };

  // throws when failed
  explicit Hdf5File(const std::string &fname, Mode mode = Mode::ReadWrite) {
    unsigned h5mode = File::ReadWrite | File::Create;

    if (mode == Mode::ReadOnly) {
      h5mode = File::ReadOnly;
    }
    _fp = std::make_unique<File>(fname, h5mode);
  }

  // get the dimensions of the dataset `dataset_name`
  // throws if failed
  std::vector<size_t> getDimensions(const std::string &dataset_name) const {
    HighFive::DataSet dataset = _fp->getDataSet(dataset_name);
    return dataset.getSpace().getDimensions();
  }

  template <typename num_t>
  void read(std::vector<std::string> data, const std::string &dataset_name) {
    auto dataset = _fp->getDataSet(dataset_name);
    dataset.read(&data);
  }

  /* Slicing reading APIs: Only support 2d data */

  // Consider the data as a matrix nXd matrix, then this API reads
  // data[0..n, ...], i.e., the first `n` rows
  template <typename num_t>
  void read(std::vector<num_t> &data, size_t n,
            const std::string &dataset_name) {
    read<num_t>(data, 0, n, 0, -1, dataset_name);
  }

  // reads data[n_start..n_end, ...], i.e., the `n_start`-th row to
  // `n_end`-th row (excluding)
  template <typename num_t>
  void read2d(std::vector<num_t> &data, size_t n_start, size_t n_end,
              const std::string &dataset_name) {
    read<num_t>(data, n_start, n_end, 0, -1, dataset_name);
  }

  // reads data[n_start..n_end, d_start..d_end]
  template <typename num_t>
  void read2d(std::vector<num_t> &data, size_t n_start, size_t n_end,
              size_t d_start, size_t d_end, const std::string &dataset_name) {
    assert(n_start < n_end && d_start < d_end);

    auto dataset = _fp->getDataSet(dataset_name);
    auto dims = dataset.getSpace().getDimensions();
    if (dims.size() != 2) {
      throw HighFive::Exception("This API only supports 2d data");
    }
    if (n_end == std::string::npos) {
      n_end = dims.front();
    } else if (n_end > dims.front()) {
      throw std::out_of_range("The dataset ONLY has " +
                              std::to_string(dims.front()) +
                              " data points. Cannot read more than that");
    }
    if (d_end == std::string::npos) {
      d_end = dims.back();
    } else if (d_end > dims.back()) {
      throw std::out_of_range(
          "The data points in this dataset have a dimension of " +
          std::to_string(dims.back()) + ", which is less than what you want");
    }
    size_t n_counts = n_end - n_start;
    size_t d_counts = d_end - d_start;
    size_t flat_dim = n_counts * d_counts;
    // allocate space
    data.resize(flat_dim);
    dataset.select({n_start, d_start}, {n_counts, d_counts}).read(&data[0]);
  }

  // reads data[n_start..n_end, ...], i.e., the `n_start`-th row to
  // `n_end`-th row (excluding)
  template <typename num_t>
  void read(std::vector<num_t> &data, size_t n_start, size_t n_end,
            const std::string &dataset_name) {
    assert(n_start < n_end);

    auto dataset = _fp->getDataSet(dataset_name);
    auto dims = dataset.getSpace().getDimensions();
    if (dims.size() != 1) {
      throw HighFive::Exception("This API only supports 1d data");
    }
    if (n_end == std::string::npos) {
      n_end = dims.front();
    } else if (n_end > dims.front()) {
      throw std::out_of_range("The dataset ONLY has " +
                              std::to_string(dims.front()) +
                              " data points. Cannot read more than that");
    }
    size_t n_counts = n_end - n_start;
    size_t flat_dim = n_counts;
    // allocate space
    data.resize(flat_dim);
    dataset.select({n_start, 0}, {n_counts, 1}).read(&data[0]);
  }

  // reads all the data (any-dimensional) in a flat vector
  // returns the dimensions of the data
  template <typename num_t>
  std::vector<size_t> read(std::vector<num_t> &data,
                           const std::string &dataset_name) {
    auto dataset = _fp->getDataSet(dataset_name);
    auto dims = dataset.getSpace().getDimensions();
    size_t flat_dim = 1ull;
    for (auto dim : dims) {
      flat_dim *= dim;
    }
    // allocate space
    data.resize(flat_dim);
    dataset.read(&data[0]);
    return dims;
  }

  void read(std::vector<std::string> &data, const std::string &dataset_name) {
    auto dataset = _fp->getDataSet(dataset_name);
    dataset.read(data);
  }

  void read(std::string &data, const std::string &dataset_name) {
    auto dataset = _fp->getDataSet(dataset_name);
    dataset.read(data);
  }

  template <typename num_t>
  void read(num_t &data, const std::string &dataset_name) {
    auto dataset = _fp->getDataSet(dataset_name);
    dataset.read(&data);
  }

  // reads 2d data as it is
  template <typename num_t>
  void read(std::vector<std::vector<num_t>> &data,
            const std::string &dataset_name) {
    auto dataset = _fp->getDataSet(dataset_name);
    auto dims = dataset.getSpace().getDimensions();
    if (dims.size() != 2) {
      // Only support reading 2d array into 2d vector,
      // as we do not have an obvious to transform
      // higher-dimensional data into 2d
      throw HighFive::Exception("Dataset dimension mismatch");
    }
    // load data directly
    dataset.read(data);
  }

  // reads data[...n, ...]
  template <typename num_t>
  void read(std::vector<std::vector<num_t>> &data, size_t n,
            const std::string &dataset_name) noexcept {
    read<num_t>(data, 0, n, 0, -1, dataset_name);
  }

  // reads data[n_start..n_end, ...]
  template <typename num_t>
  void read(std::vector<std::vector<num_t>> &data, size_t n_start, size_t n_end,
            const std::string &dataset_name) noexcept {
    read<num_t>(data, n_start, n_end, 0, -1, dataset_name);
  }

  // reads data[n_start..n_end, d_start..d_end]
  template <typename num_t>
  void read(std::vector<std::vector<num_t>> &data, size_t n_start, size_t n_end,
            size_t d_start, size_t d_end, const std::string &dataset_name) {
    assert(n_start < n_end && d_start < d_end);
    auto dataset = _fp->getDataSet(dataset_name);
    auto dims = dataset.getSpace().getDimensions();
    if (dims.size() != 2) {
      // Only support reading 2d array into 2d vector,
      // as we do not have an obvious to transform
      // higher-dimensional data into 2d
      throw HighFive::Exception("Dataset dimension mismatch");
    }
    if (n_end == -1) {
      n_end = dims.front();
    } else if (n_end > dims.front()) {
      throw std::out_of_range("The dataset ONLY has " +
                              std::to_string(dims.front()) +
                              " data points. Cannot read more than that");
    }
    if (d_end == -1) {
      d_end = dims.back();
    } else if (d_end > dims.back()) {
      throw std::out_of_range(
          "The data points in this dataset have a dimension of " +
          std::to_string(dims.back()) + ", which is less than what you want");
    }
    size_t n_counts = n_end - n_start;
    size_t d_counts = d_end - d_start;
    // load data directly
    dataset.select({n_start, d_start}, {n_counts, d_counts}).read(data);
  }

  bool exists(const std::string &dataset_name) const {
    return _fp->exist(dataset_name);
  }

  // create a dataset for data type `num_t`
  // throws if failed
  template <typename num_t>
  HighFive::DataSet createDataSet(const std::string &dataset_name,
                                  const std::vector<size_t> &dims) {
    // Only support writing the data as it is
    if (exists(dataset_name)) {
      throw HighFive::Exception("A dataset named " + dataset_name +
                                " already exists");
    }
    return _fp->createDataSet<num_t>(dataset_name, HighFive::DataSpace(dims));
  }

  // creates a dataset named as `dataset_name` and writes `data` to it
  template <typename num_t>
  void write(const std::vector<num_t> &data, const std::string &dataset_name) {
    std::vector<size_t> store_dims;
    store_dims.push_back(data.size());

    // Only support writing the data as it is
    HighFive::DataSet dataset = createDataSet<num_t>(dataset_name, store_dims);
    dataset.write(data);
  }

  // MSVC C2668: Ambiguous call to overloaded function
  // This API might not work on MSVC, details can be found at
  // https://github.com/BlueBrain/HighFive/pull/193
  template <typename num_t>
  void write(const std::vector<std::vector<num_t>> &data,
             const std::string &dataset_name) {
    std::vector<size_t> store_dims = {data.size(), data.front().size()};

    // Only support writing the data as it is
    HighFive::DataSet dataset = createDataSet<num_t>(dataset_name, store_dims);
    dataset.write(data);
  }

  template <typename num_t>
  void overWrite(const std::vector<std::vector<num_t>> &data,
                 const std::string &dataset_name) {
    // Only support writing the data as it is
    HighFive::DataSet dataset = _fp->getDataSet(dataset_name);
    dataset.write(data);
  }

  // writes data to A[n_start..n_end, d_start..d_end]
  // creates dataset if not exists, otherwise overwrites
  template <typename num_t>
  void write(const std::vector<num_t> &data, size_t n_start, size_t n_end,
             const std::string &dataset_name) {
    assert(n_start < n_end);

    auto dataset = _fp->getDataSet(dataset_name);
    size_t n_counts = n_end - n_start;
    dataset.select({n_start, 0}, {n_counts, 1}).write(data);
  }

  template <typename num_t>
  void write(const std::vector<std::vector<num_t>> &data, size_t n_start,
             size_t n_end, const std::string &dataset_name) {
    write<num_t>(data, n_start, n_end, 0, -1, dataset_name);
  }

  // writes data to A[n_start..n_end, d_start..d_end]
  // creates dataset if not exists, otherwise overwrites
  template <typename num_t>
  void write(const std::vector<std::vector<num_t>> &data, size_t n_start,
             size_t n_end, size_t d_start, size_t d_end,
             const std::string &dataset_name) {
    assert(n_start < n_end && d_start < d_end);

    auto dataset = _fp->getDataSet(dataset_name);
    auto dims = dataset.getSpace().getDimensions();
    if (n_end > dims.front()) {
      throw std::out_of_range("The dataset ONLY has " +
                              std::to_string(dims.front()) +
                              " data points. Cannot read more than that");
    }
    if (d_end == -1) {
      d_end = dims.back();
    } else if (d_end > dims.back()) {
      throw std::out_of_range(
          "The data points in this dataset have a dimension of " +
          std::to_string(dims.back()) + ", which is less than what you want");
    }
    size_t n_counts = n_end - n_start;
    size_t d_counts = d_end - d_start;
    dataset.select({n_start, d_start}, {n_counts, d_counts}).write(data);
  }

  ~Hdf5File() = default;

 private:
  std::unique_ptr<HighFive::File> _fp;
};

#endif  // __HDF5FILE_H_