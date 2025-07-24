#ifndef HDF5_LOGGER_H
#define HDF5_LOGGER_H

#include <string>
#include <vector>
#include <stdexcept>
#include <H5Cpp.h>

#include "include/system.h"
#include "utils/include/logger.h"

#include <H5Cpp.h>

class HDF5Logger {
public:
    explicit HDF5Logger(const std::string& filename);
    ~HDF5Logger(void);

    void create_group(const std::string& path);

    template<typename T>
    void write_dataset(const std::string& name, const std::vector<hsize_t>& dims, const std::vector<T>& flat_data, bool compress = false, int compression_level = 4);

    template<typename T>
    void write_attribute(const std::string& path, const std::string& attr_name, const T& value);

private:
    H5::H5File file;

    template<typename T>
    H5::PredType get_type(void);
};

#endif // HDF5_LOGGER_H
