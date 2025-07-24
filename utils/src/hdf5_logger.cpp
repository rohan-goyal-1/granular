#include "utils/include/hdf5_logger.h"
#include "utils/include/logger.h"

HDF5Logger::HDF5Logger (const std::string& filename)
    : file(filename, H5F_ACC_TRUNC) {}

HDF5Logger::~HDF5Logger () {
    file.close();
}

void HDF5Logger::create_group (const std::string& path) {
    file.createGroup(path);
}

template<typename T>
void HDF5Logger::write_dataset (const std::string& name, const std::vector<hsize_t>& dims, const std::vector<T>& flat_data, bool compress, int compression_level) {
    if (dims.empty())
        throw std::runtime_error("Empty dimensions not allowed.");

    hsize_t total_size = 1;
    for (auto d : dims) total_size *= d;
    if (flat_data.size() != total_size)
        throw std::runtime_error("Flat data size doesn't match dimensions.");

    H5::DataSpace dataspace(dims.size(), dims.data());

    H5::DSetCreatPropList plist;
    if (compress) {
        plist.setDeflate(compression_level);
        plist.setChunk(dims.size(), dims.data());
    }

    H5::DataSet dataset = file.createDataSet(name, get_type<T>(), dataspace, plist);
    dataset.write(flat_data.data(), get_type<T>());
}

template<typename T>
void HDF5Logger::write_attribute (const std::string& path, const std::string& attr_name, const T& value) {
    hsize_t dims = 1;
    H5::DataSpace attr_space(1, &dims);
    H5::Attribute attr;

    if (path == "/") {
        attr = file.createAttribute(attr_name, get_type<T>(), attr_space);
    }
    else if (file.nameExists(path)) {
        H5O_info_t info;
        H5Oget_info_by_name3(file.getId(), path.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT);

        if (info.type == H5O_TYPE_GROUP) {
            H5::Group group = file.openGroup(path);
            attr = group.createAttribute(attr_name, get_type<T>(), attr_space);
        }
        else if (info.type == H5O_TYPE_DATASET) {
            H5::DataSet dset = file.openDataSet(path);
            attr = dset.createAttribute(attr_name, get_type<T>(), attr_space);
        }
        else {
            throw std::runtime_error("Unsupported HDF5 object type.");
        }
    }
    else {
        throw std::runtime_error("Object path does not exist: " + path);
    }

    attr.write(get_type<T>(), &value);
}

template<>
H5::PredType HDF5Logger::get_type<int> ()    { return H5::PredType::NATIVE_INT; }
template<>
H5::PredType HDF5Logger::get_type<float> ()  { return H5::PredType::NATIVE_FLOAT; }
template<>
H5::PredType HDF5Logger::get_type<double> () { return H5::PredType::NATIVE_DOUBLE; }

// Explicit instantiation
template void HDF5Logger::write_dataset<int>(const std::string&, const std::vector<hsize_t>&, const std::vector<int>&, bool, int);
template void HDF5Logger::write_dataset<float>(const std::string&, const std::vector<hsize_t>&, const std::vector<float>&, bool, int);
template void HDF5Logger::write_dataset<double>(const std::string&, const std::vector<hsize_t>&, const std::vector<double>&, bool, int);

template void HDF5Logger::write_attribute<int>(const std::string&, const std::string&, const int&);
template void HDF5Logger::write_attribute<float>(const std::string&, const std::string&, const float&);
template void HDF5Logger::write_attribute<double>(const std::string&, const std::string&, const double&);

