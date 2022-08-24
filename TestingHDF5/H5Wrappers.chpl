/* These are just convenience wrappers around
   basic HDF5 routines.
*/
module H5Wrappers {
  public use HDF5;
  public use HDF5.C_HDF5;
  // Basic file operations
  class H5File {
    var id : hid_t;
    proc deinit() {
      var status = H5Fclose(id);
    }
  }
  proc createH5File(fn: string,
                    mode=H5F_ACC_EXCL,
                    create=H5P_DEFAULT,
                    access=H5P_DEFAULT) {
    var ret = new shared H5File();
    ret.id = H5Fcreate(fn.c_str(), mode, create, access);
    return ret;
  }
  proc openH5File(fn: string,
                  mode=H5F_ACC_RDONLY,
                  access=H5P_DEFAULT) {
    var ret = new shared H5File();
    ret.id = H5Fopen(fn.c_str(), mode, access);
    return ret;
  }
  // Basic group operations
  class H5Group {
    var id : hid_t;
    proc deinit() {
      var status = H5Gclose(id);
    }
  }
  proc createH5Group(loc,
                     name: string,
                     link=H5P_DEFAULT,
                     create=H5P_DEFAULT,
                     access=H5P_DEFAULT) {
    var ret = new shared H5Group();
    ret.id = H5Gcreate2(loc.id,name.c_str(), link, create, access);
    return ret;
  }
  proc openH5Group(loc,
                   name: string,
                   access=H5P_DEFAULT) {
    var ret = new shared H5Group();
    ret.id = H5Gopen2(loc.id, name.c_str(), access);
    return ret;
  }
  // Dataspaces
  class H5DataSpace {
    var id : hid_t;
    proc deinit() {
      var status = H5Sclose(id);
    }
    proc rank() {
      return H5Sget_simple_extent_ndims(id):int;
    }
    proc dims() {
      const D = {0.. #rank()};
      var tmp : [D] hsize_t;
      var ret : [D] int;
      var ndims = H5Sget_simple_extent_dims(id,c_ptrTo(tmp[0]),nil);
      [ii in D] ret[ii] = tmp[ii]:int;
      return ret;
    }
  }
  proc createSimpleDataSpace(d : domain) where isRectangularDom(d) {
    var ret = new shared H5DataSpace();
    const rank = d.rank : c_int;
    var dims : [0.. #rank] hsize_t;
    const shape = d.shape;
    for (dim1, shape1) in zip(dims, shape) do dim1 = shape1:hsize_t;
    ret.id = H5Screate_simple(rank, c_ptrTo(dims[0]), nil);
    return ret;
  }
  proc createScalarDataSpace() {
    var ret = new shared H5DataSpace();
    ret.id = H5Screate(H5S_SCALAR);
    return ret;
  }
  // Basic dataset operations
  class H5DataSet {
    var id : hid_t;
    proc deinit() {
      var status = H5Dclose(id);
    }
  }
  proc createH5DataSet(loc,
                       name: string,
                       dtype: hid_t,
                       dspace : borrowed H5DataSpace,
                       link=H5P_DEFAULT,
                       create=H5P_DEFAULT,
                       access=H5P_DEFAULT) {
    var ret = new shared H5DataSet();
    ret.id = H5Dcreate2(loc.id,name.c_str(),
                        dtype, dspace.id,
                        link, create, access);
    return ret;
  }
  proc openH5DataSet(loc,
                     name: string,
                     access=H5P_DEFAULT) {
    var ret = new shared H5DataSet();
    ret.id = H5Dopen2(loc.id, name.c_str(), access);
    return ret;
  }
  proc readH5DataSet(cptr,
                     dset : borrowed H5DataSet,
                     memtype: hid_t,
                     memspace = H5S_ALL,
                     filespace = H5S_ALL,
                     xfer = H5P_DEFAULT) {
    return H5Dread(dset.id, memtype, memspace, filespace, xfer, cptr);
  }
  proc writeH5DataSet(cptr,
                      dset : borrowed H5DataSet,
                      memtype: hid_t,
                      memspace = H5S_ALL,
                      filespace = H5S_ALL,
                      xfer = H5P_DEFAULT) {
    return H5Dwrite(dset.id, memtype, memspace, filespace, xfer, cptr);
  }
  proc getDataSpace(dset : borrowed H5DataSet) {
    var ret = new shared H5DataSpace();
    ret.id = H5Dget_space(dset.id);
    return ret;
  }
  // Basic attribute operations
  class H5Attribute {
    var id : hid_t;
    proc deinit() {
      var status = H5Aclose(id);
    }
  }
  proc createH5Attribute(loc,
                         name: string,
                         dtype: hid_t,
                         dspace : borrowed H5DataSpace,
                         create=H5P_DEFAULT,
                         access=H5P_DEFAULT) {
    var ret = new shared H5Attribute();
    ret.id = H5Acreate2(loc.id,name.c_str(),
                        dtype, dspace.id,
                        create, access);
    return ret;
  }
  proc openH5Attribute(loc,
                       name: string,
                       access=H5P_DEFAULT) {
    var ret = new shared H5Attribute();
    ret.id = H5Aopen(loc.id, name.c_str(), access);
    return ret;
  }
  proc readH5Attribute(cptr,
                       attr : borrowed H5Attribute,
                       memtype: hid_t) {
    return H5Aread(attr.id, memtype, cptr);
  }
  proc writeH5Attribute(cptr,
                        attr : borrowed H5Attribute,
                        memtype: hid_t) {
    return H5Awrite(attr.id, memtype, cptr);
  }
  proc getDataSpace(attr : borrowed H5Attribute) {
    var ret = new shared H5DataSpace();
    ret.id = H5Aget_space(attr.id);
    return ret;
  }
  // Define a hyperslab
  // This function sets the block size to 1, use
  // the overload below to access a different block.
  //
  // Uses domains to define sizes etc.
  //
  // Note that dataspaces start at 0
  proc createHyperSlab(dspace: borrowed H5DataSpace,
                       d : domain,
                       op = H5S_SELECT_SET)
    where isRectangularDom(d) {
    const rank = d.rank : c_int;
    var starts, counts, strides : [0..<rank] hsize_t;
    // Work around for first having different behaviour
    // for rank-1 and rank-n arrays.
    proc tuplify(t) {
      if !isTuple(t) then return (t,); else return t;
    }
    const first = tuplify(d.first);
    const stride = tuplify(d.stride);
    for ii in 0..<rank {
      starts[ii] = first(ii):hsize_t;
      counts[ii] = d.shape(ii):hsize_t;
      strides[ii] = stride(ii):hsize_t;
    }
    return H5Sselect_hyperslab(dspace.id,op,
                               c_ptrTo(starts[0]),
                               c_ptrTo(strides[0]),
                               c_ptrTo(counts[0]),
                               nil);
  }
  proc createHyperSlab(dspace: borrowed H5DataSpace,
                       d : domain,
                       blocks : [?Dom] hsize_t,
                       op = H5S_SELECT_SET)
    where isRectangularDom(d) && (Dom.size == d.rank)
  {
    const rank = d.rank : c_int;
    var starts, counts, strides : [0..<rank] hsize_t;
    for ii in 0..<rank {
      starts[ii] = d.first(ii):hsize_t;
      counts[ii] = d.shape(ii):hsize_t;
      strides[ii] = d.stride(ii):hsize_t;
    }
    return H5Sselect_hyperslab(dspace.id,op,
                               c_ptrTo(starts[0]),
                               c_ptrTo(strides[0]),
                               c_ptrTo(counts[0]),
                               c_ptrTo(blocks[Dom.first]));
  }
  // Basic Plist operations
  class H5Plist {
    var id : hid_t;
    proc deinit() {
      var status = H5Pclose(id);
    }
  }
  proc createH5Plist(cls : hid_t) {
    var ret = new shared H5Plist();
    ret.id = H5Pcreate(cls);
    return ret;
  }
  // A convenience complex type
  // Treat a complex number as 2-elt array
  record H5ComplexType {
    var hid : hid_t;
    proc init(std=true) {
      var dim : hsize_t = 2;
      if std {
        hid = H5Tarray_create2(H5T_IEEE_F64LE, 1, c_ptrTo(dim));
      } else {
        hid = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, c_ptrTo(dim));
      }
    }
    proc deinit() {
      H5Tclose(hid);
    }
  }
  // Type conversion routines
  proc TypeToH5(type t, std = true) {
    if std {
      select t {
          when int(32) do return H5T_STD_I32LE;
          when int(64) do return H5T_STD_I64LE;
          when uint(32) do return H5T_STD_U32LE;
          when uint(64) do return H5T_STD_U64LE;
          when real(32) do return H5T_IEEE_F32LE;
          when real(64) do return H5T_IEEE_F64LE;
          otherwise halt("Unknown type");
        }
    } else {
      select t {
          when int(32) do return H5T_NATIVE_INT32;
          when int(64) do return H5T_NATIVE_INT64;
          when uint(32) do return H5T_NATIVE_UINT32;
          when uint(64) do return H5T_NATIVE_UINT64;
          when real(32) do return H5T_NATIVE_FLOAT;
          when real(64) do return H5T_NATIVE_DOUBLE;
          otherwise halt("Unknown type");
        }
    }
  }
}
