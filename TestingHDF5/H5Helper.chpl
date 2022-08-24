/* IO Routines */
prototype module H5Helper {
  public use H5Wrappers;
  use FileSystem;
  // Define a type for reading h5py compatible complex numbers
  record H5PyComplexType {
    var hid : hid_t;
    proc init(std=true) {
      var Type = if std then H5T_IEEE_F64LE else H5T_NATIVE_DOUBLE;
      hid = H5Tcreate(H5T_COMPOUND, 16);
      H5Tinsert(hid, "r", 0, Type);
      H5Tinsert(hid, "i", 8, Type);
    }
    proc deinit() {
      H5Tclose(hid);
    }
  }
  // Read a real array -- these should be merged?????
  proc readRealArray(loc, arrName : string, ref Dom : domain, ref arr : [Dom] ?T) {
    var dset = openH5DataSet(loc, arrName);
    var dspace = getDataSpace(dset);
    if (dspace.rank() != Dom.rank) then halt("Rank of domain and dataspace do not agree");
    var dims = dspace.dims();
    var tup : Dom.rank*range;
    for (tup1, dim1) in zip(tup, dims) do tup1 = 0.. #dim1;
    Dom = {(...tup)};
    select T {
      when real(32) do readH5DataSet(c_ptrTo(arr), dset, H5T_NATIVE_FLOAT);
      when real(64) do readH5DataSet(c_ptrTo(arr), dset, H5T_NATIVE_DOUBLE);
      when uint(32) do readH5DataSet(c_ptrTo(arr), dset, H5T_NATIVE_UINT);
      otherwise halt("Unknown type");
    }
    return dset;
  }
  // Read a complex array
  proc readComplexArray(loc, arrName : string, ref Dom : domain, ref arr : [Dom] complex) {
    var dset = openH5DataSet(loc, arrName);
    var dspace = getDataSpace(dset);
    if (dspace.rank() != Dom.rank) then halt("Rank of domain and dataspace do not agree");
    var dims = dspace.dims();
    var tup : Dom.rank*range;
    for (tup1, dim1) in zip(tup, dims) do tup1 = 0.. #dim1;
    Dom = {(...tup)};
    var dtype = new H5PyComplexType(false);
    var status = readH5DataSet(c_ptrTo(arr), dset, dtype.hid);
    return dset;
  }
  // Read a complex array
  proc readComplexArray(dset, ref Dom : domain, ref arr : [Dom] complex) {
    var dspace = getDataSpace(dset);
    if (dspace.rank() != Dom.rank) then halt("Rank of domain and dataspace do not agree");
    var dims = dspace.dims();
    var tup : Dom.rank*range;
    for (tup1, dim1) in zip(tup, dims) do tup1 = 0.. #dim1;
    Dom = {(...tup)};
    var dtype = new H5PyComplexType(false);
    var status = readH5DataSet(c_ptrTo(arr), dset, dtype.hid);
  }
  // Save a complex array
  proc writeComplexArray(loc, arrName : string, arr : [?Dom] complex) {
    var memtype = new H5PyComplexType(false);
    var dtype = new H5PyComplexType(true);
    var dspace = createSimpleDataSpace(Dom);
    var dset = createH5DataSet(loc,arrName,dtype.hid,dspace);
    var status = writeH5DataSet(c_ptrTo(arr),dset,memtype.hid);
    return dset;
  }
  // Save a real array
  proc writeRealArray(loc, arrName : string, arr : [?Dom] real) {
    var dspace = createSimpleDataSpace(Dom);
    var dset = createH5DataSet(loc,arrName,H5T_IEEE_F64LE,dspace);
    var status = writeH5DataSet(c_ptrTo(arr),dset,H5T_NATIVE_DOUBLE);
    return dset;
  }
  proc writeScalarAttribute(loc, name : string, in value : ?T ) {
    var dspace = createScalarDataSpace();
    var attr = createH5Attribute(loc, name, TypeToH5(T,true), dspace);
    var status = writeH5Attribute(c_ptrTo(value), attr, TypeToH5(T));
  }
  proc readScalarAttribute(loc, name : string, type T) {
    var attr = openH5Attribute(loc, name);
    var value : T;
    var status = readH5Attribute(c_ptrTo(value), attr, TypeToH5(T));
    return value;
  }
  proc readArrayAttribute(loc, name, type T, param rank : int) {
    var attr = openH5Attribute(loc, name);
    var dspace = getDataSpace(attr);
    if (dspace.rank() != rank) then halt("Incorrect attribute rank");
    var dims = dspace.dims();
    var tup : rank*range;
    for param ii in 0..<rank {
      tup(ii) = 0.. #dims[ii];
    }
    const Dom = {(...tup)};
    var value : [Dom] T;
    var status = readH5Attribute(c_ptrTo(value), attr, TypeToH5(T));
    return value;
  }
  /* Open a file safely */
  proc safeOpen(fn : string, clobber : bool) {
    if (exists(fn)&&!clobber) {
      return openH5File(fn, H5F_ACC_RDWR);
    }
    return createH5File(fn, H5F_ACC_TRUNC);
  }
  /* Create empty dataset to fill */
  proc createEmptyDataset(ff : shared H5File, arrName : string, arr : [?Dom] ?T) {
    var dspace = createSimpleDataSpace(Dom);
    if (T == complex) {
      var dtype = new H5PyComplexType(true);
      return createH5DataSet(ff, arrName, dtype.hid, dspace);
    } else {
      return createH5DataSet(ff, arrName, H5T_IEEE_F64LE, dspace);
    }
  }
  /* Write a block distributed array to an HDF5 file.
     For simplicity, we do not use the HDF5 Parallel IO features, but just
     write this in series. This means we must continually open and close the
     same HDF5 file, so this routine takes in a filename as input.
     For now, the path to the dataset must either exist, or the dataset
     must be in the root path. No attempt is made to make parent directories.
  */
  proc writeBlockArray(fn : string, arrName : string, arr : [?Dom] ?T, clobber : bool = false,
                       createDataSpace : bool = true)
    where (Dom.hasSingleLocalSubdomain())
  {
    if (createDataSpace) {
      var ff = safeOpen(fn, clobber);
      createEmptyDataset(ff, arrName, arr);
      // Data set should close
    }
    // Now loop over all target locales, writing the piece we own.
    for loc in Dom.targetLocales() do on loc {
        const localfn = fn;
        var locDom = Dom.localSubdomain();
        var memspace = createSimpleDataSpace(locDom);
        var ff = openH5File(localfn, H5F_ACC_RDWR);
        var dset = openH5DataSet(ff,arrName);

        var dspace = getDataSpace(dset);
        var status = createHyperSlab(dspace, locDom);
        if (T==complex) {
          var dtype = new H5PyComplexType(false);
          status = writeH5DataSet(c_ptrTo(arr.localAccess[locDom.first]),
                                  dset, dtype.hid,
                                  memspace.id, dspace.id);
        } else {
          status = writeH5DataSet(c_ptrTo(arr.localAccess[locDom.first]),
                                  dset, H5T_NATIVE_DOUBLE,
                                  memspace.id, dspace.id);
        }
      }
    // End writeBlockArray
  }
  /* Read a block distributed array to an HDF5 file.
     Currently reads this in serially.
     TODO : Check the HDF5 dataset dimensions match block array
  */
  proc readBlockArray(fn : string, arrName : string, arr : [?Dom] ?T)
    where (Dom.hasSingleLocalSubdomain())
  {
    // Now loop over all target locales, writing the piece we own.
    for loc in Dom.targetLocales() do on loc {
        const localfn = fn;
        var locDom = Dom.localSubdomain();
        var memspace = createSimpleDataSpace(locDom);
        var ff = openH5File(localfn, H5F_ACC_RDONLY);
        var dset = openH5DataSet(ff,arrName);

        var dspace = getDataSpace(dset);
        var status = createHyperSlab(dspace, locDom);
        if (T==complex) {
          var dtype = new H5PyComplexType(false);
          status = readH5DataSet(c_ptrTo(arr.localAccess[locDom.first]),
                                 dset, dtype.hid,
                                 memspace.id, dspace.id);
        } else {
          status = readH5DataSet(c_ptrTo(arr.localAccess[locDom.first]),
                                 dset, H5T_NATIVE_DOUBLE,
                                 memspace.id, dspace.id);
        }
      }
    // End readBlockArray
  }
  /* The routines below are copied from the latest version of the HDF5.chpl module.
     Changes:
     - Support complex(64)
     - Changed name to nphdf5 etc to prevent a namespace collision.
  */
  /* Write the Block distributed Array `A` as an HDF5 dataset named `dsetName`
     in the file `filename` using parallel collective IO. This requires the
     hdf5-parallel library, which requires the MPI library.
     The file written by this function can be read in parallel with the
     function `hdf5ReadDistributedArray`. The write and read operations
     can use arrays distributed over different numbers of locales.
   */

    /* Read the HDF5 dataset named `dsetName` from the file `filename` into
       the distributed array `A`.  Each locale reads its local portion of
       the array from the file.
       This function can read the file that is generated by
       `hdf5WriteDistributedArray`.
       Currently only Block and Cyclic distributed arrays are supported.
     */
    proc nphdf5ReadDistributedArray(A, filename: string, dsetName: string) {
      // This function currently only supports Block and Cyclic distributed
      // arrays.  It is expected that this will be expanded to other
      // distributions in the future.
      //
      // BlockCyclic would work if it stored elements in its local blocks
      // sequentially instead of storing blocks side-by-side.
      // e.g. for 2x2 blocks A and B, store the elements of each block next
      // to each other in memory:
      // A11, A12, A21, A22, B11, B12, B21, B22
      // instead of:
      // A11, A12, B11, B12
      // A21, A22, B21, B22
      use BlockDist, CyclicDist;
      use HDF5;
      proc isBlock(D: Block) param return true;
      proc isBlock(D) param return false;
      proc isCyclic(D: Cyclic) param return true;
      proc isCyclic(D) param return false;
      if !(isBlock(A.dom.dist) || isCyclic(A.dom.dist)) then
        compilerError("hdf5ReadDistributedArray currently only supports block or cyclic distributed arrays");
      coforall loc in A.targetLocales() do on loc {
        // make sure the file name is local
        var filenameCopy = filename;
        var file_id = C_HDF5.H5Fopen(filenameCopy.c_str(),
                                     C_HDF5.H5F_ACC_RDONLY,
                                     C_HDF5.H5P_DEFAULT);
        var dataset = C_HDF5.H5Dopen(file_id, dsetName.c_str(),
                                     C_HDF5.H5P_DEFAULT);
        var dataspace = C_HDF5.H5Dget_space(dataset);
        var dsetRank: c_int;
        C_HDF5.H5LTget_dataset_ndims(file_id, dsetName.c_str(), dsetRank);
        var dims: [0.. #dsetRank] C_HDF5.hsize_t;
        {
          /*
          // This chained module call causes a compiler error.  As a
          // workaround, 'use' the module in a new scope instead.
          C_HDF5.HDF5_WAR.H5LTget_dataset_info_WAR(file_id, dsetName.c_str(),
                                                   c_ptrTo(dims), nil, nil);
          */
          use C_HDF5.HDF5_WAR;
          H5LTget_dataset_info_WAR(file_id, dsetName.c_str(),
                                   c_ptrTo(dims), nil, nil);
        }
        const wholeLow = A.domain.whole.low;
        for dom in A.localSubdomains() {
          const domIndex = getTupleIndices(dom.rank);
          // The dataset is 0-based, so unTranslate each block
          const dsetBlock = dom.chpl__unTranslate(wholeLow);
          // Arrays to represent locations with the file
          var dsetOffsetArr, dsetCountArr,
              dsetStrideArr: [domIndex] C_HDF5.hsize_t;
          // Arrays to represent locations in the distributed array
          var memOffsetArr, memCountArr,
              memStrideArr: [domIndex] C_HDF5.hsize_t;
          for i in domIndex {
            dsetOffsetArr[i] = dsetBlock.dim(i).low: C_HDF5.hsize_t;
            dsetCountArr[i]  = dsetBlock.dim(i).size: C_HDF5.hsize_t;
            dsetStrideArr[i] = dsetBlock.dim(i).stride: C_HDF5.hsize_t;
            // We'll write to a slice of the array so that offset is always 0
            memOffsetArr[i] = 0: C_HDF5.hsize_t;
            memStrideArr[i] = 1: C_HDF5.hsize_t;
            memCountArr[i] = dom.dim(i).size: C_HDF5.hsize_t;
          }
          // create the hyperslab into the dataset on disk
          C_HDF5.H5Sselect_hyperslab(dataspace, C_HDF5.H5S_SELECT_SET,
                                     c_ptrTo(dsetOffsetArr),
                                     c_ptrTo(dsetStrideArr),
                                     c_ptrTo(dsetCountArr), nil);
          // create the hyperslab into the array to read the dataset into
          var memspace = C_HDF5.H5Screate_simple(dom.rank,
                                                 c_ptrTo(memCountArr), nil);
          C_HDF5.H5Sselect_hyperslab(memspace, C_HDF5.H5S_SELECT_SET,
                                     c_ptrTo(memOffsetArr),
                                     c_ptrTo(memStrideArr),
                                     c_ptrTo(memCountArr), nil);
          ref AA = A[dom];
          const h5ComplexType = new H5PyComplexType(false);
          const hdf5Type = if (A.eltType==complex) then h5ComplexType.hid else getHDF5Type(A.eltType);
          C_HDF5.H5Dread(dataset, hdf5Type, memspace, dataspace,
                         C_HDF5.H5P_DEFAULT, c_ptrTo(AA));
          C_HDF5.H5Sclose(memspace);
        }
        // close dataspace, dataset
        C_HDF5.H5Dclose(dataset);
        C_HDF5.H5Sclose(dataspace);
        C_HDF5.H5Fclose(file_id);
      }
    }
    private inline proc getTupleIndices(param n:int) {
      var x : n*int;
      return x.indices;
    }
}
