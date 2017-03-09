import os
import subprocess
import glob
import tempfile

import numpy as np
import gdal
import osr


def find_gdal_exe(gdalcmd):
    try:
        if not gdalcmd.endswith('.exe'):
            gdalcmd += '.exe'
        pattern = os.path.join('C:\\', 'OSGeo4W*', 'bin', gdalcmd)
        cmdpath = glob.glob(pattern)[0]
    except IndexError:
        cmdpath = gdalcmd
    return cmdpath


gdal_translate_exe = find_gdal_exe('gdal_translate')
gdal_rasterize_exe = find_gdal_exe('gdal_rasterize')
gdalwarp_exe = find_gdal_exe('gdalwarp')
gdalbuildvrt_exe = find_gdal_exe('gdalbuildvrt')


dtype_np_to_gdal = {
        'i2': gdal.GDT_Int16,
        'i4': gdal.GDT_Int32,
        'u1': gdal.GDT_Byte,
        'u2': gdal.GDT_Int16,
        'f4': gdal.GDT_Float32,
        'f8': gdal.GDT_Float64}

default_nodata = {
        'f4': 9.969209968386869e+36,
        'f8': 9.969209968386869e+36,
        'i1': -127,
        'i2': -32767,
        'i4': -2147483647,
        'i8': -9223372036854775806L,
        'u1': 255,
        'u2': 65535,
        'u4': 4294967295L,
        'u8': 18446744073709551614L}

# Pixel data types
# http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4
gdt_dtype = ['', 'u1', 'u2', 'i2', 'u4', 'i4', 'f4', 'f8']


class GdalErrorHandler():
    def __init__(self):
        self.err_level = gdal.CE_None
        self.err_no = 0
        self.err_msg = ''

    def handler(self, err_level, err_no, err_msg):
        self.err_level = err_level
        self.err_no = err_no
        self.err_msg = err_msg


class gdal_handle_errors:
    def __init__(self):
        pass

    def __enter__(self):
        err = GdalErrorHandler()
        handler = err.handler
        gdal.PushErrorHandler(handler)
        gdal.UseExceptions()

    def __exit__(self, type, value, traceback):
        gdal.PopErrorHandler()


def get_default_nodata(dtype):
    dt = np.dtype(dtype)
    key = dt.str[1:]
    return np.array(default_nodata[key], dt)[()]


def get_gdal_dtype(dtype):
    dt = np.dtype(dtype)
    key = dt.str[1:]
    return dtype_np_to_gdal[key]


def get_file_nodata(intif):
    """Get nodata from first band in file"""
    with gdal_open(intif) as ds:
        b = ds.GetRasterBand(1)
        src_nodata = b.GetNoDataValue()
        if src_nodata is None:
            key = gdt_dtype[b.DataType]
            src_nodata = default_nodata[key]
    return src_nodata


def _get_startupinfo():
    """startupinfo to suppress external command windows"""
    startupinfo = None
    if os.name == 'nt':
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    return startupinfo


def array_to_gtiff(arr, outfile, projection, geotransform, banddim=0,
        tgt_nodata=None, create_options=['COMPRESS=LZW', 'BIGTIFF=IF_SAFER'], dtype=None):
    """
    Save a numpy array to GTiff

    Parameters
    ----------
    arr : ndarray
        a (n, m, b) matrix, where b is the band number,
        n and m is the row and collum size of the matrices.
    outfile : str
        filename
    projection : str
        image projection
    geotransform : gdal.GetGeoTransform()
        image geotransform
    banddim :
        swapping of axis if (n, m, b) is not in the correct order
    tgt_nodata : float, int
        target nodata value
    create_options :
        passed to gdal.Create
    dtype : numpy.dtype or parseable
        enforce this dtype
        default: array dtype
    """
    err = GdalErrorHandler()
    handler = err.handler
    gdal.PushErrorHandler(handler)
    gdal.UseExceptions()
    if dtype is None:
        dtype = arr.dtype
    gdal_dtype = get_gdal_dtype(dtype)

    if tgt_nodata is None:
        tgt_nodata = get_default_nodata(dtype)
    # GDAL for some reason fails with other types
    # TODO: Find and resolve this issue
    tgt_nodata = int(tgt_nodata)

    arr = np.ma.masked_invalid(arr).filled(tgt_nodata)

    # get array into right format
    if np.ndim(arr) == 3:
        if banddim == 1:
            arr = np.swapaxes(arr, 0, 1)
            arr = np.swapaxes(arr, 1, 2)
        if banddim == 2:
            arr = np.swapaxes(arr, 1, 2)
            arr = np.swapaxes(arr, 0, 1)
    elif np.ndim(arr) == 2:
        arr = arr[np.newaxis, :, :]
        banddim = 0
    else:
        raise NotImplementedError("Need at least 2D data.")
    nbands, nx, ny = arr.shape

    if outfile == 'MEM':
        drv = gdal.GetDriverByName('MEM')
    else:
        drv = gdal.GetDriverByName('GTiff')

    out_tif = drv.Create(outfile, ny, nx, nbands, gdal_dtype, create_options)
    if out_tif is None:
        raise IOError('Unable to create new dataset in {}.'.format(outfile))

    try:
        out_tif.SetGeoTransform(geotransform)
        out_tif.SetProjection(projection)
        for b in range(nbands):
            out_tif.GetRasterBand(b+1).WriteArray(arr[b, :, :])
        out_tif.GetRasterBand(1).SetNoDataValue(tgt_nodata)
    finally:
        if outfile != 'MEM':
            out_tif = None
        gdal.PopErrorHandler()

    return out_tif


def extent_qgis_to_gdal(extent):
    """Convert extent string or list from qgis to gdal order

    Parameters
    ----------
    extent : comma-separated str or list
        xmin, xmax, ymin, ymax

    Returns
    -------
    list of str
        xmin, ymax, xmax, ymin
    """
    if ',' in extent:
        extent = extent.split(',')
    xmin, xmax, ymin, ymax = extent
    return list(map(str, [xmin, ymax, xmax, ymin]))


def extent_qgis_to_gdal_warp(extent):
    """Convert extent string or list from qgis to gdal order

    Parameters
    ----------
    extent : comma-separated str or list
        xmin, xmax, ymin, ymax

    Returns
    -------
    list of str
        xmin, ymax, xmax, ymin
    """
    if ',' in extent:
        extent = extent.split(',')
    xmin, xmax, ymin, ymax = extent
    return list(map(str, [xmin, ymin, xmax, ymax]))


def extent_from_dataset(ds):
    gt = ds.GetGeoTransform()
    xmin = gt[0]
    ymax = gt[3]
    xmax = xmin + gt[1]*ds.RasterXSize
    ymin = ymax + gt[5]*ds.RasterYSize
    return dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

def get_extent(gdalstr):
    with gdal_open(gdalstr) as ds:
        return extent_from_dataset(ds)

def cmd_gdalwarp_cutline(intif, inshp, outtif):
    dstnodata = get_file_nodata(intif)
    cmd = [gdalwarp_exe]
    cmd += ['-cutline', inshp]
    cmd += ['-crop_to_cutline']
    cmd += ['-dstnodata', str(dstnodata)]
    cmd += [intif, outtif]
    return cmd

def cmd_gdalwarp_reproject(infile, outfile, t_srs=None, r='near',
        ot=None, of=None, extent=None,
        cmd_extra=[]):
    dstnodata = get_file_nodata(infile)
    if outfile.endswith('.vrt'):
        of = 'VRT'
    cmd = [gdalwarp_exe]
    if t_srs is not None:
        cmd += ['-t_srs', t_srs]
    cmd += ['-r', r]
    if ot is not None:
        cmd += ['-ot', ot]
    if of is not None:
        cmd += ['-of', of]
    if extent:
        cmd += ['-te'] + extent_qgis_to_gdal_warp(extent)
    cmd += ['-srcnodata', str(dstnodata)]
    cmd += ['-dstnodata', str(dstnodata)]
    cmd += ['-overwrite']
    cmd += cmd_extra
    cmd += [infile, outfile]
    return cmd

def cmd_gdalwarp_reproject_cutline(intif, inshp, outtif, t_srs, r='near', cmd_extra=[]):
    dstnodata = get_file_nodata(intif)
    cmd = [gdalwarp_exe]
    cmd += ['-t_srs', t_srs]
    cmd += ['-r', r]
    cmd += ['-cutline', inshp]
    cmd += ['-crop_to_cutline']
    cmd += ['-dstnodata', str(dstnodata)]
    cmd += ['-overwrite']
    cmd += cmd_extra
    cmd += [intif, outtif]
    return cmd

def cmd_gdalbuildvrt(outfile, infiles=[], input_file_list=[],
        resolution='average', separate=False, proj_difference=False):
    cmd = [gdalbuildvrt_exe]
    cmd += ['-q']
    if resolution != 'average':
        cmd += ['-resolution', resolution]
    if separate:
        cmd += ['-separate']
    if proj_difference:
        cmd += ['-allow_projection_difference']
    cmd += []
    cmd += [outfile]
    if infiles:
        cmd += infiles
    elif input_file_list:
        cmd += ['-input_file_list', input_file_list]
    return cmd

def cmd_gdal_translate(infile, outfile,
        of='GTiff', ot='',
        extent=[], outsize=None,
        co_dict={}, extra=[]):
    """Generate gdal_translate command

    Parameters
    ----------
    infile, outfile : str
        paths to in/out files
    of : str
        output format
    ot : str
        output data type
    extent : list or str
        extent to subset
        xmin, xmax, ymin, ymax
    outsize : list of len 2
        output size
    co_dict : dict
        creation options
        e.g. {'COMPRESS': 'JPEG', 'JPEG_QUALITY': 75}
    extra : list
        extra commands
    """
    cmd = [gdal_translate_exe]
    if of.lower() != 'gtiff':
        cmd += ['-of', of]
    if ot:
        cmd += ['-ot', ot]
    if extent:
        cmd += ['-projwin'] + extent_qgis_to_gdal(extent)
    outsize = _format_outsize_arg(outsize)
    if outsize:
        cmd += ['-outsize'] + [str(s) for s in outsize]
    if co_dict:
        for k, v in co_dict.items():
            cmd += ['-co', '{}={}'.format(k, v)]
    cmd += extra
    cmd += [infile, outfile]
    return cmd


def _format_outsize_arg(outsize):
    """Format outsize argument for GDAL translate"""
    wrong_format_err = ValueError('`outsize` parameter format not understood: {}'.format(outsize))
    if not outsize:
        return None
    elif outsize in [['100%', 0], ['100%', '100%'], '100%', 100]:
        return None
    elif isinstance(outsize, list):
        if len(outsize) == 2:
            return outsize
        else:
            raise wrong_format_err
    elif isinstance(outsize, int):
        return ['{}%'.format(outsize)] * 2
    else:
        return [outsize] * 2


def run_cmd(cmd, outfile):
    subprocess.check_output(cmd,
            stdin=subprocess.PIPE, stderr=subprocess.PIPE,  # to avoid error in pythonw
            shell=True,
            startupinfo=_get_startupinfo())
    if not os.path.isfile(outfile):
        cmdstr = subprocess.list2cmdline(cmd)
        raise RuntimeError('GDAL command failed. No output created with '
                'cmd \'{}\'.'.format(cmdstr))


def make_gdalstr(fname, group=None, varn=None):
    ext = os.path.splitext(fname)[1].lower()
    if ext == '.nc' and varn is not None:
        gdalstr = 'NETCDF:{}:{}'.format(fname, varn)
    elif ext == '.hdf' and varn is not None and group is not None:
        gdalstr = 'HDF4_EOS:EOS_GRID:"{fname}":{group}:{varn}'.format(
                fname=fname, group=group, varn=varn)
    elif ext == '.h5' and varn is not None and group is not None:
        gdalstr = 'HDF5:"{fname}"//{group}/{varn}'.format(
                fname=fname, group=group, varn=varn)
    else:
        gdalstr = fname
    return gdalstr


class gdal_open:

    def __init__(self, fname, mode=gdal.GA_ReadOnly, group=None, varn=None):
        """Open file with GDAL (tif or netCDF)

        Parameters
        ----------
        fname : str
            path to input file
            can be netCDF
        group : str
            name of group in netCDF or HDF file
        varn : str, optional
            name of variable in netCDF or HDF file
        """
        gdalstr = make_gdalstr(fname, group=group, varn=varn)
        self.file = gdal.Open(gdalstr, mode)
        if self.file is None:
            raise IOError(
                "Loading file '{}' with GDAL failed".format(gdalstr))

    def __enter__(self):
        return self.file

    def __exit__(self, type, value, traceback):
        self.file = None


def check_gdal_success(outfile, cmd):
    """Make sure GDAL command `cmd` succeeded in creating `outfile`"""
    if not os.path.isfile(outfile):
        raise RuntimeError('GDAL command \'{}\' did not produce the '
                'expected output file {}.'.format(cmd, outfile))


def cutline(intif, inshp, outtif, t_srs=None):
    """GDAL cutline with default parameters"""
    if t_srs:
        cmd = cmd_gdalwarp_reproject_cutline(intif, inshp, outtif, t_srs)
    else:
        cmd = cmd_gdalwarp_cutline(intif, inshp, outtif)
    run_cmd(cmd, outtif)


def warp(infile, outfile, **kwargs):
    """Run gdalwarp"""
    cmd = cmd_gdalwarp_reproject(infile, outfile=outfile, **kwargs)
    run_cmd(cmd, outfile)


def warp_reproject(infile, outfile, t_epsg=4326, r='near'):

    # Define target SRS
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(t_epsg)
    dst_wkt = dst_srs.ExportToWkt()

    # error threshold
    error_threshold = 0.125  # error threshold --> use same value as in gdalwarp

    # resampling
    if r == 'near':
        resampling = gdal.GRA_NearestNeighbour
    elif r == 'bilinear':
        resampling = gdal.GRA_Bilinear
    else:
        raise ValueError('Resampling `r` should be \'near\' or \'bilinear\'.')

    with gdal_open(infile) as src_ds:
        # Call AutoCreateWarpedVRT() to fetch default values for target raster dimensions and geotransform
        tmp_ds = gdal.AutoCreateWarpedVRT(src_ds,
                                          None,  # src_wkt
                                          dst_wkt,
                                          resampling,
                                          error_threshold)

    # Create the final warped raster
    if outfile.lower().endswith('.vrt'):
        driver = 'VRT'
    else:
        driver = 'GTiff'
    dst_ds = gdal.GetDriverByName(driver).CreateCopy(outfile, tmp_ds)
    dst_ds = None
    check_gdal_success(outfile, 'gdalwarp')


def buffer_extent(extent, geotransform):
    """Buffer extent by one pixel to the east and south (if necessary) to make
       sure that gdal_translate leaves no bit of AOI out when subsetting

    Parameters
    ----------
    Extent is minX, maxX, minY, maxY
    Geotransform is GDAL geotransform
    """
    if ',' in extent:
        e = extent.split(',')
    else:
        e = extent
    # From what I could figure out gdal_translate first shifts the extent north-west until UL corner aligns
    # with a UL corner of a pixel in the source layer. Then the BR corner is rounded to the BR corner of
    # the closest pixel. So if the north(west) pixel is closer then south(east) then the clipped layer will
    # not include all of the area specified in the extent.
    shiftWest = float(e[0]) % geotransform[1]
    if 0 < float(e[1]) % geotransform[1]/geotransform[1] < 0.5:
        east = float(e[1]) + geotransform[1]
    else:
        east = float(e[1])
    e[1] = str(east + shiftWest)
    shiftNorth = float(e[3]) % geotransform[5]
    if 0 < float(e[2]) % geotransform[5]/geotransform[5] < 0.5:
        south = float(e[2]) + geotransform[5]
    else:
        south = float(e[2])
    e[2] = str(south + shiftNorth)
    return e


def buildvrt(infiles, outfile, **kwargs):
    """GDAL build virtual raster

    Parameters
    ----------
    infiles : list of str
        paths to input files
    outfile : str
        path to output vrt
    kwargs : dict
        keyword arguments passed to
        cmd_gdalbuildvrt
    """
    if len(infiles) < 5:
        cmd = cmd_gdalbuildvrt(infiles=infiles, outfile=outfile)
        run_cmd(cmd, outfile)
    else:
        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            f.writelines('\n'.join(infiles))
            input_file_list = f.name
        try:
            cmd = cmd_gdalbuildvrt(input_file_list=input_file_list, outfile=outfile, **kwargs)
            run_cmd(cmd, outfile)
        finally:
            os.remove(input_file_list)


def translate(infile, outfile, **kwargs):
    """GDAL build virtual raster

    Parameters
    ----------
    infile : str
        path to input file
    outfile : str
        path to output file
    kwargs : dict
        keyword arguments passed to
        cmd_gdal_translate
    """
    cmd = cmd_gdal_translate(infile, outfile, **kwargs)
    run_cmd(cmd, outfile)


def get_target_extent(fname, varn):
    # Get extent and spatial resolution of the master data
    # and provide as extra parameters since resolution can
    # be different in X and Y dimensions
    with gdal_open(fname, varn) as masterFile:
        gt = masterFile.GetGeoTransform()
        resX = gt[1]
        resY = gt[5]
        ext = [gt[0], gt[3]+gt[5]*masterFile.RasterYSize, gt[0] + gt[1]*masterFile.RasterXSize, gt[3]]
    return ("-tr {resX} {resY} -te {ext[0]} {ext[1]} "
            "{ext[2]} {ext[3]}".format(resX=resX, resY=resY, ext=ext))
