import os
import datetime
import tempfile
import shutil
from zipfile import BadZipfile
import logging

import gdal_utils as gu

from external import copernicus_download as codo

logger = logging.getLogger('download_swi')


def download(outdir, progress, username, password,
        product='SWI',
        startdate='', enddate='', extent='',
        delete_downloaded=False):

    # PARSE
    valid_products = ['SWI', 'SWI10']
    if product not in valid_products:
        raise ValueError('Valid products are {}.'.format(valid_products))
    varn = 'SWI_100'

    download_dir = os.path.join(outdir, 'download')
    try:
        os.makedirs(download_dir)
    except OSError:
        pass

    # Parse start and end dates
    dt_days = codo.config.timestep_days[product]
    dt = datetime.timedelta(days=dt_days)
    startdate = datetime.datetime.strptime(startdate, '%Y%m%d')
    enddate = datetime.datetime.strptime(enddate, '%Y%m%d')

    extent_kw = {}
    global_extent = '-180,180,-90,90'
    if extent not in ['', None, global_extent]:
        xmin, xmax, ymin, ymax = map(float, extent.split(','))
        extent_kw = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    else:
        extent = global_extent

    output_files = []

    # DOWNLOAD
    download_attempts = 3
    date = startdate
    while date <= enddate+dt:

        # If using a composite product then the date should be at the
        # beginning of the compositing period
        if dt_days > 1:
            day = date.day
            day = int((day-1)/dt_days)*dt_days+1
            download_date = date.replace(day=day)
        else:
            download_date = date

        # Download the data
        query_url = codo.build_url(product=product,
                year=download_date.year, month=download_date.month, day=download_date.day,
                extent=extent_kw)
        logger.info('Query URL is \'{}\'.'.format(query_url))

        download_tempdir = tempfile.mkdtemp(prefix='download_', dir=download_dir)

        tempfiles = []
        h5files = []
        try:
            h5fname = None
            download_OK = False
            download_n = 1
            while not download_OK and download_n <= download_attempts:
                logger.info('Starting download attempt nr {}'.format(download_n))

                downloaded_files = codo.download_data(query_url,
                        username=username, password=password,
                        download_dir=download_tempdir)

                if len(downloaded_files) != 1:
                    logger.info('Expecting the script to download exactly '
                            'one file per date. Got {}: {}.'.format(len(downloaded_files), downloaded_files))
                    download_n += 1
                    continue

                zipFilename = downloaded_files[0]
                logger.info("Extracting archive...")
                logger.info(zipFilename)
                try:
                    h5fname = codo.extract_h5(zipFilename, download_dir)
                    h5files.append(h5fname)
                except BadZipfile:
                    logger.info("Warning: Could not successfully extract data for date {}.".format(date))
                    download_n += 1
                    try:
                        os.remove(h5fname)
                    except OSError:
                        pass
                    continue

                download_OK = True

            if not download_OK:
                # Start downloading the next file
                logger.info('WARNING: Failed to download data for date {}. Skipping.'.format(date))
                date += dt
                continue

            # read data from h5
            data = codo.read_h5(h5fname, group=product, varn=varn)
            tempnc = h5fname.replace('.h5', '.nc')
            data.to_netcdf(tempnc)
            tempfiles.append(tempnc)

            # translate to GeoTIFF
            outfile = os.path.join(outdir, '{}_{:%Y%m%d}.tif'.format(product, date))
            logger.info(outfile)

            gdalstr = gu.make_gdalstr(tempnc, varn=varn)
            with gu.gdal_open(gdalstr) as ds:
                gt = ds.GetGeoTransform()
            buffered_extent = gu.buffer_extent(extent, gt)
            gu.translate(gdalstr, outfile, extent=buffered_extent, extra=['-a_srs', 'EPSG:4326'])
            output_files.append(outfile)

        finally:
            if delete_downloaded:
                # also remove h5files
                tempfiles += h5files
            shutil.rmtree(download_tempdir)
            for filename in tempfiles:
                try:
                    os.remove(filename)
                except:
                    pass

        date += dt

    return output_files
