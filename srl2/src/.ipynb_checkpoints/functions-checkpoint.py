import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import os
import copy
from astropy.wcs import WCS
from astroquery.vizier import Vizier
import aplpy
import matplotlib.pyplot as plt
from astroquery.hips2fits import hips2fits
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import astropy.units as u
from astropy.coordinates import Longitude, Latitude, Angle
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import pandas as pd
import aplpy
#from astroquery.hips import hips2fits

def load_fits_data(file_path):
    """Load FITS data and convert to a pandas DataFrame.
     Parameters:
        file_path (str): Path to the FITS file.
    
    Returns:
        pandas.DataFrame: Data extracted from the FITS file."""
    with fits.open(file_path) as hdul:
        data = hdul[1].data
        return Table(data).to_pandas()

def filter_kids_data(df, ra_min, ra_max, dec_min, dec_max):
    """
    Filter KiDS data based on RA and DEC ranges.
    Parameters:
        df (pandas.DataFrame): DataFrame containing KiDS data.
        ra_min (float): Minimum RA value.
        ra_max (float): Maximum RA value.
        dec_min (float): Minimum DEC value.
        dec_max (float): Maximum DEC value.
    
    Returns:
        pandas.DataFrame: Filtered KiDS data within the given RA and DEC range
    """
    fil_df = df[(df['RAJ2000'] >= ra_min) & (df['RAJ2000'] <= ra_max) & 
                (df['DECJ2000'] >= dec_min) & (df['DECJ2000'] <= dec_max)]
    return fil_df.reset_index(drop=True)


def find_closest_matches(fil_df, df_mk):
    """Find closest matches between KiDS and MeerKAT sources based on angular separation.
    Parameters:
        fil_df (pandas.DataFrame): Filtered KiDS data containing RA and DEC.
        df_mk (pandas.DataFrame): MeerKAT data containing RA and DEC.
    
    Returns:
        pandas.DataFrame: DataFrame containing the closest matched sources and their separations.
    
    """
    closest_matches = []
    
    # Convert coordinates to SkyCoord
    coords1 = SkyCoord(ra=fil_df['RAJ2000'].values * u.deg, dec=fil_df['DECJ2000'].values * u.deg)
    coords2 = SkyCoord(ra=df_mk['RA'].values * u.deg, dec=df_mk['DEC'].values * u.deg)

    for idx2, mk_coord in enumerate(coords2):
        sep = mk_coord.separation(coords1)
        closest_idx = np.argmin(sep)
        closest_separation = sep[closest_idx].deg

        kids_info = fil_df.iloc[closest_idx].to_dict()
        meerkat_info = df_mk.iloc[idx2].to_dict()
        match_info = {**kids_info, **meerkat_info}
        match_info['MK_RA'] = df_mk['RA'].iloc[idx2]
        match_info['MK_DEC'] = df_mk['DEC'].iloc[idx2]
        match_info['Separation'] = closest_separation

        closest_matches.append(match_info)
    return pd.DataFrame(closest_matches)


def filter_and_save_matches(df, separation_threshold, output_file):
    """
    Filter sources with separation below threshold and save to a file.
    Parameters:
        df (pandas.DataFrame): DataFrame containing matched sources with separation values.
        separation_threshold (float): Maximum allowed separation (in degrees) between matches.
        output_file (str): File path to save the filtered results.
    
    Returns:
        pandas.DataFrame: DataFrame containing sources that did not meet the threshold.
    """
    filtered_df = df[df['Separation'] < separation_threshold]
    print(f'Amount of sources with a separation of {separation_threshold * 3600:.1f} arcsec:', len(filtered_df))

    df = df[~df.index.isin(filtered_df.index)].reset_index(drop=True)
    filtered_df = filtered_df.reset_index(drop=True)

    filtered_df.to_csv(output_file, sep=',', index=False)
    return df
#################################################################################
def load_coordinates(coord_file):
    """
    Load source coordinates from a text file.
    Parameters:
        coord_file (str): Path to the text file containing coordinates.
    
    Returns:
        tuple: ID, RA, Dec arrays, and SkyCoord object.
    """
    ID, RA, Dec = np.loadtxt(coord_file, unpack=True, usecols=[0, 1, 2], delimiter=',', skiprows=1, dtype=[('ID', '|S50'), ('MK_RA', 'f8'), ('MK_DEC', 'f8')])
    c = SkyCoord(RA, Dec, unit="deg")
    return ID, RA, Dec, c

def get_pixel_coordinates(c, wcs_header):
    """
    Convert world coordinates to pixel coordinates.
     Parameters:
        c (SkyCoord): SkyCoord object containing world coordinates.
        wcs_header (dict): WCS header from a FITS file.
    
    Returns:
        tuple: Arrays of pixel x and y coordinates.
    """
    w = WCS(wcs_header, naxis=2)
    x, y = w.world_to_pixel(c)
    pixel_x = np.rint(x).astype(int)
    pixel_y = np.rint(y).astype(int)
    return pixel_x, pixel_y
fits_file = 'D01-05_LOC22_im-di2_smallFacet.deeper.DI.int.restored.fits'
def cutout(file, xc, yc, name, xw=100, yw=100, units='pixels', JPEGFolder="WESTTEST", clobber=True):
    """
    Create a cutout from a FITS file at specified coordinates.
    
    """
    if isinstance(file, str):
        file = fits.open(file)
        opened = True
    elif isinstance(file, fits.HDUList):
        opened = False
    else:
        raise Exception("cutout: Input file is wrong type (string or HDUList are acceptable).")

    head = file[0].header
    try:
        cd1 = head['CDELT1']
        cd2 = head['CDELT2']
    except:
        try:
            cd1 = head['CD1_1']
            cd2 = head['CD2_2']
        except:
            raise Exception("No CD or CDELT keywords in header")

    if not os.path.exists(JPEGFolder):
        os.makedirs(JPEGFolder)
    
    for n, x, y in zip(name, xc, yc):
        naxis1 = head['NAXIS1']
        naxis2 = head['NAXIS2']       
        crpix1 = head['CRPIX1']
        crpix2 = head['CRPIX2']

        if units == 'pixels':
            xmin, xmax = max(0, x - xw), min(naxis1, x + xw)
            ymin, ymax = max(0, y - yw), min(naxis2, y + yw)
        elif units == 'wcs':
            xmin, xmax = max(0, x - xw / abs(cd1)), min(naxis1, x + xw / abs(cd1))
            ymin, ymax = max(0, y - yw / abs(cd2)), min(naxis2, y + yw / abs(cd2))
        else:
            raise Exception(f"Can't use units {units}.")
        
        if xmax < 0 or ymax < 0:
            raise Exception("Coordinate is outside of map.")

        img = file[0].data[0, 0, ymin:ymax, xmin:xmax]
        crpix1 -= xmin
        crpix2 -= ymin
        naxis1 = img.shape[1]
        naxis2 = img.shape[0]
    
        newhead = copy.deepcopy(head)
        newhead['CRPIX1'] = crpix1
        newhead['CRPIX2'] = crpix2
        newhead['NAXIS1'] = naxis1
        newhead['NAXIS2'] = naxis2
        newfile = fits.PrimaryHDU(data=img, header=newhead)
        
        outfile = os.path.join(JPEGFolder, n.decode('UTF-8') + '.fits')
        newfile.writeto(outfile, overwrite=clobber)

#cutout1 = cutout(file=fits_file,xc=pixel_x,yc=pixel_y,xw=1000,yw=1000,units='pixels',name=ID,clobber=True)


def get_background_variance(data, sigma_clip=5.0, tolerance=0.01):
    """Compute the variance by iteratively removing outliers greater than a given sigma
    until the mean changes by no more than tolerance.
    Parameters:
        data (numpy.ndarray): Input data array.
        sigma_clip (float, optional): Sigma threshold for clipping (default: 5.0).
        tolerance (float, optional): Convergence tolerance (default: 0.01).
    
    Returns:
        float: Computed background variance.
    """
    diff = 1
    mean = np.nanmean(data)
    data_clip = data
    while diff > tolerance:
        data_clip = data_clip[np.abs(data_clip)<mean+sigma_clip*np.nanstd(data_clip)]
        newmean = np.nanmean(data_clip)
        diff = np.abs(mean-newmean)/(mean+newmean)
        mean = newmean
    return np.nanvar(data_clip)

def download_cutout(target_RA, target_DEC, hips, fov=0.05, width=256, height=256):
    """Download the FITS cutout for a given RA/DEC from the HIPS catalog.
    
    Parameters:
        target_RA (float): Right Ascension of the target (degrees).
        target_DEC (float): Declination of the target (degrees).
        hips (str): HiPS survey name.
        fov (float, optional): Field of view in degrees (default: 0.05).
        width (int, optional): Image width in pixels (default: 256).
        height (int, optional): Image height in pixels (default: 256).
    
    Returns:
        FITS file: FITS image cutout resul"""
    result = hips2fits.query(
        hips=hips,
        width=width,
        height=height,
        ra=Longitude(target_RA),
        dec=Latitude(target_DEC),
        fov=Angle(fov * u.deg),
        projection="SIN",
        get_query_payload=False,
        format='fits',
        min_cut=0,
        max_cut=100,
        cmap='viridis'
    )
    return result
def clean_fits_header(fitsfile):
    """Removes unwanted header keys from FITS file."""
    HF_fits = fits.open(fitsfile, mode='update')
    HF_header = HF_fits[0].header
    
    remove_list = ['WCSAXES', 'CTYPE3', 'CRVAL3', 'CDELT3', 'CRPIX3', 'CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4']
    for key in remove_list:
        if key in HF_header:
            del  HF_header[key]
            #HF_header.remove(key)
    HF_fits.close()

def compute_contour_levels(fitsfile):
    """Computes contour levels based on data variance.
    
      Parameters:
        fitsfile (str): Path to the FITS file.
    
    Returns:
        np.ndarray: Contour levels for plotting."""
    fitsimage = fits.open(fitsfile)
    data = fitsimage[0].data
    sigma = np.sqrt(get_background_variance(data.flatten()))
    HF_levels = np.logspace(np.log10(3.0 * sigma), np.log10(45 * sigma), num=7)
    return HF_levels

def create_overlay_image(target_name, result, fitsfile, xG, yG, HF_levels, index):
    """Creates and saves an overlay image with contours."""
    plt.rcParams.update({'font.size': 15})
    fig = plt.figure(figsize=(8, 8))
    
    f = aplpy.FITSFigure(result, slices=[0], figure=fig)
    f.show_grayscale(invert= True)
    f.show_contour(fitsfile, slices=[0,0], levels=HF_levels, colors='b', overlap=True, smooth=1)
    f.show_markers(xG, yG, marker='s', edgecolor='r', s=50, linewidths=1)
    
    plt.title(f"{target_name}[{index}]")
    #Save the figure
    SUFFIX = '-MeerKLASS_KiDS'       
    fig.savefig(f'image_overlay/{target_name}-MeerKLASS_KiDS.png')