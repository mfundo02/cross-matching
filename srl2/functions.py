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
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
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

def filter_kids_by_fits(df, meerklass_image):
    """
    Extracts the RA and DEC range from a FITS image and filters KiDS data within that range.

    Parameters:
        df (pandas.DataFrame): DataFrame containing KiDS data with 'RAJ2000' and 'DECJ2000' columns.
        fits_file (str): Path to the FITS file.

    Returns:
        pandas.DataFrame: Filtered KiDS data within the RA and DEC range of the FITS image.
    """
    # Open the FITS file and extract WCS information
    with fits.open(meerklass_image) as hdul:
        header = hdul[0].header
        wcs = WCS(header, naxis=2)  # Use only the spatial axes

    # Get image dimensions
    ny, nx = header['NAXIS2'], header['NAXIS1']

    # Define the corner pixel coordinates
    x_pixels = np.array([0, 0, nx-1, nx-1])  # X-coordinates
    y_pixels = np.array([0, ny-1, 0, ny-1])  # Y-coordinates

    # Convert pixel coordinates to world coordinates
    ra_values, dec_values = wcs.pixel_to_world_values(x_pixels, y_pixels)

    # Get min and max RA and Dec
    ra_min, ra_max = np.min(ra_values), np.max(ra_values)
    dec_min, dec_max = np.min(dec_values), np.max(dec_values)

    # Filter the data frame based on the extracted RA and DEC range
    filtered_df = df[
        (df['RAJ2000'] >= ra_min) & (df['RAJ2000'] <= ra_max) &
        (df['DECJ2000'] >= dec_min) & (df['DECJ2000'] <= dec_max)
    ].reset_index(drop=True)

    return filtered_df


def find_and_filter_matches(fil_df, df_mk , separation_threshold):
    """
    Find closest matches between KiDS and MeerKAT sources based on angular separation,
    filter sources with separation below threshold, and save to a file.
    
    Parameters:
        fil_df (pandas.DataFrame): Filtered KiDS data containing RA and DEC.
        df_mk (pandas.DataFrame): MeerKAT data containing RA and DEC.
        separation_threshold (float): Maximum allowed separation (in degrees) between matches.
        output_file (str): File path to save the filtered results.
    
    Returns:
        pandas.DataFrame: DataFrame containing sources that did not meet the threshold.
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
    
    matched_df = pd.DataFrame(closest_matches)
    # Filter matches based on separation threshold
    filtered_df = matched_df[matched_df['Separation'] < separation_threshold/3600]
    print(f'Amount of sources with a separation of {separation_threshold * 3600:.1f} arcsec: {len(filtered_df)} out of {len(matched_df)}')

    filtered_df = filtered_df.reset_index(drop=True)
    
    # Save filtered matches to a file
    filtered_df.to_csv('closest_matches_df1.txt', sep=',', index=False)
    
    return filtered_df


def load_and_convert_coordinates(header, fits_file):
    """
    Load source coordinates from a text file and convert to pixel coordinates.
    
    Parameters:
        coord_file (str): Path to the text file containing coordinates.
        wcs_header (dict): WCS header from a FITS file.
    
    Returns:
        tuple: ID, RA, Dec arrays, SkyCoord object, and pixel x, y coordinates.
    """
    
    # Load coordinate data
    ID, RA, Dec = np.loadtxt('closest_matches_df1.txt', unpack=True, usecols=[0, 1, 2], delimiter=',', skiprows=1, 
                             dtype=[('ID', '|S50'), ('MK_RA', 'f8'), ('MK_DEC', 'f8')])
    
    # Create SkyCoord object
    c = SkyCoord(RA, Dec, unit="deg")
    f = fits.open(fits_file)
    # Convert to pixel coordinates
    #w = WCS(wcs_header, naxis=2)
    w = WCS(f[0].header, naxis=2)
    x, y = w.world_to_pixel(c)
    pixel_x = np.rint(x).astype(int)
    pixel_y = np.rint(y).astype(int)
    
    return ID, RA, Dec, c, pixel_x, pixel_y
    
#fits_file = 'D01-05_LOC22_im-di2_smallFacet.deeper.DI.int.restored.fits'

def cutout(fits_file, xc, yc, name, xw=100, yw=100, units='pixels', JPEGFolder="MeerKLASS_fits_cutouts", clobber=True):
    """
    Create a cutout from a FITS file at specified coordinates.
    
    """
    if isinstance(fits_file, str):
        file = fits.open(fits_file)
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

def create_dataframe(config):
    # File inputs and relevent qualities
    kids_fits = config['kiDS']['fits_data']
    meerkat_fits = config['meerKAT']['fits_data']
    meerklass_image = config['FITSFile']['input_fits']
    separation_threshold = float(config['General']['separation_threshold'])
    # Load data
    kids_df = load_fits_data(kids_fits)
    df_mk = load_fits_data(meerkat_fits)
    # Filter KiDS data
    fil_df = filter_kids_by_fits(kids_df, meerklass_image)
    closest_matches_df1 = find_and_filter_matches(fil_df, df_mk,separation_threshold)
    return closest_matches_df1

def process_fits_cutout(config):
    """
    #Process the FITS cutouts from the configuration and save them.
    """
    # Load parameters from the .ini file
    #fits_file = 'el_withsmear_nr0.deeper.DI.int.restored.fits'
    fits_file = config['FITSFile']['input_fits']#the mosaic image
    output_folder = config['FITSFile']['output_folder']
    xw = int(config['Cutout']['xw'])
    yw = int(config['Cutout']['yw'])
    units = config['Cutout']['units']
    clobber = config.getboolean('Cutout', 'clobber')

    # Load FITS data and coordinates
    fits_data = fits.getdata(fits_file, 0)
    header = fits.getheader(fits_file)
    image_data = fits_data[0][0]
    #fits_data, header, image_data = load_fits_data(fits_file)
    
     #Load and Convert the coordinates from the MeerKLASS catalog into pixels.
    ID, RA, Dec, c, pixel_x, pixel_y = load_and_convert_coordinates(header, fits_file)
    
    # Create cutout
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    cutout(fits_file, xc=pixel_x, yc=pixel_y, name=ID, xw=xw, yw=yw, units=units, JPEGFolder="MeerKLASS_fits_cutouts", clobber=clobber)

def process_hips_data(config):
    """
    #Process the HIPS catalog data for the given targets and save results.
    """
    hips_catalog = config['General']['hips_catalog']
    output_dir = config['General']['output_dir']
    input_file = config['General']['input_file']
    fov = float(config['General']['fov'])
    width = int(config['General']['width'])
    height = int(config['General']['height'])
    output_folder = config['FITSFile']['output_folder']

    mkcat = np.genfromtxt('closest_matches_df1.txt', 
                      delimiter=',',      # Adjust the delimiter if necessary
                      dtype=None,        # Allow mixed data types
                      encoding=None,     # Use 'utf-8' or None depending on your data
                      names=True)        # Load column names as field names
    # Create directories if they don't exist
    os.makedirs('cutout_mk', exist_ok=True)
    os.makedirs('image_overlay', exist_ok=True)
    import shutil
    shutil.copytree('./MeerKLASS_fits_cutouts', 'cutout_mk', dirs_exist_ok=True)

    
    # Process each target
    for ii in range(mkcat.shape[0]):
        target_name = str(mkcat['ID'][ii])
        target_RA = mkcat['MK_RA'][ii] * u.deg
        target_DEC = mkcat['MK_DEC'][ii] * u.deg
        print(f"Processing {target_name} at RA: {target_RA}, DEC: {target_DEC}")

        #hips = 'CDS/P/KiDS/DR5/color-gri'
        # Download the optical cutout from the HIPS catalog
        result = download_cutout(target_RA, target_DEC, hips_catalog, fov=fov, width=width, height=height)
        #extract the coordinates of the closest matches
        df_clo = pd.read_csv('closest_matches_df1.txt', delimiter=',')
        # Extract the RA and DEC for the closest matches
        xG = df_clo['RAJ2000'][ii]*u.deg
        yG = df_clo['DECJ2000'][ii]*u.deg

        fitsfile = f'cutout_mk/{target_name}.fits'
        clean_fits_header(fitsfile)
        
        HF_levels = compute_contour_levels(fitsfile)
        create_overlay_image(target_name, result, fitsfile, xG, yG, HF_levels, ii)
        