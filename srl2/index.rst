.. Multiwavelength cros-matching with MeerKLASS documentation master file, created by
   sphinx-quickstart on Tue Mar 25 17:18:17 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Multiwavelength cros-matching with MeerKLASS documentation
==========================================================
# MeerKAT and KiDS DR4 Cross-Matching Pipeline Documentation

## Overview
This project cross-matches MeerKAT radio sources with KiDS DR4 optical sources using astrometric positions (RA, Dec). The goal is to identify potential counterparts and analyze their properties across different wavebands. This is done using FITS files from both catalogs and dataframes with the physical infomation about the sources.

## Requirements

To run this pipeline, you need the following dependencies:

- Python 3.8
- `pandas`
- `astropy`
- `numpy`
- `matplotlib`
- `astroquery`
- `aplpy`
- `configparser`
- `astropy.coordinates`
- `os`

You can install these packages using:

pip install pandas astropy numpy matplotlib astroquery aplpy, configparser


## Data Sources

The pipeline requires the following input data:

- `KiDS_DR4_brightsample_LePhare.fits`: Contains redshift, colour, star formation rate (SFR), and other optical properties of galaxies from the KiDS DR4 survey.
- `D01-05_LOC22_im-di2_smallFacet.deeper.DI.int.restored.pybdsf.srl.fits`: Contains MeerKAT radio source properties such as flux and coordinates.
- `D01-05_LOC22_im-di2_smallFacet.deeper.DI.int.restored.fits`: MeerKAT radio mosaic image with source positions.
- `closest_matches_df1.txt`: Generated file listing matched sources with relevant properties.
- HiPS to FITS input parameters:

The HiPS-to-FITS service is used to extract FITS images from a HiPS (Hierarchical Progressive Surveys) dataset. The following parameters are used as input for the query:

hips (str): Name of the HiPS survey from which the FITS image is requested.

width, height (int): Dimensions of the output FITS image in pixels.

ra, dec (float, degrees): Right Ascension (RA) and Declination (Dec) of the target sky position, converted into Longitude and Latitude objects.

fov (float, degrees): Field of view (FoV) in degrees, converted into an astropy.units.Angle object.

projection (str, default: "SIN"): Sky projection method used for the output image.

get_query_payload (bool, default: False): If True, the function returns the query payload instead of executing the request.

format (str, default: "fits"): Specifies the output format of the image (FITS format is used).

min_cut, max_cut (int, optional): Intensity scaling limits for the image.

cmap (str, optional, default: "viridis"): Colormap used for visualization.

This query extracts a FITS image centered on the given sky coordinates, allowing for further analysis of the selected HiPS dataset.

## Processing Steps/Workflow

### 1. Loading Data
- The KiDS DR4 and MeerKAT data are read from their respective FITS files.
- The coner coordiantes are extracted to use the relevenat area of the KiDS catalog for  filtring and crossmatcing. 

### 2. Filtering KiDS Data
- The dataset is filtered based on specific RA and Dec ranges.
- The subset of sources is extracted and reset for indexing consistency.

### 3. Matching Sources
- The script calculates the angular separation between MeerKAT and KiDS sources.
- The closest KiDS counterpart is identified for each MeerKAT source.
- A DataFrame (`closest_matches_df`) is created containing matched sources and their separation distances.

### 4. Filtering Based on Separation
- Matched sources are categorized based on angular separation:
  - Sources with separation < 3 arcsec which can be configured in the scirpt.
- Sources are progressively removed from the dataset after each filtering step.

### 5. Generating FITS Cutouts
- The script extracts small FITS cutouts around matched sources,using the coordinates from the same catalog
- The cutout function ensures that sub-images are centered on source positions.
- The cutouts are saved for further processing.

### 6. Generating Contour plots and overplotting the plot on top of the 'optical image'
- A set of optical plot is queried using Vizier using the dataframe and the coordinates of the optical catalog that was generated in the crossmatching stage
- Once the image is generated the contour levels are calculated and they are overplotted onto the optical image to show the radio emissions captured using the MeerKLASS survey 

## Output Files
- `closest_matches_df1.txt`: Filtered match lists with the infomation about the crossmatched source, containing infomation offered by both catalogs.
- FITS cutout files for selected sources from the radio images which can be found in the 'MeeKlASS_cutouts' folder.
- Cutouts images of the sources from from the optical survey centred around the specific sources with the 'radio contours' that represent the radio emission caputured in MeerKLASS survey, super-plotted onto the source. The imagescan be found the directory named './image_overlay'.

## Usage

Make sure the required FITS and CSV files are in the working directory.

## Notes
- Modify RA/Dec filtering parameters if needed to change the region of interest.
- Adjust the cutout size in the `cutout()` function if necessary.
- Ensure proper installation of dependencies to avoid import errors.
- Ensure the values in the 'parameter.ini' file are called correctly

## Troubleshooting
- **ImportError**: Check that all required Python packages are installed.
- **FileNotFoundError**: Ensure that input FITS and CSV files are in the correct directory.
- **Incorrect Matches**: Verify the RA/Dec columns and filtering conditions to ensure accurate cross-matching.

For further assistance, refer to the documentation of `astropy`, `astroquery`, and `aplpy`.


Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

