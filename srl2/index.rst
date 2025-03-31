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

pip install pandas astropy numpy matplotlib astroquery aplpy, configparser, os etc


## Data Sources
The pipeline requires the following input data:
- **KiDS DR4 Bright Sample FITS file, contains the physical infomation collected in the optical survey such as the redshift, coulor, SFR etc.** (`KiDS_DR4_brightsample_LePhare.fits`)
- **MeerKLASS radio source FITS file, contains the physical infomation that was collected in the radio survey such the the radio flux, astrometric positio etc** (`D01-05_LOC22_im-di2_smallFacet.deeper.DI.int.restored.pybdsf.srl.fits`)
- **MeerKLASS radio source mosaic image FITS file, contains the astrometric coordinates and the actual image of the sources captured in the survey** ('D01-05_LOC22_im-di2_smallFacet.deeper.DI.int.restored.fits') 
- **Text file with coordinates for FITS cutouts, which will be created in the script after the crossmatching process** (`closest_matches_df1.txt`)

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
- `closest_matches_df1.txt`: Filtered match lists with the infomation about the crossmatched source, containing infomation offered by both catalogs
- FITS cutout files for selected sources.

## Usage

Make sure the required FITS and CSV files are in the working directory.

## Notes
- Modify RA/Dec filtering parameters if needed to change the region of interest.
- Adjust the cutout size in the `cutout()` function if necessary.
- Ensure proper installation of dependencies to avoid import errors.

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

