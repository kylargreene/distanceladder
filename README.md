distanceladder, putting the local observations into MontePython!
================================================================

Main Developer: Kylar Greene <kygreene@unm.edu>

To make this tool as understandable and usable as possible by the broader cosmology community, this readme will go into extensive detail on how the code functions and how it is implemented into MontePython. It will be divided into three sections; installation, file and function definitions, and computation.

A flowchart is provided below describing the process flow of the package. The light pink area indicates processes involving the local distance ladder. The light blue background indicates processes involving cosmology. Boxes found on the border of the two backgrounds depend on both. The pentagons indicate user inputs. The circles indicate data sets used in the function. The soft rectangles indicate processes done by distanceladder.

![distanceladder Flowchart](https://i.imgur.com/btRsdwP.png)

INSTALLATION
------------

The distanceladder package exists in two parts; the python code, which MontePython calls and the datasets, which distanceladder calls. The python code is contained in a folder "likelihood", and within it has the two files "__init__.py" and "distanceladder.data". "__init__.py" is the code MontePython calls to calculate the chi-squared value of the likelihood. "distancealdder.data" contains the user set options which control the calibration scheme used by __init__.py. 

The datasets are in a folder "data", and within it has nine text files which will be enumerated in the next section.

To install the distanceladder package, copy the likelihood folder into your montepython/likelihoods/ directory and __RENAME THE FOLDER 'DISTANCELADDER'__.
Then, copy the data folder into your /data/ directory and __RENAME THE FOLDER 'DISTANCELADDER'__.

USAGE
-----

Once installed, the distanceladder likelihood is easy to use. The calibration options are chosen in the /likelihoods/distanceladder/distanceladder.data file, with options listed here.

* "distanceladder.data_directory = data.path['data']" defines where the tabulated data files are located and should not be changed if installed in the folder described above.

* "distanceladder.calib = " is where the user chooses which calibration scheme to use. Options include 'cepheid', 'TRGB', 'concordant'.

* "distanceladder.anchor = " is where the user chooses which anchors to use. TRGB and Concordant assume LMC only. Cepheid calibration is sensitive to this choice with options, ['N4258'], ['LMC'], or ['N4258','LMC'].

* "distanceladder.supernova_data = " is a future flag to allow the user to select between the pantheon or CSP datasets for Hubble flow supernova. However, only 'panth' is supported in this release.

FILE AND FUNCTION DEFINITIONS
-----------------------------

DATA FILES
----------

* "anchor_data.txt" contains the distance and name information to the local anchors used to calibrate the Cepheid period to luminosity relationship. It contains three columns; (1) Anchor Name (2) Distance to Anchor in Parsecs (3) Error in Distance to Anchor in Parsecs.

* "SHOES_cepheid_data.txt" contains data from https://arxiv.org/abs/1607.08658 pertaining to the Cepheid populations in both the anchor measurements, and supernova host galaxies. It contains ten columns; (1) Host Name (2) Cepheid ID (3) Period in Days (4) V band magnitude (5) error in V band magnitude (6) I band magnitude (7) error in I band magnitude (8) NIR band magnitude (9) error in NIR band magnitude (10) metalicity.

* "LMC_ceph_data.txt" contains data from https://arxiv.org/abs/1510.03682 pertaining to the Cepheid population in the LMC. A cut has been made to select only type II Cepheids as Riess and previous studies have done. It contains twentythree columnes; (1) Cepheid ID (2) Fundamental Resonation Mode (3) Period in log(days) (4) Weisenheit j-h mag (5) Weisenheit j-k mag (6) Weisenheit h-k mag (7) Weisenheit v-j (8) Weisenheit v-h (9) Weisenheit v-k (10) Weisenheit i-j mag (11) Weisenheit H magnitude (12-23) are the errors associated with the above magnitudes.

* "SHOES_sn_data.txt" contains data from table 5 of https://arxiv.org/pdf/1604.01424.pdf pertaining to the local supernovae observed in host galaxies which also contain observable Cepheids. It contains four columns; (1) Host Name (2) Supernova ID (3) B magnitude (4) error in B Magnitude.

* "b_TRGB_sn_data.txt" contains data from table 3 of https://arxiv.org/abs/1907.05922 pertaining to the distance measurements and magnitude measurements of supernovae in host galaxies calibrated by TRGB method. It contains six columns; (1) Host (2) Supernova ID (3) distance modulus (4) error in distance modulus (5) B filter magnitude *NOT B PRIME* (6) error in B filter magnitude.

* "combined_localsn_data.txt" contains data from table 5 and table 3 of Riess 2016 and Freedman 2019 in which the distance measurements agree to one sigma. IT contains six columns; (1) Host (2) Supernova ID (3) distance modulus (4) error in distance modulus (5) B filter magnitude *NOT B PRIME* (6) error in B filter magnitude.

* "panth_sn_data.txt" contains data from https://arxiv.org/abs/2005.07707 and https://arxiv.org/abs/1710.00845 pertaining to the Hubble flow supernova. It contains six columns; (1) Supernova ID (2) the Heliocentric redshift (3) the CMB frame redshift (4) the error in redshift measurement (5) B filter magnitude (6) error in B filter magnitude

functions and classes in order of appearance
--------------------------------------------

* "lin_reg()" is a simple linear regression used as a parameter estimator for the york fit. Its inputs are an x and y array. Its outputs are an intercept (B_0) and the slope (B_1).

* "weighted_fixed_slope()" is a more detailed linear regression that uses the error terms as a weighting scheme and allows the user to define a fixed slope. Its inputs are an x and y array, a slope (m), the error in the slope (dm), and the error in the x and y array (dx) (dy). Its outputs are the intercept (B0), the slope (B1), the error in the intercept (dB0), and the error in the slope (dB1).

* "york_fit()" is a more robust linear regression that seeks to minimize the fit error in both the x and y-direction. Its inputs are and x and y array, the error in x (sigma_x), the error in y (sigma_y), the correlation (r), the tolerance of the fit to know when to end the iterations (tol), and the number of maximum steps allowed (n_max). Its outputs are the intercept (B_0), the slope (B_1), the error in the intercept (sigma_a_new), the error in the slope (sigma_b_new), the history of the intercept as the iterations progressed (b_hist), the simple linear regression intercept (B_0_simple), and the simple linear regression slope (B_1_simple).

* "Anchor" class contains the information from the geometric distance anchor measurements. The __init__() function inputs a name, distance measurement in parsecs, and error in distance measurement in parsecs. It computes the distance modulus and error in distance modulus. The Compute_Abs_ceph_mag() function inputs are a period in days, error in period in days, Weisenheit H band magnitude, and error in Weisenheit H band magnitude. It computes the absolute Cepheid magnitude Mceph and error associated dMceph.

* "SHOES_ceph_data" class contains the information on the Cepheids found in both the anchor measurements and supernova host galaxies. The __init__() function inputs are the host name, the cepheid ID, the period in days, the V band mag, the error in V band mag, the I band mag, the error in I band mag, the NIR band mag, the error in NIR band mag, and the metalicity information. It computes the Weisenheit H band magnitude assuming R = 0.386 as done by Riess, and the associated error. The proto_Compute_mu() inputs are the Weisenheit H band mag, the Weisenheit H band error mag, the period in days, the absolute magnitude of Cepheids, the error in the absolute magnitude of Cepheids, the slope of the P-L relationship for Cepheids, and the error in slope for the P-L relationship for Cepheids. It computes the distance modulus and associated error to Cepheids.

* "LMC_ceph_data" class contains the information on the Cepheids found in the LMC. A different class was required as the V, I, and NIR band information was not readily available, but the direct Weisenheit magnitude was. The __init__() input requires host name, Cepheid ID, period in days, Weisenheit H band magnitude, and error in Weisenheit H band magnitude. The proto_compute_mu is identical to the case in SHOES_ceph_data.

* "TRGB_sn_data" class contains the supernova in galaxies with observable TRGB features. The __init__() function inputs are a host name, supernova ID, distance modulus to SN, error in distance modulus, B band magnitude of supernova, and error in B band mag.

* "Local_SN_data" class contains the supernova in galaxies with observable Cepheid calibrators. The __init__() inputs are host name, supernova ID, B band mag, and error in B band mag. The Compute_abs_sn_mag() inputs are B band magnitude, error in B band magnitude, distance modulus, and error in distance modulus. It computes the absolute magnitude of supernova and associated error through linear regression with fixed slope. The slope is fixed to one with zero error as we define the relationship between magnitude and distance modulus to be so.

* "Hubble_SN_data" class contains the Hubble flow supernova. The __init__() function inputs are supernova ID, B band mag, error in B band mag, heliocentric redshift, CMB frame redshift, and error in redshift. The Compute_hubble_mu() inputs are B band mag, error in B band mag, absolute supernova mag, and error in absolute supernova mag. It computes the distance modulus and associated error to Hubble flow supernova using the distance modulus equation.

* "distanceladder" class is the default class to initialize the package from MontePython.

TABLE READERS
-------------

All table files take the same inputs. (1) file: filename to be read (2) index_array: columns to be read in file (3) header_length: how many lines to skip before data starts (4) delim: how the data is separated in file.

* "anchor_reader()" reads the anchor_data.txt file and seeks to import host name, distance in parsecs, and error in parsecs.

* "SHOES_ceph_table_reader()" reads the SHOES_cepheid_data.txt file and seeks to import host name, Cepheid ID, period in days, V band mag, error in V band mag, I band mag, error in I band mag, NIR band mag, error in NIR band mag, and metalicity.

* "LMC_ceph_reader()" reads in the LMC_ceph_data.txt file and seeks to import host name, Cepheid ID, period in days, Weisenheit H band mag, and error in Weisenheit H band mag.

* "TRGB_sn_reader()" reads in the b_TRGB_sn_data.txt file and seeks to import host name, Supernova ID, distance modulus, error in distance modulus, B band magnitude, and error in B band magnitude.

* "SHOES_sn_table_reader()" reads in the SHOES_sn_data.txt file and seeks to import host name, supernova ID, B band magnitude, and error in B band magnitude.

* "panth_sn_table_reader()" reads in the panth_sn_data.txt file and seeks to import supernova ID, heliocentric redshift, CMB frame redshift, B band magnitude, and error in B band magnitude.

COMPUTATION
-----------

The tabular data is read in according to the user inputs from "distanceladder.data".

The dictionaries are then created using the tabulated data.

Then the actual calculations begin. The following section will be divided into calibration choices and not be in the order of appearance.

CEPHEID CALIBRATION
-------------------

First, the absolute magnitude of Cepheids is calculated by combining distance information for the anchors and the Cepheids which exist in them. If multiple anchors are selected, the absolute cepheid magnitude is calculated for both anchors then combined using a weighted average where the weight is the error associated with each anchor's dMceph calculation. Then, the slope of the P-L relationship for Cepheids is calculated and combined using a weighted average.

With the slope, Mceph, and associated errors, the distance to all Cepheids is then computed using proto_Compute_mu. Then, the absolute supernova magnitude (Msn) is calculated by combining the distance measurements above with photometric measurement of the apparent B band magnitude of supernova.

Once Msn is calculated, the distance modulus to Hubble flow supernova in the pantheon data set is determined.

TRGB or CONCORDANT CALIBRATION
------------------------------

In this calibration scheme, we do not reproduce the TRGB tip absolute magnitude calculation but instead use the results of Freedman 2019 table 3. First, the absolute supernova magnitude is calculated by a linear regression of the distances to the TRGB host galaxies and apparent B band magnitude measurements of type Ia supernova in those host galaxies.

Once Msn is calculated, the distance modulus to Hubble flow supernova in the pantheon data set are determined.

LIKELIHOOD CALCULATION
----------------------

Once the distance modulus to Hubble flow supernova is calculated above using the distance ladder, the distance modulus to those same supernovae is calculated from the luminosity-redshift relationship derived from the cosmological model presented.

The difference between the local distance ladder and theory in regards to distance is then used as the chi-squared statistic.







