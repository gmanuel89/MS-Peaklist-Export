# MS PEAKLIST EXPORT

***

## Programming language
[R: The Comprehensive R Archive Network](https://www.r-project.org/)

***

## Scope of the software
The software imports mass spectra (stored in different file formats), performs spectral preprocessing onto the imported spectra (with different parameters), computes the peak picking and peak alignment onto the preprocessed spectra and generates a peaklist matrix.

The peaklist matrix is built in such a way that each row corresponds to a spectrum and each column to a mass peak in the spectral dataset. Further columns are added: the "Sample" column and the "Class" column. The sample entry is determined by the name of the spectrum file, while the class corresponds to the name of the folder the spectral files are stored into.

Further data can be manually added to the matrix, such as demographical data (age, sex, etc...).

The peaklist matrix is a classical statistical matrix, which can therefore be imported into any program that can perform statistical analysis, such as [Orange Canvas](https://orange.biolab.si/) and [Weka](http://www.cs.waikato.ac.nz/ml/weka/).

***

## Type of data and its organization

#### Type of data
The software can import spectral files in different file formats: [imzML](https://ms-imaging.org/wp/introduction/) (imaging datasets), [Xmass](https://www.bruker.com/products/mass-spectrometry-and-separations/ms-software.html) (Bruker Daltonics Xmass format), TXT,  CSV and MSD files. 


#### Organization of the data
The spectral files should be organized in such a way that each file corresponds to a patient's spectrum (or imaging dataset in the case of imzML format) and all the files of patients belonging to the same class (e.g. pathology) are put in the same folder, named according to the class itself.

In this way, the software will read at first the folders in the provided sample folder and then the spectral files in each folder, therefore storing the information regarding the sample name and the belonging class.


#### Output data
The software generates a CSV file corresponding to the peaklist intensity matrix, in which each row is a spectrum (pixel or patient) and each column is a peak,with additional "Sample" and "Class" columns.

The spectra files are placed in a folder with the same name as the peaklist file: if "MSD" is selected as the file format, one MSD file for each spectrum is generated, with the peaks embedded in the same file; if "TXT" is selected as the file format, two TXT files are generated, one forthe spectrum and one for the peaks.

***

## Buttons, entries and operations

* **Spectra format**: selects the format of the input spectral files (imzML, Xmass, TXT, CSV or MSD).

* **Peak picking algorithm**: selects the algorithm to be employed for peak picking ("Friedman's Super Smoother" or "Median Absolute Deviation").

* **Signal-to-noise ratio**: defines the signal-to-noise ratio to be used as a threshold for peak picking, after noise estimation.

* **Peak picking mode**: defines if "all" the peaks should be kept or only the "most intense" for each spectrum.

* **Most intense signals to take**: defines the number of most intense signals to preserve (if "most intense" is selected as "Peak picking mode").

* **Peak deisotoping enveloping**: defines if "Peak deisotoping" (preserve only the monoisotipic peak of the isotope cluster) or "Peak enveloping" (preserve only the most intense peak of the isotope cluster) should be performed after peak picking.

* **Low-intensity peak removal percentage threshold**: defines the percentage threshold according to which all the peaks that have an intensity below the selected percentage of the base peak are discarded.

* **Intensity threshold method**: selects the method for base peak identification, according to which the low-intensity peak removal is performed for each spectrum evaluating the intensity of the base peak of that spectrum ("element-wise") or the base peak of the whole spectral dataset ("whole").

* **Peak filtering threshold percentage**: defines the percentage threshold according to which all the peaks that are not present in at least the selected percentage of the spectra are discarded.

* **Peak filtering mode**: selects the method for peak filtering, according to which the percentage threshold is relative to the whole spectral dataset ("whole dataset") or relative to the classes of spectra ("class-wise").

* **Average replicates**: selects if an average spectrum is to be generated for each patient ("YES") or all the individual spectra corresponding to pixels are preserved for each patient ("NO").

* **Allow parallel computing**: enables the parallel computation (multi-CPU computation).

* **Spectra preprocessing parameters**: sets the parameters for spectral preprocessing.
    * **_Mass range_**: defines the mass range to which the imported spectra should be cut.
    * **_TOF mode_**: defines if the TOF has been used in the "Linear" or "Reflectron" mode (to adjust the parameters for spectral preprocessing and peak picking).
    * **_Transform the data_**: selects if data transformation should be performed (applies a mathematical operation to all the intensities, such as square root or logarithm with base e, 2 and 10).
    * **_Smoothing_**: defines the algorithm for the smoothing ("Savitzky-Golay", "Moving Average", "None") and the strength of the smoothing ("medium", "strong", "stronger").
    * **_Baseline subtraction_**: defines the algorithm for baseline subtraction ("SNIP", "TopHat", "ConvexHull", "median"). Before selecting the baseline subtraction algorithm, a number defining the value of the specific parameter for the algorithm can be inserted, and the program will read it while setting the baseline subtraction algorithm.
    * **_Normalization_**: defines the algorithm for normalization ("TIC", "RMS", "PQN", "median", "None"). Before selecting the normalization algorithm, a number defining the normalization mass range can be inserted, and the program will read it while setting the normalization algorithm.
    * **_Align spectra_**: selects if alignment of spectra should be performed, by employing an automatically generated peaklist ("auto") as reference or by taking the peaks of the "average spectrum" or of the "skyline spectrum" as reference.
    * **_Preprocess spectra in packages of_**: defines the number of spectra to be taken at a time for preprocessing, when the computer resources are limited (taking all the spectra in RAM could cause the computer to freeze).
    * **_Commit preprocessing_**: stores the preprocessing parameters to be applied for analysis.

* **Browse spectra**: selects the single spectral "file" (e.g. single imzML file) or the "folder" in which all the spectral files are stored.

* **Browse output folder**: selects the folder in which all the output files should be saved.

* **File type export**: sets the type of the exported files ("CSV", "XLS" or "XLSX").

* **Import and preprocess spectra...**: imports the spectra and computes preprocessing (according to the specified parameters), averaging and alignment.

* **Peak picking...**: performs the peak picking on the imported spectra.

* **Export peaklist...**: saves the peaklist matrix file.

* **Mean +/- sd of number of signals...**: displays the mean number of signals of the spectra in the dataset, along with the standard deviation and the coefficient of variation of the signal number.

* **Dump spectra files...**: generates a folder (named as the peaklist matrix file) in which all the spectral files are saved (MSD or TXT).

* **Quit**: close the program and the R.

***

## Example

#### Organization of the data
Example of folder hierarchy:

/Input folder/Class1/Spectra (imzML, Xmass, TXT, CSV, MSD) files

/Input folder/Class2/Spectra (imzML, Xmass, TXT, CSV, MSD) files

/Input folder/Class3/Spectra (imzML, Xmass, TXT, CSV, MSD) files
