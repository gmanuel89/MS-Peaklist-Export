# MS PEAKLIST EXPORT

***

## Program version
In order for this WIKI to be applicable, the version of the program must be equal to or higher than **2017.06.19.0**.

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

Additionally, a file with the list of the parameters is generated and saved next to the peaklist file.

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
    * **_Data transformation_**: selects if data transformation should be performed (applies a mathematical operation to all the intensities, among "Square root", "Natural logarithm", "Decimal Logarithm" and "Binary Logarithm").
    * **_Smoothing_**: defines the algorithm for the spectral smoothing ("Savitzky-Golay", "Moving Average", "None") and the strength of the smoothing ("medium", "strong", "stronger").
        * *Savitzky-Golay*: It is a process known as convolution, by fitting successive sub-sets of adjacent data points with a low-degree polynomial by the method of linear least squares. (A. Savitzky and M. J. Golay. 1964. Smoothing and differentiation of data by simplified least squares procedures. Analytical chemistry, 36(8), 1627-1639).
        * *Moving Average*: Given a series of numbers and a fixed subset size, the first element of the moving average is obtained by taking the average of the initial fixed subset of the number series. Then the subset is modified by "shifting forward"; that is, excluding the first number of the series and including the next value in the subset. (Booth et al., San Francisco Estuary and Watershed Science, Volume 4, Issue 2, 2006).
    * **_Baseline subtraction_**: defines the algorithm for baseline subtraction ("SNIP", "TopHat", "ConvexHull", "median", "None"). Before selecting the baseline subtraction algorithm, a number defining the value of the specific parameter for the algorithm can be inserted, and the program will read it while setting the baseline subtraction algorithm.
        * *SNIP*: Statistics-sensitive Non-linear Iterative Peak-clipping algorithm (C.G. Ryan, E. Clayton, W.L. Griffin, S.H. Sie, and D.R. Cousens. 1988. Snip, a statistics-sensitive background treatment for the quantitative analysis of pixe spectra in geoscience applications. Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms, 34(3): 396-402).
        * *TopHat*: This algorithm applies a moving minimum (erosion filter) and subsequently a moving maximum (dilation filter) filter on the intensity values (M. van Herk. 1992. A Fast Algorithm for Local Minimum and Maximum Filters on Rectangular and Octagonal Kernels. Pattern Recognition Letters 13.7: 517-521).
        * *ConvexHull*: The baseline estimation is based on a convex hull constructed below the spectrum (Andrew, A. M. 1979. Another efficient algorithm for convex hulls in two dimensions. Information Processing Letters, 9(5), 216-219).
        * *Median*: This baseline estimation uses a moving median.
    * **_Normalization_**: defines the algorithm for normalization ("TIC", "RMS", "PQN", "median", "None"). Before selecting the normalization algorithm, a number defining the normalization mass range can be inserted, and the program will read it while setting the normalization algorithm.
        * *TIC (Total Ion Current)*: It divides the intensities of the spectrum by the sum of all the intensity values of the spectrum itself (the sum of all the intensities being the spectrum's total current). It becomes less suitable when very intense peak(s) (compared to the others) are present in the spectrum.
        * *RMS (Root Mean Square)*: It divides the intensities of the spectrum by the square root of the sum of all the intensity values of the spectrum itself squared. Like the TIC, it becomes less suitable when very intense peak(s) (compared to the others) are present in the spectrum.
        * *PQN (Probabilistic Quotient Normalization)*: It calibrates the spectra using the TIC normalization; then, a median reference spectrum is obtained; the quotients of all intensities of the spectra with those of the reference spectrum are calculated; the median of these quotients is calculated for each spectrum; finally, all the intensity values of each spectrum are divided by the median of the quotients for the spectrum (F. Dieterle, A. Ross, G. Schlotterbeck, and Hans Senn. 2006. Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures. Application in 1H NMR metabonomics. Analytical Chemistry 78(13): 4281-4290).
        * *Median*: It divides the intensities of the spectrum by the median of all the intensity values of the spectrum itself. It has been proved to be the most robust normalization method.
    * **_Align spectra_**: selects if alignment of spectra should be performed, by generating a calibration curve ("cubic", "quadratic", "linear", "lowess") employing an automatically generated peaklist ("auto") as reference or by taking the peaks of the "average spectrum" or of the "skyline spectrum" as reference.
    * **_Preprocess spectra in packages of_**: defines the number of spectra to be taken at a time for preprocessing, when the computer resources are limited (taking all the spectra in RAM could cause the computer to freeze).
    * **_Tolerance (in ppm)_**: defines the tolerance (in ppm, parts per million) for the spectral alignment, peak alignment and database-sample signal match. For linear TOF mode the tolerance should be set to 1000 ppm (0.1%, 4 Da at 4000 Da), while for reflectron TOF mode the tolerance should be set to 100 ppm (0.01%, 0.2 Da at 2000 Da).
    * **_Commit preprocessing_**: stores the preprocessing parameters to be applied for analysis.

* **Browse spectra**: selects the single spectral "file" (e.g. single imzML file) or the "folder" in which all the spectral files are stored.

* **Browse output folder**: selects the folder in which all the output files should be saved.

* **File type export**: sets the type of the exported files ("CSV", "XLS" or "XLSX").

* **Import and preprocess spectra...**: imports the spectra and computes preprocessing (according to the specified parameters), averaging and alignment.

* **Peak picking...**: performs the peak picking on the imported spectra.

* **Export peaklist...**: saves the peaklist matrix file.

* **Mean +/- sd of number of signals...**: displays the mean number of signals of the spectra in the dataset, along with the standard deviation and the coefficient of variation of the signal number.

* **Dump spectra files...**: generates a folder (named as the peaklist matrix file) in which all the spectral files are saved (MSD or TXT).

* **Quit**: close the program and the R session.

***

## Example

#### Organization of the data
Example of folder hierarchy:

/Input folder/Class1/Spectra (imzML, Xmass, TXT, CSV, MSD) files

/Input folder/Class2/Spectra (imzML, Xmass, TXT, CSV, MSD) files

/Input folder/Class3/Spectra (imzML, Xmass, TXT, CSV, MSD) files
