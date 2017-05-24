# MS PEAKLIST EXPORT



## Programming language
R: The Comprehensive R Archive Network

https://www.r-project.org/



## Scope of the software
The software imports mass spectra (stored in different file formats), performs spectral preprocessing onto the imported spectra (with different parameters), computes the peak picking and peak alignment onto the preprocessed spectra and generates a peaklist matrix.

The peaklist matrix is built in such a way that each row corresponds to a spectra and each column to a mass peak in the spectral dataset. Further columns are added: the "Sample" column and the "Class" column. The sample entry is determined by the name of the spectrum file, while the class corresponds to the name of the folder the spectral files are stored into.
Further data can be manually added to the matrix, such as demographical data (age, sex, etc...).

The peaklist matrix is a classical statistical matrix, which can therefore be imported into any program that can perform statistical analysis, such as Orange Canvas (https://orange.biolab.si/) and Weka (http://www.cs.waikato.ac.nz/ml/weka/).



## Type of data and its organization

#### Type of data
The software can import spectral files in different file formats: imzML (imaging datasets), Xmass (Bruker Daltonics Xmass format), TXT and CSV files. 


#### Organization of the data
The spectral files should be organized in such a way that each file corresponds to a patient's spectrum (or imaging dataset in the case of imzML) and all the files of patients belonging to the same class (e.g. pathology) are put in the same folder, named according to the class itself.
In this way, the software will read at first the folders in the provided sample folder and then the spectral files in each folder, therefore storing the information regarding the sample name and the belonging class.
