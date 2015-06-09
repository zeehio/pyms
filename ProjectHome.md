# PyMS #

PyMS is a modular software for processing of chromatography-mass spectrometry data, currently primarily tested on gas chromatography-mass spectrometry (GC-MS) data. PyMS is written in Python, an object oriented, general purpose programming language widely used in scientific computing.

**To download PyMS use the tab "Downloads" at the top of this page**

PyMS provides functions for reading the raw data files in ANDI-MS and JCAMP-DX formats,
noise smoothing, baseline correction, peak detection, peak integration, parsing of peak list files, and manipulation of these objects in Python. PyMS includes an implementation of  dynamic programming approach for peak alignment described earlier ([Robinson et al., BMC Bioinformatics 2007, 8:419](http://www.biomedcentral.com/1471-2105/8/419/)). PyMS implements MPI (Message Passing Interface) based parallel processing, thus enabling processing of metabolomic data on multiple CPUs or computing clusters.

PyMS has been extensively tested on gas chromatography--mass spectrometry (GC-MS) data. We are currently working on extending PyMS processing capabilities to liquid chromatography--mass spectrometry (LC-MS) data.

The PyMS project consists of three sub-projects:

  * [pyms](http://code.google.com/p/pyms/) -- the PyMS project (**this project**)
  * [pyms-docs](http://code.google.com/p/pyms-docs/) -- Documentation, including the PDF of the user guide with examples and tutorials
  * [pyms-test](http://code.google.com/p/pyms-test/) -- Tests & examples discussed in the user guide

The three parts of PyMS are hosted on Google Code as separate projects, and are freely available for download.

The data files for examples can be downloaded from
[here](http://bioinformatics.bio21.unimelb.edu.au/pyms/data/).

The PyMS APIs documentation is available [here](http://bioinformatics.bio21.unimelb.edu.au/pyms/api/index.html) (for programmers).

For pre-packaged PyMS dependencies (libraries required by PyMS) for Linux and Windows visit [University of Melbourne PyMS page](http://bioinformatics.bio21.unimelb.edu.au/pyms.html).