# MSFileParser
Perl library to parse Mass Spectrometry text based file formats

# Supported Formats
|File Format |File extension|Standard Type |File Type |Format Type |Module|
|:----|:----|:----|:----|:----|:----|
|MGF (Mascot genetic format)|.mgf|De facto|Text| Input MS/MS spectra|MSFileParser.pm|
|mzML|.mzML|HUPO-PSI|XML| Input MS/MS spectra (HUPO standard)|mzMLParser.pm|
|Tandem XML|.t.xml|De facto X!Tandem output|XML| Search result|PSMFileIO.pm|
|pepXML|.pep.xml|De facto TPP output|XML| Search result|PepXMLParser.pm|
|mzIdentML|.mzid|HUPO-PSI|XML| Search result (HUPO standard) |mzIdentMLParser.pm|
|Pkl (Waters MassLynx peak list format)|.pkl|De facto|Text| Input MS/MS spectra|MSFileParser|
|DTA|.dta|De facto|Text| Input MS/MS spectra|MSFileParser|
|mzQuantML|.mzq|HUPO-PSI|XML| Quantitation result|mzMLQuant.pm|

# Example Use cases
Please note that the library will provide data strcutures for further processing. The choice of algorithms (third-party or user-developed) is out of the current scope, even thou7gh some examples are provided.
## Before search (uninterpreted spectra)
* Deisotoping
* Peak picking from raw MS/MS spectra
* Noise removal
* Spectra filtering and selection based on disgnostic ions
* de novo sequencing
* PTM-specific spectra filtering

## Spectra search (MS/MS based matching)
* MS/MS peptide matching
* PTM site localization
* PTM identification
* iTRAQ quantitation
* TMT quantitation

## Post-search analysis
* FDR estimation 
* integrating search results from multiple algorithms
* visualization
