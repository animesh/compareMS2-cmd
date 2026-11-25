# compareMS2-cmd

compareMS2-cmd is the command line software of compareMS2. compareMS2 calculates the global similarity between tandem mass spectrometry datasets.
A GUI for this software is available at [compareMS2](https://github.com/524D/compareMS2).

[1. Introduction](#1-Introduction)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.1 What is compareMS2?](#11-What-is-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.2 How does compareMS2 differ from other tools?](#12-How-does-compareMS2-differ-from-other-tools)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.3 What can compareMS2 be used for?](#13-What-can-compareMS2-be-used-for)  
[2. Installing compareMS2](#2-Installing-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.1 Building from source](#21-Building-from-source)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.2 Running tests](#22-Running-tests)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.3 Release](#23-Release)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.4 Recent Changes](#24-Recent-Changes)  
[3. Using compareMS2](#3-Using-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1 Configuring compareMS2](#31-Configuring-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2 Calculating distance matrices](#32-Calculating-distance-matrices)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.3 Running compareMS2](#33-Running-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.4 Molecular phylogenetics](#34-Molecular-phylogenetics)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.5 Data quality control](#35-Data-quality-control)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.6 Experimental features](#36-Experimental-features)  
[4. Acknowledgements ](#4-Acknowledgements)  
[5. Further reading](#5-Further-reading)  

## 1. Introduction

### 1.1 What is compareMS2?

compareMS2 is a tool for direct comparison of tandem mass spectrometry datasets, typically from liquid chromatography-tandem mass spectrometry (LC-MS/MS), defining similarity as a function of shared (similar) spectra and distance as the inverse of this similarity. Data with identical spectral content thus have similarity 1 and distance 0. The similarity of datasets with no similar spectra tend to 0 (distance +∞) as the size of the sets go to infinity. The extremes of none or all spectra being similar between two LC-MS/MS datasets are unlikely to occur in reality.

### 1.2 How does compareMS2 differ from other tools?

Though compareMS2 is not limited to tandem mass spectra of peptides, it has seen most application to this type of data. There are four broad categories of tools for the analysis of peptide tandem mass spectra in mass spectrometry-based proteomics based on what prior information they utilize. compareMS2 belongs to a class of tools that do not use existing sequence data or libraries of spectra assigned to a specific peptide sequence, but compare tandem mass spectra directly with other tandem mass spectra:

<table>
    <tbody>
        <tr>
            <th rowspan=2 colspan=2></th>
            <th colspan=2>(translated) genome sequences available</th>
        </tr>
        <tr>
            <td>+</td>
            <td>-</td>
        </tr>
        <tr>
            <th rowspan=2>prior/other tandem <br>mass spectra available</th>
            <td>+</td>
            <td>spectral libraries (<a href="https://doi.org/10.1021/ac060279n">BiblioSpec</a>,
                <a href="https://doi.org/10.1002/pmic.200600625">SpectraST</a>, ...)</td>
            <td>direct comparison (<a href="https://doi.org/10.1002/rcm.6162"><b>compareMS2</b></a>,
                <a href="https://doi.org/10.1186/s12859-017-1514-2">DISMS2</a>, ...)</td>
        </tr>
        <tr>
            <td>-</td>
            <td>database search (<a href="https://doi.org/10.1002/(SICI)1522-2683(19991201)20:18<3551::AID-ELPS3551>3.0.CO;2-2">Mascot</a>,
                <a href="https://doi.org/10.1038/msb4100024">Comet</a>, ...)</td>
            <td><i>de novo</i> sequencing (<a href="https://doi.org/10.1002/pro.5560010902">LUTEFISK</a>,
                <a href="https://doi.org/10.1021/ac048788h">PepNovo</a>, ...)</td>
        </tr>
    </tbody>
</table>

### 1.3 What can compareMS2 be used for?

compareMS2 (and similar tools) have extremely broad utility, but have so far seen most utility in data quality control, food/feed species identification and molecular phylogenetics. Molecular phylogenetics is the study of evolution and relatedness of organisms, genes or proteins. The field dates back to [1960](https://doi.org/10.1073/pnas.46.10.1349) using patterns of tryptic peptides separated by paper chromatography. compareMS2 is a [21st-century analogue](https://doi.org/10.1021/acs.jproteome.1c00528), comparing patterns of tryptic peptides as analyzed by tandem mass spectrometry, with the difference that it can use thousands of peptides and that the tandem mass spectra are highly peptide-specific.

However, not only the amino acid sequences of the peptides affect the distance metric in compareMS2, but also the abundance (or coverage) of the proteins. compareMS2 can also be used to quantify the similarity of proteomes from different cell lines or tissues from the same species, before and independently of any protein identification by database or spectral library search.

## 2. Installing the compareMS2 command line software

To build compareMS2 and compareMS2_to_distance_matrices from source:

- [Install Git](https://git-scm.com/install/)
- [Install Go](https://go.dev/doc/install)
- In a CMD/PowerShell/Bash window, run: `go install github.com/524D/compareMS2-cmd/compareMS2@latest` and `go install github.com/524D/compareMS2-cmd/compareMS2_to_distance_matrices@latest`.
After building, the programs will be located in a location set by the `GOBIN` environment variable (typically `$HOME/go/bin`)

### 2.1 Building from source

`compareMS2` and `compareMS2_to_distance_matrices` are organized as separate nested Go modules (requiring Go 1.26.1+ and `CGO_ENABLED=0`). Build each binary from the repository root (`compareMS2-cmd`):

**Linux / macOS (Bash):**

```bash
cd compareMS2 && CGO_ENABLED=0 go build -o ../compareMS2 .
cd ../compareMS2_to_distance_matrices && CGO_ENABLED=0 go build -o ../compareMS2_to_distance_matrices .
cd ..
```

**Windows (PowerShell):**

```powershell
$env:CGO_ENABLED = "0"
cd compareMS2
go build -o ..\compareMS2.exe .
cd ..\compareMS2_to_distance_matrices
go build -o ..\compareMS2_to_distance_matrices.exe .
cd ..
```

### 2.2 Running tests

Unit tests are located in `compareMS2/compareMS2_test.go` and `compareMS2_to_distance_matrices/compareMS2_to_distance_matrices_test.go`. Run `go vet` and `go test` inside each module directory:

**Linux / macOS (Bash):**

```bash
cd compareMS2 && go vet ./... && go test -v ./...
cd ../compareMS2_to_distance_matrices && go vet ./... && go test -v ./...
cd ..
```

**Windows (PowerShell):**

```powershell
cd compareMS2
go vet ./...
go test -v ./...
cd ..\compareMS2_to_distance_matrices
go vet ./...
go test -v ./...
cd ..
```

### 2.3 Release

New releases are build with [GoReleaser](https://goreleaser.com/) using GitHub actions.
A pre-release can be build locally with:

```bash
goreleaser release --snapshot --clean
```

### 2.4 Recent Changes

- **Two-pass `-N` tie handling (`compareMS2`):** Fixed `PreCheckMGF` and `StoreTopNSpectrum` so tied boundary intensities (including datasets where all spectra have identical summed intensities) fill up to `TopN` spectra without displacing strictly higher-intensity spectra, and properly handle `-N 0`, `-N < 0`, and `-N >= spectrumCount`.
- **Two-pass blank MGF line handling (`compareMS2`):** Updated `PreCheckMGF` to skip empty lines during `-two-pass` scanning, preventing an index-out-of-range panic on blank lines.
- **Empty-dataset metric 2 guard (`compareMS2`):** Updated `WriteJSON` and `WriteTextOutput` to emit `"INF"` instead of `NaN` when either dataset has `0` spectra under `-d 2`.
- **Robust sample-to-species mapping (`compareMS2_to_distance_matrices`):** Hardened `readSampleToSpecies` and `sampleNameToSpeciesIndex` (`compareMS2_to_distance_matrices/compareMS2_to_distance_matrices.go`) by reading files with leading UTF-8 BOM (`\uFEFF`) removal and CRLF/CR/LF normalization, collapsing blank/consecutive tab fields per line (`len(parts) >= 2` non-empty trimmed tokens while guarding empty species names), extracting filenames across both POSIX (`/`) and Windows (`\`) paths via `sampleBaseName`, and storing a secondary `sampleLowerToIdx` map so lookups check exact basename first with case-insensitive fallback.
- **Unit test suites:** Added `compareMS2/compareMS2_test.go` and `compareMS2_to_distance_matrices/compareMS2_to_distance_matrices_test.go`.

## 3. Using compareMS2

This section describes how to run the compareMS2 command line software. Every compareMS2 analysis consists of two phases: (1) pairwise comparison of all LC-MS/MS datasets and (2) calculating a distance matrix from all pairwise comparisons. The pairwise comparison of the datasets provides a measure of their similarity. The distance matrix can then be used for molecular phylogenetics, data quality control or other purposes.

### 3.1 Configuring compareMS2  

The compareMS2 CLI has a small number of parameters, which are:

-A *First dataset filename (required)*  
-B *Second dataset filename (required)*  
-I *Minimum peaks (default 20)*  
-J *JSON output filename*  
-L *Minimum m/z (default 109)*  
-N *Top N spectra (default -1=all)*  
-R *RT range (start,end)*  
-U *Maximum m/z (default 2000)*  
-W *Scan range (start,end)*  
-X *Experimental output filename*  
-b *Bin size (default 0.2)*  
-c *Score cutoff (default 0.8)*  
-d *Version of distance metric (0, 1 or 2=default)*  
-f *Spectrum metric*  
-m *Minimum intensity (basepeak,total)*  
-n *Noise threshold*  
-o *Output filename*  
-p *Maximum difference in precursor mass (default 2.05)*  
-profile *CPU profiling output filename*  
-q *QC measure*  
-r *Maximum retention time difference (default 60)*  
-s *Scaling power (default 0.5)*  
-two-pass *Use two-pass processing (first pass to determine top-N cutoff, second pass to read spectra)*  
-w *Maximum scan number difference (default 10000)*  
-x *Experimental features*  

The compareMS2 GUI exposes some of these, and determine others automatically, e.g. the dataset filenames from a specified directory.

### 3.2 Calculating distance matrices

Distance matrices are calculated using a separate executable, `compareMS2_to_distance_matrices`. The `-i` input list accepts both `.json` and `.txt` pairwise comparison files produced by `compareMS2`. This tool can also average the distances for multiple replicates per species using a tab-delimited sample-to-species mapping file (`-x`), tolerating UTF-8 BOM headers, CRLF/CR/LF line endings, trailing or consecutive tabs, extra metadata columns beyond the second column, path prefixes (`sampleBaseName`), and Windows case differences (`sampleLowerToIdx` fallback after exact match) while ignoring rows with an empty species field. If no mapping file or valid row is provided, species names are inferred from the dataset filenames. Supported output formats are NEXUS (`-n`, default), MEGA < 12 (`-m`, written with CRLF line endings for Windows compatibility), MEGA 12+ (`-m2`), and JSON (`-J`). [MEGA](https://www.megasoftware.net/) is recommended for creating trees from compareMS2 results.

### 3.3 Running compareMS2

Run pairwise comparisons between MGF files (optionally enabling `-two-pass -N` for top-N spectrum filtering and `-J` for JSON output):

```bash
compareMS2 -A sampleA.mgf -B sampleB.mgf -o AB.txt -J AB.json
compareMS2 -A sampleA.mgf -B sampleC.mgf -two-pass -N 2000 -o AC.txt -J AC.json
```

Then combine pairwise results into distance matrices:

```bash
compareMS2_to_distance_matrices -i results_list.txt -o distance_matrix -x species_map.txt -n -m -m2 -J
```

### 3.4 Molecular phylogenetics

We recommend [MEGA](https://www.megasoftware.net/) creating phylogenetic trees from compareMS2 results. However, most phylogenetic software can take distance matrices as input for UPGMA analysis. This was the original use for which compareMS2 was developed, see the [2012 paper](https://doi.org/10.1002/rcm.6162).

### 3.6 Experimental features

Starting in version 2.0, we have begun to include experimental features in compareMS2. These are only available on the command line, but allow extraction of additional information from the comparisons, such as the distribution of similarity between tandem mass spectra as function of precursor mass measurement error, allowing identification of isotope errors and charge state distributions *before* any database search:

![Experimental feature](./pictures/experimental_features.png)  
Figure 2. Similarity (spectral angle from 0 to 1) of tandem mass spectra plotted against precursor *m*/*z* difference, revealing isotope errors up to at least 2 (corresponding to bands at *m*/*z* difference 2/3 and 2/5) and charge states up to 6 (corresponding to the band at *m*/*z* difference 1/6).

## 4. Acknowledgements

The developers wish to thank Dr. Michael Dondrup at the University of Bergen for providing changes and additions to make compareMS2 work under macOS. All users and beta testers are also acknowledged for their valuable feedback that helped to improve compareMS2.

## 5. Further reading

compareMS2 and related applications have been described or used in a number of papers:

compareMS2 2.0: An Improved Software for Comparing Tandem Mass Spectrometry Datasets, Marissen M, Varunjikar MS, Laros JFJ, Rasinger JD, Neely BA and Palmblad M, *J. Proteome Res.*, *in press*, 2022, [doi.org/10.1021/acs.jproteome.2c00457](https://doi.org/10.1021/acs.jproteome.2c00457)

Shotgun proteomics approaches for authentication, biological analyses, and allergen detection in feed and food-grade insect species, Varunjikar MS, Belghit I, Gjerde J, Palmblad M, Oveland E and Rasinger JD, *Food Control* **131**, 2022, [doi.org/10.1016/j.foodcont.2022.108888](https://doi.org/10.1016/j.foodcont.2022.108888)

Comparing novel shotgun DNA sequencing and state-of-the-art proteomics approaches for authentication of fish species in mixed samples, Varunjikar MS, Moreno-Ibarguen C, Andrade-Martinez JS, Tung HS, Belghit I, Palmblad M, Olsvik PA, Reyes A, Rasinger JD and Lie KK, *Food Control* **131**:108417, 2022, [doi.org/10.1016/j.foodcont.2021.108417](https://doi.org/10.1016/j.foodcont.2021.108417)

Rewinding the molecular clock: looking at pioneering molecular phylogenetics experiments in the light of proteomics, Neely B and Palmblad M, *J. Proteome Res.* **20(10)**:4640-4645, 2021, [doi.org/10.1021/acs.jproteome.1c00528](https://doi.org/10.1021/acs.jproteome.1c00528)

Future feed control – Tracing banned bovine material in insect meal. Belghit I, Varunjikar M, Lecrenier MC, Steinhilber A, Niedzwiecka A, Wang YV, Dieu M, Azzollini D, Lie K, Lock EJ, Berntssen MHG, Renard P, Zagon J, Fumière O, van Loon JJA, Larsen T, Poetz O, Braeuning A, Palmblad M and Rasinger JD, *Food Control* **128**:108183, 2021, [doi.org/10.1016/j.foodcont.2021.108183](https://doi.org/10.1016/j.foodcont.2021.108183)

Species-Specific Discrimination of Insect Meals for Aquafeeds by Direct Comparison of Tandem Mass Spectra. Belghit I, Lock EJ, Fumière O, Lecrenier MC, Renard P, Dieu M, Berntssen MHG, Palmblad M and Rasinger JD, *Animals* **9(5)**:222, 2019 [doi.org/10.3390/ani9050222](https://doi.org/10.3390/ani9050222)

Palaeoproteomics of bird bones for taxonomic classification. Horn IR, Kenens Y, Palmblad M, van der Plas-Duivesteijn SJ, Langeveld BW, Meijer HJM, Dalebout H, Marissen RJ, Fischer A, Vincent Florens FB, Niemann J, Rijsdijk KF, Schulp AS, Laros JFJ and Gravendeel B, *Zoological Journal of the Linnean Society* **186(3)**:650–665, 2019, [doi.org/10.1093/zoolinnean/zlz012](https://doi.org/10.1093/zoolinnean/zlz012)

Species and tissues specific differentiation of processed animal proteins in aquafeeds using proteomics tools. Rasinger JD, Marbaix H, Dieu M, Fumière O, Mauro S, Palmblad M, Raes M and Berntssen MHG, *J. Proteomics* **147**:125-131, 2016, [doi.org/10.1016/j.jprot.2016.05.036](https://doi.org/10.1016/j.jprot.2016.05.036)   

Authentication of closely related fish and derived fish products using tandem mass spectrometry and spectral library matching. Nessen M, van der Zwaan D, Greevers S, Dalebout H, Staats M, Kok E and Palmblad M, *J. Agric. Food Chem.* **64(18)**:3669-3677, 2016, [doi.org/10.1021/acs.jafc.5b05322](https://doi.org/10.1021/acs.jafc.5b05322)

Identification of meat products by shotgun spectral matching. Ohana D, Dalebout H, Marissen RJ, Wulff J, Bergquist J, Deelder AM and Palmblad M, *Food Chem.* **203**:28-34, 2016, [doi.org/10.1016/j.foodchem.2016.01.138](https://doi.org/10.1016/j.foodchem.2016.01.138)   

Differentiating samples and experimental protocols by direct comparison of tandem mass spectra. van der Plas-Duivesteijn SJ, Wulff T, Klychnikov O, Ohana D, Dalebout H, van Veelen PA, de Keijzer J, Nessen MA, van der Burgt YEM, Deelder AM and Palmblad M, *Rapid Commun. Mass Spectrom.* **30**:731-738, 2016, [doi.org/10.1002/rcm.7494](https://doi.org/10.1002/rcm.7494)   

Authentication of Fish Products by Large-Scale Comparison of Tandem Mass Spectra. Wulff  T, Nielsen ME, Deelder AM, Jessen F and Palmblad M, *J. Proteome Res.* **12(11)**:5253-5259, 2013, [doi.org/10.1021/pr4006525](https://doi.org/10.1021/pr4006525)

Molecular phylogenetics by direct comparison of tandem mass spectra. Palmblad M and Deelder AM, *Rapid Commun. Mass Spectrom.* **26(7)**:728-732, 2012, [doi.org/10.1002/rcm.6162](https://doi.org/10.1002/rcm.6162)
