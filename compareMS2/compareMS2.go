package main

// compareMS2 - comparison of tandem mass spectra datasets, e.g. from LC-MS/MS
//
// MIT License
//
// Copyright (c) Magnus Palmblad and Rob Marissen
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Contact information: n.m.palmblad@lumc.nl

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
)

const (
	dotprodHistogramBins   = 200
	massdiffHistogramBins  = 320
	massdiffHistogramRange = 3.2
)

type Peak struct {
	Mz        float64
	Intensity float64
}

// Spectrum represents a single MS2 spectrum
type Spectrum struct {
	Scan  int
	RT    float64
	Peaks []Peak
	// FIXME: Peaks Not used in ReadMGF2
	nPeaks            int // Number of peaks, used after removing the actual peaks to save memory
	Charge            int
	Bin               []float32
	PrecursorMz       float64
	BasepeakIntensity float64
	TotalIonCurrent   float64
	SumIntensity      float64
	SpecIdx           int // Index of the spectrum in the file (starts at 0)
}

// Dataset represents a collection of spectra from an MGF file
type Dataset struct {
	ID                     rune
	Filename               string
	ScanNumbersCouldBeRead bool
	RTsCouldBeRead         bool
	Spectra                []Spectrum
	Intensities            []float64 // FIXME: Not used in ReadMGF2
	IntensityCutoff        float64
	IntensityTieAllowance  int
	TwoPassFilterActive    bool
	SpectrumCount          int // Total number of spectra, including the onces that were filtered out
}

// Parameters holds all command-line and processing parameters
type Parameters struct {
	DatasetAFilename           string
	DatasetBFilename           string
	OutputFilename             string
	JSONFilename               string
	ExperimentalOutputFilename string
	MinBasepeakIntensity       float64
	MinTotalIonCurrent         float64
	MaxScanNumberDifference    float64
	MaxRTDifference            float64
	MaxPrecursorDifference     float64
	StartScan                  int
	EndScan                    int
	StartRT                    float64
	EndRT                      float64
	Cutoff                     float64
	Scaling                    float64
	Noise                      float64
	Metric                     int
	SpectrumMetric             int
	QC                         int
	BinSize                    float64
	MinPeaks                   int
	PeakCount                  int // Keep track of the highest peak count, only for backward compatibility
	MinMz                      float64
	MaxMz                      float64
	NBins                      int
	TopN                       int
	ExperimentalFeatures       int
	CPUProfile                 string
	TwoPass                    bool
}

const (
	NOT_IN_RANGE = iota
	IN_RANGE
	IN_RANGE_AND_FIRST_PEAK
)

// NewParameters creates a Parameters struct with default values
func NewParameters() *Parameters {
	return &Parameters{
		OutputFilename:             "output.txt",
		JSONFilename:               "",
		ExperimentalOutputFilename: "experimental_output.txt",
		MinBasepeakIntensity:       0,
		MinTotalIonCurrent:         0,
		MaxScanNumberDifference:    10000,
		MaxRTDifference:            60,
		MaxPrecursorDifference:     2.05,
		StartScan:                  1,
		EndScan:                    1000000,
		StartRT:                    0,
		EndRT:                      100000,
		Cutoff:                     0.8,
		Scaling:                    0.5,
		Noise:                      0,
		Metric:                     2,
		SpectrumMetric:             0,
		QC:                         0,
		BinSize:                    0.2,
		MinPeaks:                   20,
		PeakCount:                  0,
		MinMz:                      109,
		MaxMz:                      2000,
		NBins:                      9455,
		TopN:                       -1,
		ExperimentalFeatures:       0,
		CPUProfile:                 "",
		TwoPass:                    false,
	}
}

// ParseFlags parses command-line flags
func (p *Parameters) ParseFlags() error {
	flag.StringVar(&p.DatasetAFilename, "A", "", "First dataset filename (required)")
	flag.StringVar(&p.DatasetBFilename, "B", "", "Second dataset filename (required)")
	flag.StringVar(&p.OutputFilename, "o", p.OutputFilename, "Output filename")
	flag.StringVar(&p.JSONFilename, "J", p.JSONFilename, "JSON output filename")
	flag.StringVar(&p.ExperimentalOutputFilename, "X", p.ExperimentalOutputFilename, "Experimental output filename")
	flag.Float64Var(&p.Cutoff, "c", p.Cutoff, "Score cutoff")
	flag.Float64Var(&p.Scaling, "s", p.Scaling, "Scaling power")
	flag.Float64Var(&p.Noise, "n", p.Noise, "Noise threshold")
	flag.Float64Var(&p.MaxScanNumberDifference, "w", p.MaxScanNumberDifference, "Maximum scan number difference")
	flag.Float64Var(&p.MaxRTDifference, "r", p.MaxRTDifference, "Maximum retention time difference")
	flag.Float64Var(&p.MaxPrecursorDifference, "p", p.MaxPrecursorDifference, "Maximum difference in precursor mass")
	flag.IntVar(&p.Metric, "d", p.Metric, "Distance metric (0, 1 or 2)")
	flag.IntVar(&p.SpectrumMetric, "f", p.SpectrumMetric, "Spectrum metric")
	flag.IntVar(&p.QC, "q", p.QC, "QC measure")
	flag.Float64Var(&p.BinSize, "b", p.BinSize, "Bin size")
	flag.IntVar(&p.MinPeaks, "I", p.MinPeaks, "Minimum peaks")
	flag.Float64Var(&p.MinMz, "L", p.MinMz, "Minimum m/z")
	flag.Float64Var(&p.MaxMz, "U", p.MaxMz, "Maximum m/z")
	flag.IntVar(&p.TopN, "N", p.TopN, "Top N spectra")
	flag.IntVar(&p.ExperimentalFeatures, "x", p.ExperimentalFeatures, "Experimental features")
	flag.StringVar(&p.CPUProfile, "profile", "", "CPU profiling output file")
	flag.BoolVar(&p.TwoPass, "two-pass", false, "Use two-pass processing (first pass to determine top-N cutoff, second pass to read spectra)")

	// Custom flags for ranges
	var scanRange string
	var rtRange string
	var intensityRange string
	flag.StringVar(&scanRange, "W", "", "Scan range (start,end)")
	flag.StringVar(&rtRange, "R", "", "RT range (start,end)")
	flag.StringVar(&intensityRange, "m", "", "Minimum intensity (basepeak,total)")

	flag.Parse()

	// Validate required parameters
	if p.DatasetAFilename == "" || p.DatasetBFilename == "" {
		return fmt.Errorf("both -A and -B dataset filenames are required")
	}

	// Parse scan range
	if scanRange != "" {
		parts := strings.Split(scanRange, ",")
		if len(parts) == 2 {
			if v, err := strconv.Atoi(parts[0]); err == nil {
				p.StartScan = v
			}
			if v, err := strconv.Atoi(parts[1]); err == nil {
				p.EndScan = v
			}
		}
	}

	// Parse RT range
	if rtRange != "" {
		parts := strings.Split(rtRange, ",")
		if len(parts) == 2 {
			if v, err := strconv.ParseFloat(parts[0], 64); err == nil {
				p.StartRT = v
			}
			if v, err := strconv.ParseFloat(parts[1], 64); err == nil {
				p.EndRT = v
			}
		}
	}

	// Parse intensity range
	if intensityRange != "" {
		parts := strings.Split(intensityRange, ",")
		if len(parts) == 2 {
			if v, err := strconv.ParseFloat(parts[0], 64); err == nil {
				p.MinBasepeakIntensity = v
			}
			if v, err := strconv.ParseFloat(parts[1], 64); err == nil {
				p.MinTotalIonCurrent = v
			}
		}
	}

	return nil
}

// quickSelect finds the k-th largest element (for filtering)
func quickSelect(a []float64, left, right, k int) float64 {
	if left == right {
		return a[left]
	}
	pivotIndex := left + rand.Intn(right-left+1)
	pivotIndex = partition(a, left, right, pivotIndex)
	if k == pivotIndex {
		return a[k]
	} else if k < pivotIndex {
		return quickSelect(a, left, pivotIndex-1, k)
	} else {
		return quickSelect(a, pivotIndex+1, right, k)
	}
}

func partition(a []float64, left, right, pivotIndex int) int {
	pivot := a[pivotIndex]
	a[pivotIndex], a[right] = a[right], a[pivotIndex]
	storeIndex := left
	for i := left; i < right; i++ {
		if a[i] >= pivot {
			a[i], a[storeIndex] = a[storeIndex], a[i]
			storeIndex++
		}
	}
	a[storeIndex], a[right] = a[right], a[storeIndex]
	return storeIndex
}

// parseNumber extracts the first number from a string
func parseNumber(s string) string {
	// Find first digit
	start := -1
	for i, c := range s {
		if c >= '0' && c <= '9' {
			start = i
			break
		}
	}
	if start == -1 {
		return ""
	}

	// Find end of number (including decimal point)
	end := start
	for i := start; i < len(s); i++ {
		c := s[i]
		if (c >= '0' && c <= '9') || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E' {
			end = i + 1
		} else {
			break
		}
	}
	return s[start:end]
}

// PreCheckMGF computer the intensity cutoff for top-N filtering
func PreCheckMGF(par *Parameters, dataset *Dataset) error {
	if par.TopN < 0 {
		dataset.IntensityCutoff = 0.0 // When TopN filtering is disabled, the cutoff is zero
		dataset.IntensityTieAllowance = 0
		dataset.TwoPassFilterActive = false
		return nil
	}
	if par.TopN == 0 {
		dataset.IntensityCutoff = math.Inf(1)
		dataset.IntensityTieAllowance = 0
		dataset.TwoPassFilterActive = true
		fmt.Printf("done (ion intensity threshold %.3f)\n", dataset.IntensityCutoff)
		return nil
	}

	file, err := os.Open(dataset.Filename)
	if err != nil {
		return fmt.Errorf("error opening dataset %c %s for reading: %v", dataset.ID, dataset.Filename, err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	specStatus := IN_RANGE
	var intensities []float64
	nPeaks := 0
	SumIntensity := 0.0

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		// Check if line is peak data
		if line[0] >= '0' && line[0] <= '9' {
			nPeaks++
			if specStatus == NOT_IN_RANGE {
				continue
			}
			parts := strings.Fields(line)
			if len(parts) >= 2 {
				mz, err1 := strconv.ParseFloat(parts[0], 64)
				intensity, err2 := strconv.ParseFloat(parts[1], 64)
				if err1 == nil && err2 == nil {
					if mz >= par.MinMz && mz <= par.MaxMz {
						// Compute sum of intensities for top-N filtering
						SumIntensity += intensity
					}
					// FIXME: log errors in parsing mz or intensity, but continue processing
				}
			}
		} else if strings.HasPrefix(line, "BEGIN") {
			specStatus = IN_RANGE
		} else if strings.HasPrefix(line, "SCANS=") || strings.HasPrefix(line, "###MSMS:") {
			numStr := parseNumber(line)
			if scan, err := strconv.Atoi(numStr); err == nil {
				if scan < par.StartScan || scan > par.EndScan {
					specStatus = NOT_IN_RANGE
				}
			}
		} else if strings.HasPrefix(line, "TITLE=") && strings.Contains(line, "scan=") {
			// Extract scan from TITLE
			scanIdx := strings.Index(line, "scan=")
			if scanIdx != -1 {
				numStr := parseNumber(line[scanIdx+5:])
				if scan, err := strconv.Atoi(numStr); err == nil {
					if scan < par.StartScan || scan > par.EndScan {
						specStatus = NOT_IN_RANGE
					}
				}
			}
		} else if strings.HasPrefix(line, "RTINSECONDS=") {
			numStr := parseNumber(line)
			if rt, err := strconv.ParseFloat(numStr, 64); err == nil {
				if rt < par.StartRT || rt > par.EndRT {
					specStatus = NOT_IN_RANGE
				}
			}
		} else if strings.HasPrefix(line, "END") {
			if specStatus != NOT_IN_RANGE {
				intensities = append(intensities, SumIntensity)
			}
			SumIntensity = 0.0
		}
	}
	// When we don't have enough spectra to fill the top-N, disable filtering
	intensityCutoff := 0.0
	if len(intensities) > par.TopN {
		intensityCutoff = quickSelect(intensities, 0, len(intensities)-1, par.TopN-1)
		nGreater := 0
		for _, val := range intensities {
			if val > intensityCutoff {
				nGreater++
			}
		}
		dataset.IntensityCutoff = intensityCutoff
		dataset.IntensityTieAllowance = par.TopN - nGreater
		dataset.TwoPassFilterActive = true
	} else {
		dataset.IntensityCutoff = 0.0
		dataset.IntensityTieAllowance = 0
		dataset.TwoPassFilterActive = false
	}
	fmt.Printf("done (ion intensity threshold %.3f)\n", intensityCutoff)
	return scanner.Err()
}

// StoreTopNSpectrum processes and stores a spectrum with TopN filtering
func StoreTopNSpectrum(par *Parameters, dataset *Dataset, spec Spectrum, peaks []Peak, nPeaks int, minIdx *int) {
	// Keep track of the spectrum with the highest peakcount
	// This is only for backward compatibility, it is unused
	if nPeaks > par.PeakCount {
		par.PeakCount = nPeaks
	}
	if par.TwoPass {
		keep := !dataset.TwoPassFilterActive || spec.SumIntensity > dataset.IntensityCutoff
		if !keep && spec.SumIntensity == dataset.IntensityCutoff && dataset.IntensityTieAllowance > 0 {
			dataset.IntensityTieAllowance--
			keep = true
		}
		if keep {
			ScaleNormalizeBinSpectrum(&spec, par, peaks)
			dataset.Spectra = append(dataset.Spectra, spec)
		}
	} else {
		// if TopN filtering is enabled, keep track of the topN intensities and only store spectra that are in the topN
		if par.TopN == 0 {
			return
		}
		if par.TopN > -1 {
			if len(dataset.Spectra) < par.TopN { // if we haven't reached topN spectra yet, just add the spectrum
				// Pre-process spectrum.
				ScaleNormalizeBinSpectrum(&spec, par, peaks)
				dataset.Spectra = append(dataset.Spectra, spec)
				// Keep track of the spectrum with the lowest intensity in the current topN
				if dataset.IntensityCutoff == 0 || spec.SumIntensity < dataset.IntensityCutoff {
					dataset.IntensityCutoff = spec.SumIntensity
					*minIdx = len(dataset.Spectra) - 1
				}
			} else if spec.SumIntensity > dataset.IntensityCutoff { // if we have reached topN spectra, only add the spectrum if its intensity is above the cutoff
				// Pre-process spectrum.
				ScaleNormalizeBinSpectrum(&spec, par, peaks)
				// Delete the spectrum with the lowest intensity and add the new spectrum
				dataset.Spectra[*minIdx] = spec
				// Update intensity cutoff to the new lowest intensity in the topN
				// Find the spectrum with the lowest intensity in the current topN
				dataset.IntensityCutoff = dataset.Spectra[0].SumIntensity
				*minIdx = 0
				for i := 1; i < len(dataset.Spectra); i++ {
					if dataset.Spectra[i].SumIntensity < dataset.IntensityCutoff {
						dataset.IntensityCutoff = dataset.Spectra[i].SumIntensity
						*minIdx = i
					}
				}
			}
		} else { // No TopN filtering, just add the spectrum
			// Pre-process spectrum.
			// We do this here, since the data is probably still in cache
			ScaleNormalizeBinSpectrum(&spec, par, peaks)
			dataset.Spectra = append(dataset.Spectra, spec)
		}
	}
}

// ReadMGF Reads the MGF file and populates the dataset, with top-N filtering
func ReadMGF(par *Parameters, dataset *Dataset) error {
	if par.TwoPass {
		err := PreCheckMGF(par, dataset)
		if err != nil {
			return err
		}
	}

	file, err := os.Open(dataset.Filename)
	if err != nil {
		return fmt.Errorf("error opening dataset %c %s for reading: %v", dataset.ID, dataset.Filename, err)
	}
	defer file.Close()

	// Initialize dataset
	dataset.Spectra = []Spectrum{}

	scanner := bufio.NewScanner(file)
	specCount := 0
	minIdx := 0 // Index of the spectrum with the lowest intensity in the current topN, used for efficient replacement when a new spectrum with higher intensity is found

	specStatus := IN_RANGE
	var spec Spectrum
	var peaks []Peak
	spec.SpecIdx = 0
	nPeaks := 0

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		// Check if line is peak data
		if line[0] >= '0' && line[0] <= '9' {
			nPeaks++
			if specStatus == NOT_IN_RANGE {
				continue
			}
			parts := strings.Fields(line)
			if len(parts) >= 2 {
				mz, err1 := strconv.ParseFloat(parts[0], 64)
				intensity, err2 := strconv.ParseFloat(parts[1], 64)
				if err1 == nil && err2 == nil {
					if mz >= par.MinMz && mz <= par.MaxMz {
						// Append peak data to current spectrum
						peaks = append(peaks, Peak{Mz: mz, Intensity: intensity})
						// Compute sum of intensities for top-N filtering
						spec.SumIntensity += intensity
					}
					// For BasepeakIntensity and TotalIonCurrent, don't use mz filter
					// this is equal to the original C code
					if intensity > spec.BasepeakIntensity {
						spec.BasepeakIntensity = intensity
					}
					spec.TotalIonCurrent += intensity

					// FIXME: log errors in parsing mz or intensity, but continue processing
				}
			}
		} else if strings.HasPrefix(line, "BEGIN") {
			specStatus = IN_RANGE
		} else if strings.HasPrefix(line, "SCANS=") || strings.HasPrefix(line, "###MSMS:") {
			numStr := parseNumber(line)
			if scan, err := strconv.Atoi(numStr); err == nil {
				if scan < par.StartScan || scan > par.EndScan {
					specStatus = NOT_IN_RANGE
				}
				spec.Scan = scan
				dataset.ScanNumbersCouldBeRead = true
			}

		} else if strings.HasPrefix(line, "TITLE=") && strings.Contains(line, "scan=") {
			// Extract scan from TITLE
			scanIdx := strings.Index(line, "scan=")
			if scanIdx != -1 {
				numStr := parseNumber(line[scanIdx+5:])
				if scan, err := strconv.Atoi(numStr); err == nil {
					if scan < par.StartScan || scan > par.EndScan {
						specStatus = NOT_IN_RANGE
					}
					spec.Scan = scan
					dataset.ScanNumbersCouldBeRead = true
				}
			}
		} else if strings.HasPrefix(line, "RTINSECONDS=") {
			numStr := parseNumber(line)
			if rt, err := strconv.ParseFloat(numStr, 64); err == nil {
				if rt < par.StartRT || rt > par.EndRT {
					specStatus = NOT_IN_RANGE
				}
				spec.RT = rt
				dataset.RTsCouldBeRead = true
			}
		} else if strings.HasPrefix(line, "PEPMASS=") {
			numStr := parseNumber(line)
			if precursorMz, err := strconv.ParseFloat(numStr, 64); err == nil {
				spec.PrecursorMz = precursorMz
			}
		} else if strings.HasPrefix(line, "CHARGE=") {
			numStr := parseNumber(line)
			if charge, err := strconv.Atoi(numStr); err == nil {
				spec.Charge = charge
			}
		} else if strings.HasPrefix(line, "END") {
			if specStatus != NOT_IN_RANGE {
				StoreTopNSpectrum(par, dataset, spec, peaks, nPeaks, &minIdx)
			}
			// Start a new spectrum
			spec = Spectrum{}
			// Clear peaks, but keep the underlying array to avoid unnecessary allocations
			peaks = peaks[:0]
			nPeaks = 0
			specCount++
			spec.SpecIdx = specCount
			// Set default RT and Scan, in case they are not in the MGF file
			spec.RT = par.StartRT
			spec.Scan = par.StartScan
		}
	}
	if par.TopN > -1 && !par.TwoPass {
		// After TopN filtering, the spectra are not in the same order as in the file
		// Re-order spectra according to SpecIdx
		sort.Slice(dataset.Spectra, func(i, j int) bool {
			return dataset.Spectra[i].SpecIdx < dataset.Spectra[j].SpecIdx
		})
		fmt.Printf("done (ion intensity threshold %.3f)\n", dataset.IntensityCutoff)
	}
	dataset.SpectrumCount = specCount

	return scanner.Err()
}

// ScaleNormalizeBinSpectrum processes a single spectrum: scaling, normalizing, and binning
func ScaleNormalizeBinSpectrum(spec *Spectrum, par *Parameters, peaks []Peak) {
	spec.Bin = make([]float32, par.NBins)

	// Apply scaling and noise removal
	for k := range peaks {
		intensity := peaks[k].Intensity
		if intensity > par.Noise {
			peaks[k].Intensity = math.Pow(intensity, par.Scaling)
		} else {
			peaks[k].Intensity = 0
		}
	}

	// Normalize against spectral vector magnitude
	squareSum := 0.0
	for _, peak := range peaks {
		squareSum += peak.Intensity * peak.Intensity
	}
	rootSquareSum := math.Sqrt(squareSum)

	// Normalize and bin. Bins are auto-initialized to zero in Go, so we can skip that step.
	for _, peak := range peaks {
		intensity := peak.Intensity
		if rootSquareSum > 0 {
			intensity /= rootSquareSum
		}
		mz := peak.Mz
		if mz >= par.MinMz && mz < par.MaxMz {
			binIdx := int(math.Round((mz - par.MinMz) / par.BinSize))
			if binIdx >= 0 && binIdx < par.NBins {
				spec.Bin[binIdx] += float32(intensity)
			}
		}
	}

	// Normalize binned spectra
	squareSum = 0
	for _, bin := range spec.Bin {
		squareSum += float64(bin * bin)
	}
	rootSquareSum = math.Sqrt(squareSum)
	if rootSquareSum > 0 {
		for k := range spec.Bin {
			spec.Bin[k] /= float32(rootSquareSum)
		}
	}
	spec.nPeaks = len(peaks)
}

// ComputeDotProdHistogram compares spectra and computes histograms
func ComputeDotProdHistogram(par *Parameters, datasetA, datasetB *Dataset,
	dotprodHistogram []int, massDiffDotProductHistogram [][]int,
	nComparisons, greaterThanCutoff, sAB, actualCompared *int, asymmetric bool) {

	for i := 0; i < len(datasetA.Spectra); i++ {
		A := &datasetA.Spectra[i]

		if A.nPeaks < par.MinPeaks {
			continue
		}
		if A.BasepeakIntensity < par.MinBasepeakIntensity {
			continue
		}
		if A.TotalIonCurrent < par.MinTotalIonCurrent {
			continue
		}

		maxDotProd := 0.0
		(*actualCompared)++

		for j := 0; j < len(datasetB.Spectra); j++ {
			B := &datasetB.Spectra[j]

			if B.nPeaks < par.MinPeaks {
				continue
			}
			if B.BasepeakIntensity < par.MinBasepeakIntensity {
				continue
			}
			if B.TotalIonCurrent < par.MinTotalIonCurrent {
				continue
			}

			// Check scan number difference
			if math.Abs(float64(A.Scan-B.Scan)) > par.MaxScanNumberDifference {
				if B.Scan > A.Scan {
					break
				}
				continue
			}

			// Check RT difference
			if math.Abs(A.RT-B.RT) > par.MaxRTDifference {
				if B.RT > A.RT {
					break
				}
				continue
			}

			// Check precursor m/z difference
			if math.Abs(B.PrecursorMz-A.PrecursorMz) < par.MaxPrecursorDifference {
				// Compute dot product
				dotProd := 0.0
				for k := 0; k < par.NBins; k++ {
					dotProd += float64(A.Bin[k] * B.Bin[k])
				}

				if math.Abs(dotProd) <= 1.00 {
					// Apply spectrum metric if specified
					if par.SpectrumMetric == 1 {
						// Use spectral angle (SA)
						dotProd = 1 - 2*(math.Acos(dotProd)/math.Pi)
					}

					// Update histogram
					histIdx := dotprodHistogramBins/2 + int(math.Floor(dotProd*float64(dotprodHistogramBins/2-1e-9)))
					if histIdx >= 0 && histIdx < dotprodHistogramBins {
						dotprodHistogram[histIdx]++
					}

					// Update mass difference histogram (only for asymmetric comparison)
					if par.ExperimentalFeatures == 1 && asymmetric {
						massDiffBin := massdiffHistogramBins/2 + int(math.Floor((B.PrecursorMz-A.PrecursorMz)*massdiffHistogramBins/massdiffHistogramRange))
						if massDiffBin >= 0 && massDiffBin < massdiffHistogramBins {
							dotProdBin := dotprodHistogramBins/2 + int(math.Floor(dotProd*float64(dotprodHistogramBins/2-1e-9)))
							if dotProdBin >= 0 && dotProdBin < dotprodHistogramBins {
								massDiffDotProductHistogram[massDiffBin][dotProdBin]++
							}
						}
					}
				}

				(*nComparisons)++

				if dotProd > maxDotProd {
					maxDotProd = dotProd
				}
			}
		}

		if maxDotProd > par.Cutoff {
			(*greaterThanCutoff)++
			(*sAB)++
		}
	}
}

// JSONOutput represents the JSON output structure
type JSONOutput struct {
	CommandLine              string                 `json:"commandLine"`
	DatasetA                 string                 `json:"datasetA"`
	DatasetB                 string                 `json:"datasetB"`
	SetDistance              interface{}            `json:"setDistance"`
	SetMetric                int                    `json:"setMetric"`
	ScanRange                [2]int                 `json:"scanRange"`
	MaxScanDiff              float64                `json:"maxScanDiff"`
	MaxMzDiff                float64                `json:"maxMzDiff"`
	ScalingPower             float64                `json:"scalingPower"`
	NoiseThreshold           float64                `json:"noiseThreshold"`
	MinBasepeakIntensity     float64                `json:"minBasepeakIntensity"`
	MinTotalIonCurrent       float64                `json:"minTotalIonCurrent"`
	DatasetAQC               int                    `json:"datasetAQC"`
	DatasetBQC               int                    `json:"datasetBQC"`
	DatasetACompared         int                    `json:"datasetACompared"`
	DatasetBCompared         int                    `json:"datasetBCompared"`
	NrGtCutoff               int                    `json:"nrGtCutoff"`
	NrComparisons            int                    `json:"nrComparisons"`
	MinPeaks                 int                    `json:"minPeaks"`
	MaxPeaks                 int                    `json:"maxPeaks"`
	MzRange                  [2]float64             `json:"mzRange"`
	MzBinSize                float64                `json:"mzBinSize"`
	NrMzBins                 int                    `json:"nrMzBins"`
	DotProdHistogram         HistogramData          `json:"dotProdHistogram"`
	MassDiffDotProdHistogram *MassDiffHistogramData `json:"massDiffDotProdHistogram,omitempty"`
}

type HistogramData struct {
	NrBins       int        `json:"nrBins"`
	DotProdRange [2]float64 `json:"dotProdRange"`
	Count        []int      `json:"count"`
}

type MassDiffHistogramData struct {
	MzRange      [2]float64 `json:"mzRange"`
	DotProdRange [2]float64 `json:"dotProdRange"`
	Count        [][]int    `json:"count"`
}

// WriteJSON writes results to JSON file
func WriteJSON(par *Parameters, datasetA, datasetB *Dataset,
	dotprodHistogram []int, massDiffDotProductHistogram [][]int,
	nComparisons, greaterThanCutoff, sAB, sBA, datasetACompared, datasetBCompared int) error {

	if par.JSONFilename == "" {
		return nil
	}

	// Calculate set distance
	var setDistance interface{}
	sizeA := len(datasetA.Spectra)
	sizeB := len(datasetB.Spectra)

	switch par.Metric {
	case 0: // original metric
		if greaterThanCutoff > 0 {
			setDistance = float64(nComparisons) / 2.0 / float64(greaterThanCutoff)
		} else {
			setDistance = "INF"
		}
	case 1: // symmetric metric
		if greaterThanCutoff > 0 {
			setDistance = float64(sizeA+sizeB) / float64(greaterThanCutoff)
		} else {
			setDistance = "INF"
		}
	case 2: // compareMS2 2.0 symmetric metric
		if sizeA == 0 || sizeB == 0 {
			setDistance = "INF"
		} else if (sAB + sBA) > 0 {
			setDistance = 1.0/(float64(sAB)/float64(2*sizeA)+float64(sBA)/float64(2*sizeB)) - 1.0
		} else {
			setDistance = (4.0*float64(sizeA)*float64(sizeB))/float64(sizeA+sizeB) - 1.0
		}
	}

	output := JSONOutput{
		CommandLine:          strings.Join(os.Args, " "),
		DatasetA:             datasetA.Filename,
		DatasetB:             datasetB.Filename,
		SetDistance:          setDistance,
		SetMetric:            par.Metric,
		ScanRange:            [2]int{par.StartScan, par.EndScan},
		MaxScanDiff:          par.MaxScanNumberDifference,
		MaxMzDiff:            par.MaxPrecursorDifference,
		ScalingPower:         par.Scaling,
		NoiseThreshold:       par.Noise,
		MinBasepeakIntensity: par.MinBasepeakIntensity,
		MinTotalIonCurrent:   par.MinTotalIonCurrent,
		DatasetACompared:     datasetACompared,
		DatasetBCompared:     datasetBCompared,
		NrGtCutoff:           greaterThanCutoff,
		NrComparisons:        nComparisons,
		MinPeaks:             par.MinPeaks,
		MaxPeaks:             par.PeakCount,
		MzRange:              [2]float64{par.MinMz, par.MaxMz},
		MzBinSize:            par.BinSize,
		NrMzBins:             par.NBins,
		DotProdHistogram: HistogramData{
			NrBins:       dotprodHistogramBins,
			DotProdRange: [2]float64{-1.0, 1.0},
			Count:        dotprodHistogram,
		},
	}

	if par.QC == 0 {
		output.DatasetAQC = datasetA.SpectrumCount
		output.DatasetBQC = datasetB.SpectrumCount
	}

	if par.ExperimentalFeatures == 1 {
		output.MassDiffDotProdHistogram = &MassDiffHistogramData{
			MzRange:      [2]float64{-massdiffHistogramRange / 2, massdiffHistogramRange / 2},
			DotProdRange: [2]float64{-1.0, 1.0},
			Count:        massDiffDotProductHistogram,
		}
	}

	file, err := os.Create(par.JSONFilename)
	if err != nil {
		return fmt.Errorf("error opening JSON output file %s: %v", par.JSONFilename, err)
	}
	defer file.Close()

	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "\t")
	return encoder.Encode(output)
}

// WriteTextOutput writes results to text file
func WriteTextOutput(par *Parameters, datasetA, datasetB *Dataset,
	dotprodHistogram, massDiffHistogram []int,
	nComparisons, greaterThanCutoff, sAB, sBA int) error {

	file, err := os.Create(par.OutputFilename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", par.OutputFilename, err)
	}
	defer file.Close()

	fmt.Fprintf(file, "dataset_A\t%s\n", datasetA.Filename)
	fmt.Fprintf(file, "dataset_B\t%s\n", datasetB.Filename)

	sizeA := len(datasetA.Spectra)
	sizeB := len(datasetB.Spectra)

	// Calculate and write set distance
	switch par.Metric {
	case 0: // original metric
		if greaterThanCutoff > 0 {
			fmt.Fprintf(file, "set_distance\t%.10f\n", float64(nComparisons)/2.0/float64(greaterThanCutoff))
		} else {
			fmt.Fprintf(file, "set_distance\tINF\n")
		}
	case 1: // symmetric metric
		if greaterThanCutoff > 0 {
			fmt.Fprintf(file, "set_distance\t%.10f\n", float64(sizeA+sizeB)/float64(greaterThanCutoff))
		} else {
			fmt.Fprintf(file, "set_distance\tINF\n")
		}
	case 2: // compareMS2 2.0 symmetric metric
		if sizeA == 0 || sizeB == 0 {
			fmt.Fprintf(file, "set_distance\tINF\n")
		} else if (sAB + sBA) > 0 {
			distance := 1.0/(float64(sAB)/float64(2*sizeA)+float64(sBA)/float64(2*sizeB)) - 1.0
			fmt.Fprintf(file, "set_distance\t%.10f\n", distance)
		} else {
			distance := (4.0*float64(sizeA)*float64(sizeB))/float64(sizeA+sizeB) - 1.0
			fmt.Fprintf(file, "set_distance\t%.10f\n", distance)
		}
	}

	fmt.Fprintf(file, "set_metric\t%d\n", par.Metric)
	fmt.Fprintf(file, "scan_range\t%d\t%d\n", par.StartScan, par.EndScan)
	fmt.Fprintf(file, "max_scan_diff\t%.5f\n", par.MaxScanNumberDifference)
	fmt.Fprintf(file, "max_m/z_diff\t%.5f\n", par.MaxPrecursorDifference)
	fmt.Fprintf(file, "scaling_power\t%.5f\n", par.Scaling)
	fmt.Fprintf(file, "noise_threshold\t%.5f\n", par.Noise)
	fmt.Fprintf(file, "min_basepeak_intensity\t%.2f\n", par.MinBasepeakIntensity)
	fmt.Fprintf(file, "min_total_ion_current\t%.2f\n", par.MinTotalIonCurrent)

	if par.QC == 0 {
		fmt.Fprintf(file, "dataset_A_QC\t%.4f\n", float64(sizeA))
		fmt.Fprintf(file, "dataset_B_QC\t%.4f\n", float64(sizeB))
	}

	fmt.Fprintf(file, "n_gt_cutoff\t%d\n", greaterThanCutoff)
	fmt.Fprintf(file, "n_comparisons\t%d\n", nComparisons)
	fmt.Fprintf(file, "min_peaks\t%d\n", par.MinPeaks)
	fmt.Fprintf(file, "m/z_range\t%.4f\t%.4f\n", par.MinMz, par.MaxMz)
	fmt.Fprintf(file, "m/z_bin_size\t%.4f\n", par.BinSize)
	fmt.Fprintf(file, "n_m/z_bins\t%d\n", par.NBins)

	// Write histogram
	for i := 0; i < dotprodHistogramBins; i++ {
		fmt.Fprintf(file, "histogram\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
			float64(i-100)/100.0,
			float64(i+1-100)/100.0,
			(float64(i-100)+0.5)/100.0,
			dotprodHistogram[i],
			massDiffHistogram[i])
	}

	return nil
}

// WriteExperimentalOutput writes experimental histogram data
func WriteExperimentalOutput(par *Parameters, massDiffDotProductHistogram [][]int) error {
	if par.ExperimentalFeatures != 1 {
		return nil
	}

	file, err := os.Create(par.ExperimentalOutputFilename)
	if err != nil {
		return fmt.Errorf("error opening experimental output file %s: %v", par.ExperimentalOutputFilename, err)
	}
	defer file.Close()

	for i := 0; i < dotprodHistogramBins; i++ {
		for j := 0; j < massdiffHistogramBins-1; j++ {
			fmt.Fprintf(file, "%d\t", massDiffDotProductHistogram[j][i])
		}
		fmt.Fprintf(file, "%d\n", massDiffDotProductHistogram[massdiffHistogramBins-1][i])
	}

	return nil
}

func main() {
	// Parse parameters
	par := NewParameters()
	if err := par.ParseFlags(); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		flag.Usage()
		os.Exit(1)
	}

	if par.CPUProfile != "" {
		fmt.Printf("writing CPU profile. Run go tool pprof compareMS2 %s view\n", par.CPUProfile)
		f, err := os.Create(par.CPUProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// Calculate number of bins
	if math.Floor((par.MaxMz-par.MinMz)/par.BinSize) == (par.MaxMz-par.MinMz)/par.BinSize {
		par.NBins = int((par.MaxMz - par.MinMz) / par.BinSize)
	} else {
		par.NBins = int((par.MaxMz-par.MinMz)/par.BinSize) + 1
	}

	fmt.Printf("spectrum bin size %.3f Th -> %d bins in [%.3f,%.3f]\n",
		par.BinSize, par.NBins, par.MinMz, par.MaxMz)

	fmt.Printf("scan range=[%d,%d]\n", par.StartScan, par.EndScan)
	fmt.Printf("max scan difference=%.2f\n", par.MaxScanNumberDifference)
	fmt.Printf("max m/z difference=%.4f\n", par.MaxPrecursorDifference)
	fmt.Printf("scaling=^%.2f\n", par.Scaling)
	fmt.Printf("noise=%.1f\n", par.Noise)
	fmt.Printf("min basepeak intensity=%.2f\n", par.MinBasepeakIntensity)
	fmt.Printf("min total ion current=%.2f\n", par.MinTotalIonCurrent)

	// Initialize datasets
	datasetA := &Dataset{ID: 'A', Filename: par.DatasetAFilename}
	datasetB := &Dataset{ID: 'B', Filename: par.DatasetBFilename}

	// Read MGF files
	if err := ReadMGF(par, datasetA); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	if err := ReadMGF(par, datasetB); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
	// Check if scan numbers and RTs could be read
	if !datasetA.ScanNumbersCouldBeRead {
		fmt.Printf("warning: scan numbers could not be read from dataset A (%s)\n", datasetA.Filename)
	}
	if !datasetB.ScanNumbersCouldBeRead {
		fmt.Printf("warning: scan numbers could not be read from dataset B (%s)\n", datasetB.Filename)
	}
	if !datasetA.ScanNumbersCouldBeRead || !datasetB.ScanNumbersCouldBeRead {
		par.MaxScanNumberDifference = 10000
		fmt.Println("scan filters will be ignored")
	}

	if !datasetA.RTsCouldBeRead {
		fmt.Printf("warning: retention times could not be read from dataset A (%s)\n", datasetA.Filename)
	}
	if !datasetB.RTsCouldBeRead {
		fmt.Printf("warning: retention times could not be read from dataset B (%s)\n", datasetB.Filename)
	}
	if !datasetA.RTsCouldBeRead || !datasetB.RTsCouldBeRead {
		par.MaxRTDifference = 60
		fmt.Println("retention time filters will be ignored")
	}

	// Initialize histograms
	fmt.Println("matching spectra and computing set distance.")
	dotprodHistogram := make([]int, dotprodHistogramBins)
	massDiffHistogram := make([]int, dotprodHistogramBins)

	var massDiffDotProductHistogram [][]int
	if par.ExperimentalFeatures == 1 {
		massDiffDotProductHistogram = make([][]int, massdiffHistogramBins)
		for i := 0; i < massdiffHistogramBins; i++ {
			massDiffDotProductHistogram[i] = make([]int, dotprodHistogramBins)
		}
	}

	var nComparisons, greaterThanCutoff, sAB, sBA int
	var datasetAActualCompared, datasetBActualCompared int

	// Compare A to B
	ComputeDotProdHistogram(par, datasetA, datasetB, dotprodHistogram,
		massDiffDotProductHistogram, &nComparisons, &greaterThanCutoff,
		&sAB, &datasetAActualCompared, true)
	fmt.Print(".")

	// Compare B to A
	ComputeDotProdHistogram(par, datasetB, datasetA, dotprodHistogram,
		massDiffDotProductHistogram, &nComparisons, &greaterThanCutoff,
		&sBA, &datasetBActualCompared, false)
	fmt.Printf(".done (compared %d (|S_AB|=%d) spectra from dataset A with %d (|S_BA|=%d) spectra from dataset B)\n",
		datasetAActualCompared, sAB, datasetBActualCompared, sBA)

	// Write results
	fmt.Println("writing results to file...")

	if err := WriteTextOutput(par, datasetA, datasetB, dotprodHistogram, massDiffHistogram,
		nComparisons, greaterThanCutoff, sAB, sBA); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	if err := WriteExperimentalOutput(par, massDiffDotProductHistogram); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	if err := WriteJSON(par, datasetA, datasetB, dotprodHistogram, massDiffDotProductHistogram,
		nComparisons, greaterThanCutoff, sAB, sBA, datasetAActualCompared, datasetBActualCompared); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	fmt.Println("done")
}
