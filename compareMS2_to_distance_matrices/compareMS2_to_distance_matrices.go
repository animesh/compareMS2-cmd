package main

// SPDX-License-Identifier: MIT
// Copyright 2022 Magnus Palmblad
// Copyright 2022 Rob Marissen.
// Go conversion with help from GitHub Copilot with Claude Sonnet 4.5

// compareMS2_to_distance_matrices - collects output from compareMS2 to distance matrices
// Contact information: n.m.palmblad@lumc.nl

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"unicode"
)

// Species contains all species data
type Species struct {
	name             []string       // species names indexed by species index
	speciesHasData   map[int]bool   // tracks which species have data
	nameToIdx        map[string]int // species name -> index in name
	sampleToIdx      map[string]int // sample name -> species index
	sampleLowerToIdx map[string]int // lowercased sample name -> species index
}

// DistanceData holds distance values and sample counts
type DistanceData struct {
	Value   float64
	Samples int
}

// QCData holds QC values and sample counts
type QCData struct {
	Value   float64
	Samples int
}

// Average returns the average QC value
func (qc QCData) Average() float64 {
	if qc.Samples > 0 {
		return qc.Value / float64(qc.Samples)
	}
	return 0.0
}

// DistanceMatrix is a 2D map storing distance data between species pairs
type DistanceMatrix map[int]map[int]DistanceData

// ComparisonResult for JSON parsing
type ComparisonResult struct {
	DatasetA    string  `json:"datasetA"`
	DatasetB    string  `json:"datasetB"`
	SetDistance float64 `json:"setDistance"`
	DatasetAQC  float64 `json:"datasetAQC"`
	DatasetBQC  float64 `json:"datasetBQC"`
}

// DistanceMatrixOutput for JSON output
type DistanceMatrixOutput struct {
	DistanceMatrix struct {
		Title     string      `json:"title"`
		Cutoff    float64     `json:"cutoff"`
		Species   []SpeciesQC `json:"species"`
		Distances [][]float64 `json:"distances"`
	} `json:"distanceMatrix"`
}

// SpeciesQC for JSON output
type SpeciesQC struct {
	Name string  `json:"name"`
	QC   float64 `json:"qc"`
}

// Add a species to the species list
func (s *Species) addSpecies(speciesName string) {
	idx := len(s.name)
	s.name = append(s.name, speciesName)
	s.nameToIdx[speciesName] = idx
}

// Find species index by name, add if not found
func (s *Species) speciesNameToIdx(speciesName string) int {
	if idx, ok := s.nameToIdx[speciesName]; ok {
		return idx
	}
	s.addSpecies(speciesName)
	return len(s.name) - 1
}

// sampleBaseName returns the filename segment for either POSIX or Windows paths.
func sampleBaseName(name string) string {
	idx := strings.LastIndexAny(name, `/\`)
	if idx >= 0 {
		return name[idx+1:]
	}
	return name
}

// Add a sample to a species
func (s *Species) addSampleSpecies(sampleName, speciesName string) int {
	idx := s.speciesNameToIdx(speciesName)
	s.sampleToIdx[sampleName] = idx
	if s.sampleLowerToIdx == nil {
		s.sampleLowerToIdx = make(map[string]int)
	}
	s.sampleLowerToIdx[strings.ToLower(sampleName)] = idx
	return idx
}

// Map sample name to species index
func (s *Species) sampleNameToSpeciesIndex(sampleName string) int {
	base := sampleBaseName(sampleName)

	// Check exact basename first, then case-insensitive fallback
	if idx, ok := s.sampleToIdx[base]; ok {
		s.speciesHasData[idx] = true
		return idx
	}
	if s.sampleLowerToIdx != nil {
		if idx, ok := s.sampleLowerToIdx[strings.ToLower(base)]; ok {
			s.speciesHasData[idx] = true
			return idx
		}
	}

	// Not found, create new species from sample name
	speciesName := strings.TrimSuffix(base, filepath.Ext(base))
	idx := s.addSampleSpecies(base, speciesName)
	s.speciesHasData[idx] = true
	return idx
}

// Read sample-to-species mapping file
func readSampleToSpecies(filename string, species *Species) error {
	data, err := os.ReadFile(filename)
	if err != nil {
		return fmt.Errorf("error opening sample-species mapping file %s: %v", filename, err)
	}

	content := strings.TrimPrefix(string(data), "\uFEFF")
	content = strings.ReplaceAll(content, "\r\n", "\n")
	content = strings.ReplaceAll(content, "\r", "\n")

	for _, line := range strings.Split(content, "\n") {
		if strings.TrimSpace(line) == "" {
			continue
		}

		var parts []string
		for _, field := range strings.Split(line, "\t") {
			if trimmed := strings.TrimSpace(field); trimmed != "" {
				parts = append(parts, trimmed)
			}
		}
		if len(parts) >= 2 {
			sampleName := sampleBaseName(parts[0])
			speciesName := parts[1]
			if sampleName != "" && speciesName != "" {
				species.addSampleSpecies(sampleName, speciesName)
			}
		}
	}

	return nil
}

// Read file list
func readFileList(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening filename list %s: %v", filename, err)
	}
	defer file.Close()

	var files []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			files = append(files, line)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return files, nil
}

// Get indices of used species
func (s *Species) usedIndices() []int {
	var indices []int
	for i := range s.name {
		if s.speciesHasData[i] {
			indices = append(indices, i)
		}
	}
	return indices
}

// Create NEXUS format output
func createNexus(outputStem string, distances DistanceMatrix, species *Species, metric int) error {
	filename := outputStem + ".nexus"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in NEXUS format...")

	used := species.usedIndices()
	nSpecies := len(used)

	fmt.Fprintln(file, "#NEXUS")
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintln(file, "[ NEXUS file generated by compareMS2_to_distance_matrices from compareMS2 data ]")
	fmt.Fprintln(file, "[ compareMS2 and compareMS2_to_distance_matrices are (c) Magnus Palmblad 2010- ]")
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintln(file, "[ note that compareMS2_to_distance_matrices outputs a full matrix and that the ]")
	fmt.Fprintln(file, "[ upper (right) and lower (left) triangular matrices can be slightly different ]")
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintf(file, "[ distances are based on metric %d                                              ]\n", metric)
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintln(file, "BEGIN taxa;")
	fmt.Fprintf(file, "   DIMENSIONS ntax=%d;\n", nSpecies)
	fmt.Fprint(file, "   TAXLABELS")

	for _, i := range used {
		fmt.Fprintf(file, "   %s", species.name[i])
	}
	fmt.Fprintln(file, ";")
	fmt.Fprintln(file, "END;")
	fmt.Fprintln(file, "")
	fmt.Fprintln(file, "BEGIN distances;")
	fmt.Fprintf(file, "   DIMENSIONS ntax=%d;\n", nSpecies)
	fmt.Fprintln(file, "   FORMAT")
	fmt.Fprintln(file, "       triangle=both")
	fmt.Fprintln(file, "       labels=left")
	fmt.Fprintln(file, "     diagonal")
	fmt.Fprintln(file, "     missing=?")
	fmt.Fprintln(file, "   ;")
	fmt.Fprintln(file, "MATRIX")

	for _, i := range used {
		fmt.Fprintf(file, "%s\t", species.name[i])
		for _, j := range used {
			if i == j {
				fmt.Fprintf(file, "%1.5f\t", 0.0)
			} else if i < j {
				fmt.Fprintf(file, "%1.5f\t", distances[i][j].Value)
			} else {
				fmt.Fprintf(file, "%1.5f\t", distances[j][i].Value)
			}
		}
		fmt.Fprintln(file, "")
	}

	fmt.Fprintln(file, ";")
	fmt.Fprintln(file, "END;")

	return nil
}

// Create MEGA format (version < 12)
func createMega(outputStem string, distances DistanceMatrix, species *Species, cutoff float64, qcData map[int]QCData) error {
	filename := outputStem + ".meg"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in MEGA version < 12 format...")

	fmt.Fprint(file, "#mega\r\n")
	fmt.Fprintf(file, "TITLE: %s (lower-left triangular matrix, cutoff=%.4f)\r\n\r\n", filename, cutoff)

	used := species.usedIndices()

	// Output QC values
	for _, i := range used {
		fmt.Fprintf(file, "QC\t%s\t%.3f\r\n", species.name[i], qcData[i].Average())
	}
	fmt.Fprint(file, "\r\n")

	// Output species names
	for _, i := range used {
		fmt.Fprintf(file, "#%s\r\n", species.name[i])
	}
	fmt.Fprint(file, "\r\n\r\n")

	// Output distance matrix (lower triangular)
	for ui := 1; ui < len(used); ui++ {
		y := used[ui]
		for ux := 0; ux < ui; ux++ {
			fmt.Fprintf(file, "%1.5f\t", distances[used[ux]][y].Value)
		}
		fmt.Fprint(file, "\r\n")
	}

	return nil
}

// Fix taxonomy name for MEGA format
func fixTaxonomyName(name string) string {
	const firstAllowed = "-+."
	const restAllowed = "-+._*:()|\\//"

	runes := []rune(name)
	for i, r := range runes {
		if unicode.IsLetter(r) || unicode.IsDigit(r) {
			continue
		}
		allowed := restAllowed
		if i == 0 {
			allowed = firstAllowed
		}
		if !strings.ContainsRune(allowed, r) {
			runes[i] = '_'
		}
	}
	return string(runes)
}

// Create MEGA format version 12
func createMega12(outputStem string, distances DistanceMatrix, species *Species, cutoff float64) error {
	filename := outputStem + ".meg"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in MEGA 12 format...")

	used := species.usedIndices()
	nSpecies := len(used)

	fmt.Fprintln(file, "#mega")
	fmt.Fprintf(file, "!Title: %s;\n", outputStem)
	fmt.Fprintf(file, "!Format DataType=Distance DataFormat=LowerLeft NTaxa=%d;\n", nSpecies)
	fmt.Fprintln(file, "!Description")
	fmt.Fprintf(file, "  No. of Groups : %d\n", nSpecies)
	fmt.Fprintf(file, "  Cutoff : %.4f\n", cutoff)
	fmt.Fprintln(file, ";")
	fmt.Fprintln(file, "")

	// Output numbered species list
	for idx, i := range used {
		fmt.Fprintf(file, "[%d] #%s\n", idx+1, fixTaxonomyName(species.name[i]))
	}
	fmt.Fprintln(file, "")

	// Output column header
	fmt.Fprint(file, "[")
	for idx := range used {
		fmt.Fprintf(file, "%9d ", idx+1)
	}
	fmt.Fprintln(file, " ]")

	// Output distance matrix
	for yIdx, y := range used {
		fmt.Fprintf(file, "[%d] ", yIdx+1)
		for xIdx, x := range used {
			if xIdx < yIdx {
				fmt.Fprintf(file, "%9.5f ", distances[x][y].Value)
			} else if xIdx == yIdx {
				fmt.Fprint(file, "          ")
			}
		}
		fmt.Fprintln(file, "")
	}
	fmt.Fprintln(file, "")

	return nil
}

// Create JSON format output
func createJSON(outputStem string, distances DistanceMatrix, species *Species, cutoff float64, qcData map[int]QCData) error {
	filename := outputStem + ".json"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in JSON format...")

	var output DistanceMatrixOutput
	output.DistanceMatrix.Title = outputStem
	output.DistanceMatrix.Cutoff = cutoff

	used := species.usedIndices()

	// Build species list with QC values
	for _, i := range used {
		output.DistanceMatrix.Species = append(output.DistanceMatrix.Species, SpeciesQC{
			Name: species.name[i],
			QC:   qcData[i].Average(),
		})
	}

	// Build distance matrix (lower triangular)
	for ui := 1; ui < len(used); ui++ {
		y := used[ui]
		row := make([]float64, 0, ui)
		for ux := 0; ux < ui; ux++ {
			row = append(row, distances[used[ux]][y].Value)
		}
		output.DistanceMatrix.Distances = append(output.DistanceMatrix.Distances, row)
	}

	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "  ")
	if err := encoder.Encode(output); err != nil {
		return fmt.Errorf("error encoding JSON: %v", err)
	}

	return nil
}

// updateDistanceAndQC updates distance matrix and QC data for a comparison
func updateDistanceAndQC(species *Species, distances DistanceMatrix, qcData map[int]QCData,
	datasetA, datasetB string, setDistance, qcA, qcB float64) {

	a := species.sampleNameToSpeciesIndex(datasetA)
	b := species.sampleNameToSpeciesIndex(datasetB)

	// Ensure lower triangular: x < y
	x, y := a, b
	if a > b {
		x, y = b, a
	}

	// Update distance
	if distances[x] == nil {
		distances[x] = make(map[int]DistanceData)
	}
	d := distances[x][y]
	d.Value += setDistance
	d.Samples++
	distances[x][y] = d

	// Update QC for dataset A
	qcDataA := qcData[a]
	qcDataA.Value += qcA
	qcDataA.Samples++
	qcData[a] = qcDataA

	// Update QC for dataset B
	qcDataB := qcData[b]
	qcDataB.Value += qcB
	qcDataB.Samples++
	qcData[b] = qcDataB
}

// Parse JSON comparison file
func parseComparisonJSON(filename string, compIdx int, species *Species, distances DistanceMatrix, qcData map[int]QCData) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	fmt.Printf("\nreading pairwise comparison %d (JSON)", compIdx)

	var result ComparisonResult
	decoder := json.NewDecoder(file)
	if err := decoder.Decode(&result); err != nil {
		return fmt.Errorf("error decoding JSON file %s: %v", filename, err)
	}

	updateDistanceAndQC(species, distances, qcData,
		result.DatasetA, result.DatasetB, result.SetDistance, result.DatasetAQC, result.DatasetBQC)

	return nil
}

// Parse text comparison file
func parseComparisonTxt(filename string, compIdx int, species *Species, distances DistanceMatrix, qcData map[int]QCData) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	var result ComparisonResult

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		parts := strings.Split(line, "\t")
		if len(parts) < 2 {
			continue
		}

		field := parts[0]
		value := parts[1]

		switch field {
		case "dataset_A":
			result.DatasetA = value
		case "dataset_B":
			result.DatasetB = value
		case "set_distance":
			result.SetDistance, _ = strconv.ParseFloat(value, 64)
		case "dataset_A_QC":
			result.DatasetAQC, _ = strconv.ParseFloat(value, 64)
		case "dataset_B_QC":
			result.DatasetBQC, _ = strconv.ParseFloat(value, 64)
		}
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	fmt.Printf("\nread pairwise comparison %d (%s and %s)", compIdx, result.DatasetA, result.DatasetB)

	updateDistanceAndQC(species, distances, qcData,
		result.DatasetA, result.DatasetB, result.SetDistance, result.DatasetAQC, result.DatasetBQC)

	return nil
}

func usage() {
	fmt.Println("usage: compareMS2_to_distance_matrices -i <input file> -o <output stem> [options]")
	fmt.Println()
	fmt.Println("Required arguments:")
	fmt.Println("  -i <file>    List of compareMS2 results files to process")
	fmt.Println("  -o <stem>    Output file stem (output filename without extension)")
	fmt.Println()
	fmt.Println("Optional arguments:")
	fmt.Println("  -x <file>    Sample to species mapping file")
	fmt.Println("  -c <value>   Score cutoff (default: 0.80)")
	fmt.Println("  -n           Output in NEXUS format (default if no format specified)")
	fmt.Println("  -m           Output in MEGA format (version < 12)")
	fmt.Println("  -m2          Output in MEGA format (version 12+)")
	fmt.Println("  -J           Output in JSON format")
	fmt.Println("  -h           Display this help message")
	fmt.Println()
	fmt.Println("Note: Multiple format flags can be combined to generate multiple outputs.")
	fmt.Println()
	fmt.Println("Examples:")
	fmt.Println("  compareMS2_to_distance_matrices -i filelist.txt -o results -m2")
	fmt.Println("  compareMS2_to_distance_matrices -i filelist.txt -o results -m2 -J")
}

func main() {
	// Define flags
	inputFile := flag.String("i", "", "List of compareMS2 results files to process")
	outputStem := flag.String("o", "", "Output file stem (output filename without extension)")
	mappingFile := flag.String("x", "", "Sample to species mapping file")
	cutoff := flag.Float64("c", 0.80, "Score cutoff")
	formatNexus := flag.Bool("n", false, "Output in NEXUS format")
	formatMega := flag.Bool("m", false, "Output in MEGA format (version < 12)")
	formatMega12 := flag.Bool("m2", false, "Output in MEGA format (version 12+)")
	formatJSON := flag.Bool("J", false, "Output in JSON format")
	help := flag.Bool("h", false, "Display help message")

	flag.Parse()

	if *help {
		fmt.Println("compareMS2_to_distance_matrices - (c) Magnus Palmblad")
		fmt.Println()
		usage()
		os.Exit(0)
	}

	if *inputFile == "" || *outputStem == "" {
		usage()
		os.Exit(1)
	}

	// Determine output formats
	hasFormat := *formatNexus || *formatMega || *formatMega12 || *formatJSON
	if !hasFormat {
		// Default to NEXUS if no format specified
		*formatNexus = true
	}

	fmt.Print("reading list of compareMS2 results files...")

	// Read comparison file list
	comparisons, err := readFileList(*inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "\n%v\n", err)
		os.Exit(1)
	}

	species := &Species{
		name:             make([]string, 0),
		speciesHasData:   make(map[int]bool),
		nameToIdx:        make(map[string]int),
		sampleToIdx:      make(map[string]int),
		sampleLowerToIdx: make(map[string]int),
	}

	// Read sample-to-species mapping if provided
	if *mappingFile != "" {
		fmt.Print("done\nreading sample to species file...")
		if err := readSampleToSpecies(*mappingFile, species); err != nil {
			fmt.Fprintf(os.Stderr, "\n%v\n", err)
			os.Exit(1)
		}
	}

	fmt.Println("done")

	// Initialize data structures using maps
	distances := make(DistanceMatrix)
	qcData := make(map[int]QCData)

	// Process comparison files
	for i, compFile := range comparisons {
		isJSON := strings.HasSuffix(compFile, ".json")

		var err error
		if isJSON {
			err = parseComparisonJSON(compFile, i+1, species, distances, qcData)
		} else {
			err = parseComparisonTxt(compFile, i+1, species, distances, qcData)
		}

		if err != nil {
			fmt.Fprintf(os.Stderr, "\nerror opening MS2compare results file %s: %v\n", compFile, err)
			os.Exit(1)
		}
	}

	fmt.Printf("\ndone (read %d pairwise comparisons)\n", len(comparisons))

	// Compute averages for distances with multiple samples
	for x, row := range distances {
		for y, d := range row {
			if d.Samples > 0 {
				d.Value /= float64(d.Samples)
				distances[x][y] = d
			}
		}
	}

	// Generate output in requested formats
	metric := 0 // Not used in current implementation

	if *formatNexus {
		if err := createNexus(*outputStem, distances, species, metric); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}

	if *formatMega {
		if err := createMega(*outputStem, distances, species, *cutoff, qcData); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}

	if *formatMega12 {
		// When both MEGA and MEGA12 are requested, add suffix to avoid overwriting
		stem := *outputStem
		if *formatMega {
			stem = *outputStem + "_MEGA_12_"
		}
		if err := createMega12(stem, distances, species, *cutoff); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}

	if *formatJSON {
		if err := createJSON(*outputStem, distances, species, *cutoff, qcData); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}
}
