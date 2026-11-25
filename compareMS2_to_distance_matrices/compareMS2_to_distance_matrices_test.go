package main

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func newTestSpecies() *Species {
	return &Species{
		name:           make([]string, 0),
		speciesHasData: make(map[int]bool),
		nameToIdx:      make(map[string]int),
		sampleToIdx:    make(map[string]int),
	}
}

func averageDistances(distances DistanceMatrix) {
	for x, row := range distances {
		for y, d := range row {
			if d.Samples > 0 {
				d.Value /= float64(d.Samples)
				distances[x][y] = d
			}
		}
	}
}

func TestParseComparisonJSONAndTxt_MatchesAndAverages(t *testing.T) {
	dir := t.TempDir()

	mapPath := filepath.Join(dir, "species_map.txt")
	mapContent := "s1a.mgf\tSpecies_1\ns1b.mgf\tSpecies_1\ns2.mgf\tSpecies_2\n"
	if err := os.WriteFile(mapPath, []byte(mapContent), 0644); err != nil {
		t.Fatalf("failed writing mapping file: %v", err)
	}

	json1 := filepath.Join(dir, "c1.json")
	json2 := filepath.Join(dir, "c2.json")
	if err := os.WriteFile(json1, []byte(`{"datasetA":"/tmp/s1a.mgf","datasetB":"/tmp/s2.mgf","setDistance":1.0,"datasetAQC":10,"datasetBQC":20}`), 0644); err != nil {
		t.Fatalf("failed writing json1: %v", err)
	}
	if err := os.WriteFile(json2, []byte(`{"datasetA":"/tmp/s1b.mgf","datasetB":"/tmp/s2.mgf","setDistance":3.0,"datasetAQC":14,"datasetBQC":22}`), 0644); err != nil {
		t.Fatalf("failed writing json2: %v", err)
	}

	txt1 := filepath.Join(dir, "c1.txt")
	txt2 := filepath.Join(dir, "c2.txt")
	if err := os.WriteFile(txt1, []byte("dataset_A\t/tmp/s1a.mgf\ndataset_B\t/tmp/s2.mgf\nset_distance\t1.0\ndataset_A_QC\t10.0\ndataset_B_QC\t20.0\n"), 0644); err != nil {
		t.Fatalf("failed writing txt1: %v", err)
	}
	if err := os.WriteFile(txt2, []byte("dataset_A\t/tmp/s1b.mgf\ndataset_B\t/tmp/s2.mgf\nset_distance\t3.0\ndataset_A_QC\t14.0\ndataset_B_QC\t22.0\n"), 0644); err != nil {
		t.Fatalf("failed writing txt2: %v", err)
	}

	spJSON := newTestSpecies()
	if err := readSampleToSpecies(mapPath, spJSON); err != nil {
		t.Fatalf("readSampleToSpecies(json) failed: %v", err)
	}
	distJSON := make(DistanceMatrix)
	qcJSON := make(map[int]QCData)
	for i, f := range []string{json1, json2} {
		if err := parseComparisonJSON(f, i+1, spJSON, distJSON, qcJSON); err != nil {
			t.Fatalf("parseComparisonJSON(%s) failed: %v", f, err)
		}
	}
	averageDistances(distJSON)

	spTxt := newTestSpecies()
	if err := readSampleToSpecies(mapPath, spTxt); err != nil {
		t.Fatalf("readSampleToSpecies(txt) failed: %v", err)
	}
	distTxt := make(DistanceMatrix)
	qcTxt := make(map[int]QCData)
	for i, f := range []string{txt1, txt2} {
		if err := parseComparisonTxt(f, i+1, spTxt, distTxt, qcTxt); err != nil {
			t.Fatalf("parseComparisonTxt(%s) failed: %v", f, err)
		}
	}
	averageDistances(distTxt)

	if len(spJSON.usedIndices()) != 2 || len(spTxt.usedIndices()) != 2 {
		t.Fatalf("expected 2 used species in both modes, got json=%d txt=%d", len(spJSON.usedIndices()), len(spTxt.usedIndices()))
	}
	if math.Abs(distJSON[0][1].Value-2.0) > 1e-9 || math.Abs(distTxt[0][1].Value-2.0) > 1e-9 {
		t.Fatalf("expected averaged distance 2.0, got json=%v txt=%v", distJSON[0][1].Value, distTxt[0][1].Value)
	}
	if math.Abs(qcJSON[0].Average()-12.0) > 1e-9 || math.Abs(qcTxt[0].Average()-12.0) > 1e-9 {
		t.Fatalf("expected Species_1 QC average 12.0, got json=%v txt=%v", qcJSON[0].Average(), qcTxt[0].Average())
	}
	if math.Abs(qcJSON[1].Average()-21.0) > 1e-9 || math.Abs(qcTxt[1].Average()-21.0) > 1e-9 {
		t.Fatalf("expected Species_2 QC average 21.0, got json=%v txt=%v", qcJSON[1].Average(), qcTxt[1].Average())
	}
}

func TestCreateAllMatrixFormats(t *testing.T) {
	dir := t.TempDir()
	stem := filepath.Join(dir, "matrix")

	species := newTestSpecies()
	species.addSampleSpecies("a.mgf", "Species_A")
	species.addSampleSpecies("b.mgf", "Species_B")
	species.addSampleSpecies("c.mgf", "Species_C")

	distances := make(DistanceMatrix)
	qcData := make(map[int]QCData)
	updateDistanceAndQC(species, distances, qcData, "a.mgf", "b.mgf", 0.0, 3.0, 3.0)
	updateDistanceAndQC(species, distances, qcData, "a.mgf", "c.mgf", 2.0, 3.0, 3.0)
	updateDistanceAndQC(species, distances, qcData, "b.mgf", "c.mgf", 2.0, 3.0, 3.0)

	if err := createNexus(stem, distances, species, 0); err != nil {
		t.Fatalf("createNexus failed: %v", err)
	}
	if err := createMega(stem, distances, species, 0.80, qcData); err != nil {
		t.Fatalf("createMega failed: %v", err)
	}
	if err := createMega12(stem+"_MEGA_12_", distances, species, 0.80); err != nil {
		t.Fatalf("createMega12 failed: %v", err)
	}
	if err := createJSON(stem, distances, species, 0.80, qcData); err != nil {
		t.Fatalf("createJSON failed: %v", err)
	}

	nexusBytes, err := os.ReadFile(stem + ".nexus")
	if err != nil {
		t.Fatalf("failed reading nexus: %v", err)
	}
	nexus := string(nexusBytes)
	if !strings.HasPrefix(nexus, "#NEXUS") || !strings.Contains(nexus, "DIMENSIONS ntax=3;") {
		t.Fatalf("unexpected NEXUS header:\n%s", nexus)
	}
	if !strings.Contains(nexus, "0.00000") || !strings.Contains(nexus, "2.00000") {
		t.Fatalf("expected 0.00000 and 2.00000 in NEXUS matrix:\n%s", nexus)
	}

	megaBytes, err := os.ReadFile(stem + ".meg")
	if err != nil {
		t.Fatalf("failed reading mega: %v", err)
	}
	mega := string(megaBytes)
	lfCount := strings.Count(mega, "\n")
	crlfCount := strings.Count(mega, "\r\n")
	if lfCount == 0 || lfCount != crlfCount {
		t.Fatalf("expected CRLF on all lines in MEGA < 12 output (LF=%d, CRLF=%d)", lfCount, crlfCount)
	}
	if !strings.Contains(mega, "#mega") || !strings.Contains(mega, "QC\tSpecies_A\t3.000") {
		t.Fatalf("unexpected MEGA < 12 content:\n%s", mega)
	}

	mega12Bytes, err := os.ReadFile(stem + "_MEGA_12_.meg")
	if err != nil {
		t.Fatalf("failed reading mega12: %v", err)
	}
	if !strings.Contains(string(mega12Bytes), "!Format DataType=Distance DataFormat=LowerLeft NTaxa=3;") {
		t.Fatalf("unexpected MEGA 12 content:\n%s", string(mega12Bytes))
	}

	jsonBytes, err := os.ReadFile(stem + ".json")
	if err != nil {
		t.Fatalf("failed reading json: %v", err)
	}
	var out DistanceMatrixOutput
	if err := json.Unmarshal(jsonBytes, &out); err != nil {
		t.Fatalf("invalid distance matrix JSON: %v", err)
	}
	if len(out.DistanceMatrix.Species) != 3 {
		t.Fatalf("expected 3 species in JSON matrix, got %d", len(out.DistanceMatrix.Species))
	}
	if len(out.DistanceMatrix.Distances) != 2 ||
		out.DistanceMatrix.Distances[0][0] == nil || math.Abs(*out.DistanceMatrix.Distances[0][0]) > 1e-9 ||
		out.DistanceMatrix.Distances[1][0] == nil || math.Abs(*out.DistanceMatrix.Distances[1][0]-2.0) > 1e-9 ||
		out.DistanceMatrix.Distances[1][1] == nil || math.Abs(*out.DistanceMatrix.Distances[1][1]-2.0) > 1e-9 {
		t.Fatalf("unexpected lower-triangular distances: %+v", out.DistanceMatrix.Distances)
	}
}

func TestFixTaxonomyName(t *testing.T) {
	cases := []struct {
		in   string
		want string
	}{
		{in: "Species-1.2", want: "Species-1.2"},
		{in: "*bad name!", want: "_bad_name_"},
		{in: "_lead_ok", want: "_lead_ok"},
	}
	for _, tc := range cases {
		if got := fixTaxonomyName(tc.in); got != tc.want {
			t.Errorf("fixTaxonomyName(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func TestParseUnavailableDistances_JSONAndTxt(t *testing.T) {
	dir := t.TempDir()

	// JSON inputs: one finite (a vs b = 1.5), one "INF" (a vs c), and mixed ("INF" + 2.5 for b vs c)
	j1 := filepath.Join(dir, "ab.json")
	j2 := filepath.Join(dir, "ac.json")
	j3a := filepath.Join(dir, "bc_inf.json")
	j3b := filepath.Join(dir, "bc_fin.json")
	_ = os.WriteFile(j1, []byte(`{"datasetA":"a.mgf","datasetB":"b.mgf","setDistance":1.5,"datasetAQC":10,"datasetBQC":20}`), 0644)
	_ = os.WriteFile(j2, []byte(`{"datasetA":"a.mgf","datasetB":"c.mgf","setDistance":"INF","datasetAQC":10,"datasetBQC":0}`), 0644)
	_ = os.WriteFile(j3a, []byte(`{"datasetA":"b.mgf","datasetB":"c.mgf","setDistance":"INF","datasetAQC":20,"datasetBQC":0}`), 0644)
	_ = os.WriteFile(j3b, []byte(`{"datasetA":"b.mgf","datasetB":"c.mgf","setDistance":2.5,"datasetAQC":20,"datasetBQC":5}`), 0644)

	spJSON := newTestSpecies()
	distJSON := make(DistanceMatrix)
	qcJSON := make(map[int]QCData)
	for i, p := range []string{j1, j2, j3a, j3b} {
		if err := parseComparisonJSON(p, i+1, spJSON, distJSON, qcJSON); err != nil {
			t.Fatalf("parseComparisonJSON(%s) failed: %v", p, err)
		}
	}
	averageDistances(distJSON)

	stem := filepath.Join(dir, "inf_out")
	if err := createNexus(stem, distJSON, spJSON, 2); err != nil {
		t.Fatalf("createNexus failed: %v", err)
	}
	if err := createMega(stem, distJSON, spJSON, 0.80, qcJSON); err != nil {
		t.Fatalf("createMega failed: %v", err)
	}
	if err := createMega12(stem+"_MEGA12", distJSON, spJSON, 0.80); err != nil {
		t.Fatalf("createMega12 failed: %v", err)
	}
	if err := createJSON(stem, distJSON, spJSON, 0.80, qcJSON); err != nil {
		t.Fatalf("createJSON failed: %v", err)
	}

	nex, _ := os.ReadFile(stem + ".nexus")
	if !strings.Contains(string(nex), "?\t") {
		t.Fatalf("expected '?' for unavailable distance in NEXUS:\n%s", string(nex))
	}
	meg, _ := os.ReadFile(stem + ".meg")
	if !strings.Contains(string(meg), "?\t") {
		t.Fatalf("expected '?' for unavailable distance in MEGA:\n%s", string(meg))
	}
	meg12, _ := os.ReadFile(stem + "_MEGA12.meg")
	if !strings.Contains(string(meg12), "        ? ") {
		t.Fatalf("expected '?' for unavailable distance in MEGA12:\n%s", string(meg12))
	}

	jsonBytes, _ := os.ReadFile(stem + ".json")
	var out DistanceMatrixOutput
	if err := json.Unmarshal(jsonBytes, &out); err != nil {
		t.Fatalf("createJSON output invalid: %v\n%s", err, string(jsonBytes))
	}
	// Rows: [0] = [a-b = 1.5], [1] = [a-c = nil, b-c = 2.5]
	if out.DistanceMatrix.Distances[0][0] == nil || math.Abs(*out.DistanceMatrix.Distances[0][0]-1.5) > 1e-9 {
		t.Fatalf("expected [0][0]=1.5, got %v", out.DistanceMatrix.Distances[0][0])
	}
	if out.DistanceMatrix.Distances[1][0] != nil {
		t.Fatalf("expected [1][0]=nil (null) for INF distance, got %v", *out.DistanceMatrix.Distances[1][0])
	}
	if out.DistanceMatrix.Distances[1][1] == nil || math.Abs(*out.DistanceMatrix.Distances[1][1]-2.5) > 1e-9 {
		t.Fatalf("expected [1][1]=2.5 ignoring INF in replicate average, got %v", out.DistanceMatrix.Distances[1][1])
	}

	// TXT input with INF and finite distances
	t1 := filepath.Join(dir, "ab.txt")
	t2 := filepath.Join(dir, "ac.txt")
	_ = os.WriteFile(t1, []byte("dataset_A\ta.mgf\ndataset_B\tb.mgf\nset_distance\t1.25\ndataset_A_QC\t10\ndataset_B_QC\t20\n"), 0644)
	_ = os.WriteFile(t2, []byte("dataset_A\ta.mgf\ndataset_B\tc.mgf\nset_distance\tINF\ndataset_A_QC\t10\ndataset_B_QC\t0\n"), 0644)

	spTxt := newTestSpecies()
	distTxt := make(DistanceMatrix)
	qcTxt := make(map[int]QCData)
	for i, p := range []string{t1, t2} {
		if err := parseComparisonTxt(p, i+1, spTxt, distTxt, qcTxt); err != nil {
			t.Fatalf("parseComparisonTxt(%s) failed: %v", p, err)
		}
	}
	averageDistances(distTxt)
	if v, ok := lookupFiniteDistance(distTxt, 0, 1); !ok || math.Abs(v-1.25) > 1e-9 {
		t.Fatalf("expected TXT finite distance 1.25, got (%v, %v)", v, ok)
	}
	if _, ok := lookupFiniteDistance(distTxt, 0, 2); ok {
		t.Fatalf("expected TXT INF distance to be unavailable")
	}
}
