package main

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

type testSpectrumSpec struct {
	scan      int
	rt        float64
	pepMass   float64
	startMz   float64
	baseInten float64
}

func writeTestMGFFile(t *testing.T, dir, name string, specs []testSpectrumSpec, withBlankLines bool) string {
	t.Helper()
	var b strings.Builder
	for idx, s := range specs {
		if withBlankLines && idx == 0 {
			b.WriteString("\n")
		}
		b.WriteString("BEGIN IONS\n")
		b.WriteString("TITLE=spec\n")
		if withBlankLines {
			b.WriteString("\n")
		}
		b.WriteString("SCANS=")
		b.WriteString(itoa(s.scan))
		b.WriteString("\n")
		b.WriteString("RTINSECONDS=")
		b.WriteString(ftoa(s.rt))
		b.WriteString("\n")
		b.WriteString("PEPMASS=")
		b.WriteString(ftoa(s.pepMass))
		b.WriteString("\n")
		b.WriteString("CHARGE=2+\n")
		for i := 0; i < 25; i++ {
			mz := s.startMz + float64(i)*10.0
			inten := s.baseInten + float64(i)*50.0
			b.WriteString(ftoa(mz))
			b.WriteString(" ")
			b.WriteString(ftoa(inten))
			b.WriteString("\n")
		}
		b.WriteString("END IONS\n")
		if withBlankLines {
			b.WriteString("\n")
		}
	}
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(b.String()), 0644); err != nil {
		t.Fatalf("failed writing %s: %v", p, err)
	}
	return p
}

func itoa(v int) string {
	if v == 0 {
		return "0"
	}
	var buf [20]byte
	i := len(buf)
	for v > 0 {
		i--
		buf[i] = byte('0' + v%10)
		v /= 10
	}
	return string(buf[i:])
}

func ftoa(v float64) string {
	// Compact integer-like formatting sufficient for synthetic test peaks/metadata
	iv := int(v)
	if float64(iv) == v {
		return itoa(iv) + ".0"
	}
	whole := int(v)
	frac := int((v-float64(whole))*100 + 0.5)
	if frac < 10 {
		return itoa(whole) + ".0" + itoa(frac)
	}
	return itoa(whole) + "." + itoa(frac)
}

func readMGFWithParams(t *testing.T, mgfPath string, topN int) (*Parameters, *Dataset) {
	t.Helper()
	par := NewParameters()
	par.TwoPass = true
	par.TopN = topN
	ds := &Dataset{ID: 'A', Filename: mgfPath}
	if err := ReadMGF(par, ds); err != nil {
		t.Fatalf("ReadMGF(TopN=%d) failed: %v", topN, err)
	}
	return par, ds
}

func TestReadMGF_TwoPass_AllTiedIntensities(t *testing.T) {
	dir := t.TempDir()
	specs := []testSpectrumSpec{
		{scan: 1, rt: 10, pepMass: 500.25, startMz: 150, baseInten: 1000},
		{scan: 2, rt: 20, pepMass: 600.25, startMz: 450, baseInten: 1000},
		{scan: 3, rt: 30, pepMass: 700.25, startMz: 750, baseInten: 1000},
	}
	mgfPath := writeTestMGFFile(t, dir, "tied.mgf", specs, false)

	_, ds := readMGFWithParams(t, mgfPath, 2)
	if len(ds.Spectra) != 2 {
		t.Fatalf("expected 2 spectra for TopN=2 with all tied intensities, got %d", len(ds.Spectra))
	}
	if ds.Spectra[0].Scan != 1 || ds.Spectra[1].Scan != 2 {
		t.Fatalf("expected scans [1, 2], got [%d, %d]", ds.Spectra[0].Scan, ds.Spectra[1].Scan)
	}
}

func TestReadMGF_TwoPass_BoundaryTiePreservesHigher(t *testing.T) {
	dir := t.TempDir()
	specs := []testSpectrumSpec{
		{scan: 1, rt: 10, pepMass: 500.25, startMz: 150, baseInten: 1000},
		{scan: 2, rt: 20, pepMass: 600.25, startMz: 450, baseInten: 2000},
		{scan: 3, rt: 30, pepMass: 700.25, startMz: 750, baseInten: 2000},
		{scan: 4, rt: 40, pepMass: 800.25, startMz: 1050, baseInten: 3000},
	}
	mgfPath := writeTestMGFFile(t, dir, "boundary_tie.mgf", specs, false)

	_, ds := readMGFWithParams(t, mgfPath, 2)
	if len(ds.Spectra) != 2 {
		t.Fatalf("expected 2 spectra for TopN=2, got %d", len(ds.Spectra))
	}
	if ds.Spectra[0].Scan != 2 || ds.Spectra[1].Scan != 4 {
		t.Fatalf("expected scans [2, 4] preserving file order and highest spectrum, got [%d, %d]",
			ds.Spectra[0].Scan, ds.Spectra[1].Scan)
	}
}

func TestReadMGF_TwoPass_EmptyLines(t *testing.T) {
	dir := t.TempDir()
	specs := []testSpectrumSpec{
		{scan: 1, rt: 10, pepMass: 500.25, startMz: 150, baseInten: 1000},
		{scan: 2, rt: 20, pepMass: 600.25, startMz: 450, baseInten: 2000},
		{scan: 3, rt: 30, pepMass: 700.25, startMz: 750, baseInten: 3000},
	}
	mgfPath := writeTestMGFFile(t, dir, "blanks.mgf", specs, true)

	_, ds := readMGFWithParams(t, mgfPath, 2)
	if len(ds.Spectra) != 2 {
		t.Fatalf("expected 2 spectra with blank lines and TopN=2, got %d", len(ds.Spectra))
	}
	if ds.Spectra[0].Scan != 2 || ds.Spectra[1].Scan != 3 {
		t.Fatalf("expected scans [2, 3], got [%d, %d]", ds.Spectra[0].Scan, ds.Spectra[1].Scan)
	}
}

func TestReadMGF_TwoPass_TopNDisabledAndZero(t *testing.T) {
	dir := t.TempDir()
	specs := []testSpectrumSpec{
		{scan: 1, rt: 10, pepMass: 500.25, startMz: 150, baseInten: 1000},
		{scan: 2, rt: 20, pepMass: 600.25, startMz: 450, baseInten: 2000},
		{scan: 3, rt: 30, pepMass: 700.25, startMz: 750, baseInten: 3000},
	}
	mgfPath := writeTestMGFFile(t, dir, "modes.mgf", specs, false)

	_, dsAll := readMGFWithParams(t, mgfPath, -1)
	if dsAll.TwoPassFilterActive || len(dsAll.Spectra) != 3 {
		t.Fatalf("TopN=-1 expected 3 spectra and TwoPassFilterActive=false, got len=%d active=%v",
			len(dsAll.Spectra), dsAll.TwoPassFilterActive)
	}

	_, dsZero := readMGFWithParams(t, mgfPath, 0)
	if !dsZero.TwoPassFilterActive || len(dsZero.Spectra) != 0 {
		t.Fatalf("TopN=0 expected 0 spectra and TwoPassFilterActive=true, got len=%d active=%v",
			len(dsZero.Spectra), dsZero.TwoPassFilterActive)
	}

	_, dsGT := readMGFWithParams(t, mgfPath, 10)
	if dsGT.TwoPassFilterActive || len(dsGT.Spectra) != 3 {
		t.Fatalf("TopN=10 expected 3 spectra and TwoPassFilterActive=false, got len=%d active=%v",
			len(dsGT.Spectra), dsGT.TwoPassFilterActive)
	}
}

func TestWriteOutputs_EmptyDatasetNoNaN(t *testing.T) {
	dir := t.TempDir()
	hist := make([]int, dotprodHistogramBins)

	cases := []struct {
		name string
		aLen int
		bLen int
	}{
		{name: "both_empty", aLen: 0, bLen: 0},
		{name: "a_empty", aLen: 0, bLen: 1},
		{name: "b_empty", aLen: 1, bLen: 0},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			par := NewParameters()
			par.Metric = 2
			par.OutputFilename = filepath.Join(dir, tc.name+".txt")
			par.JSONFilename = filepath.Join(dir, tc.name+".json")

			dsA := &Dataset{ID: 'A', Filename: "a.mgf", Spectra: make([]Spectrum, tc.aLen)}
			dsB := &Dataset{ID: 'B', Filename: "b.mgf", Spectra: make([]Spectrum, tc.bLen)}

			if err := WriteJSON(par, dsA, dsB, hist, nil, 0, 0, 0, 0, tc.aLen, tc.bLen); err != nil {
				t.Fatalf("WriteJSON failed: %v", err)
			}
			if err := WriteTextOutput(par, dsA, dsB, hist, hist, 0, 0, 0, 0); err != nil {
				t.Fatalf("WriteTextOutput failed: %v", err)
			}

			rawJSON, err := os.ReadFile(par.JSONFilename)
			if err != nil {
				t.Fatalf("failed reading JSON: %v", err)
			}
			var parsed map[string]any
			if err := json.Unmarshal(rawJSON, &parsed); err != nil {
				t.Fatalf("invalid JSON: %v", err)
			}
			if parsed["setDistance"] != "INF" {
				t.Fatalf("expected setDistance == INF in JSON, got %v", parsed["setDistance"])
			}

			rawTxt, err := os.ReadFile(par.OutputFilename)
			if err != nil {
				t.Fatalf("failed reading TXT: %v", err)
			}
			if !strings.Contains(string(rawTxt), "set_distance\tINF\n") {
				t.Fatalf("expected set_distance\\tINF in text output, got:\n%s", string(rawTxt))
			}
		})
	}
}

func TestReadMGF_ExcludesOutOfRangeSpectra(t *testing.T) {
	dir := t.TempDir()
	specs := []testSpectrumSpec{
		{scan: 1, rt: 10, pepMass: 500.25, startMz: 150, baseInten: 1000},
		{scan: 5, rt: 50, pepMass: 600.25, startMz: 450, baseInten: 5000},
		{scan: 10, rt: 100, pepMass: 700.25, startMz: 750, baseInten: 3000},
		{scan: 15, rt: 150, pepMass: 800.25, startMz: 1050, baseInten: 2000},
	}
	mgfPath := writeTestMGFFile(t, dir, "range.mgf", specs, false)

	// Scan range [4, 12] + RT range [40, 120] keeps only scans 5 and 10, while SpectrumCount remains 4
	for _, twoPass := range []bool{false, true} {
		par := NewParameters()
		par.TwoPass = twoPass
		par.TopN = 1
		par.StartScan = 4
		par.EndScan = 12
		par.StartRT = 40
		par.EndRT = 120
		ds := &Dataset{ID: 'A', Filename: mgfPath}
		if err := ReadMGF(par, ds); err != nil {
			t.Fatalf("ReadMGF(twoPass=%v) failed: %v", twoPass, err)
		}
		if ds.SpectrumCount != 4 {
			t.Fatalf("twoPass=%v expected total SpectrumCount=4, got %d", twoPass, ds.SpectrumCount)
		}
		if len(ds.Spectra) != 1 || ds.Spectra[0].Scan != 5 {
			t.Fatalf("twoPass=%v expected only highest in-range scan [5], got %+v", twoPass, ds.Spectra)
		}
	}

	// Range excluding all spectra leaves len(Spectra)==0 so metric 2 emits INF
	par := NewParameters()
	par.Metric = 2
	par.StartScan = 100
	par.EndScan = 200
	par.JSONFilename = filepath.Join(dir, "excluded.json")
	ds := &Dataset{ID: 'A', Filename: mgfPath}
	if err := ReadMGF(par, ds); err != nil {
		t.Fatalf("ReadMGF(all excluded) failed: %v", err)
	}
	if ds.SpectrumCount != 4 || len(ds.Spectra) != 0 {
		t.Fatalf("expected SpectrumCount=4 and 0 stored spectra, got %d / %d", ds.SpectrumCount, len(ds.Spectra))
	}
	hist := make([]int, dotprodHistogramBins)
	if err := WriteJSON(par, ds, ds, hist, nil, 0, 0, 0, 0, 0, 0); err != nil {
		t.Fatalf("WriteJSON failed: %v", err)
	}
	var parsed map[string]any
	raw, _ := os.ReadFile(par.JSONFilename)
	if err := json.Unmarshal(raw, &parsed); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	if parsed["setDistance"] != "INF" {
		t.Fatalf("expected setDistance == INF when all spectra excluded, got %v", parsed["setDistance"])
	}
}

func TestWriteJSON_EligibleComparedCounts(t *testing.T) {
	dir := t.TempDir()
	// 4 spectra in file:
	// scan 1: out of scan range (-W 2,10)
	// scan 2: too few peaks (only 2 peaks, -I 20)
	// scan 3: below basepeak threshold (500, -m 2000,0)
	// scan 4: eligible (25 peaks, basepeak 5000, scan 4)
	specsInRange := []testSpectrumSpec{
		{scan: 1, rt: 10, pepMass: 500.25, startMz: 150, baseInten: 5000},
		{scan: 3, rt: 30, pepMass: 600.25, startMz: 450, baseInten: 500},
		{scan: 4, rt: 40, pepMass: 700.25, startMz: 750, baseInten: 5000},
	}
	mgfA := writeTestMGFFile(t, dir, "eligA.mgf", specsInRange, false)
	// Append one low-peak spectrum (scan 2, 2 peaks) to eligA.mgf
	f, err := os.OpenFile(mgfA, os.O_APPEND|os.O_WRONLY, 0644)
	if err != nil {
		t.Fatalf("failed opening %s: %v", mgfA, err)
	}
	_, _ = f.WriteString("BEGIN IONS\nSCANS=2\nRTINSECONDS=20.0\nPEPMASS=550.25\n300.0 5000.0\n310.0 5000.0\nEND IONS\n")
	f.Close()

	// File B has 3 spectra: 1 out of range (scan 20), 2 eligible (scans 4, 5)
	mgfB := writeTestMGFFile(t, dir, "eligB.mgf", []testSpectrumSpec{
		{scan: 4, rt: 40, pepMass: 700.25, startMz: 750, baseInten: 5000},
		{scan: 5, rt: 50, pepMass: 800.25, startMz: 1050, baseInten: 5000},
		{scan: 20, rt: 200, pepMass: 900.25, startMz: 1350, baseInten: 5000},
	}, false)

	par := NewParameters()
	par.StartScan = 2
	par.EndScan = 10
	par.MinPeaks = 20
	par.MinBasepeakIntensity = 2000
	par.JSONFilename = filepath.Join(dir, "elig.json")

	dsA := &Dataset{ID: 'A', Filename: mgfA}
	dsB := &Dataset{ID: 'B', Filename: mgfB}
	if err := ReadMGF(par, dsA); err != nil {
		t.Fatalf("ReadMGF A failed: %v", err)
	}
	if err := ReadMGF(par, dsB); err != nil {
		t.Fatalf("ReadMGF B failed: %v", err)
	}

	hist := make([]int, dotprodHistogramBins)
	var nComp, gt, sAB, sBA, cmpA, cmpB int
	ComputeDotProdHistogram(par, dsA, dsB, hist, nil, &nComp, &gt, &sAB, &cmpA, true)
	ComputeDotProdHistogram(par, dsB, dsA, hist, nil, &nComp, &gt, &sBA, &cmpB, false)

	if err := WriteJSON(par, dsA, dsB, hist, nil, nComp, gt, sAB, sBA, cmpA, cmpB); err != nil {
		t.Fatalf("WriteJSON failed: %v", err)
	}

	raw, err := os.ReadFile(par.JSONFilename)
	if err != nil {
		t.Fatalf("read JSON failed: %v", err)
	}
	var out JSONOutput
	if err := json.Unmarshal(raw, &out); err != nil {
		t.Fatalf("unmarshal failed: %v", err)
	}
	if out.DatasetAQC != 4 || out.DatasetBQC != 3 {
		t.Fatalf("expected raw SpectrumCount QC 4/3, got %d/%d", out.DatasetAQC, out.DatasetBQC)
	}
	if out.DatasetACompared != 1 || out.DatasetBCompared != 2 {
		t.Fatalf("expected filter-eligible compared counts 1/2, got %d/%d", out.DatasetACompared, out.DatasetBCompared)
	}
}
