# Two-Pass & Empty-Dataset Edge Case MGF Reproducers

These minimal synthetic MGF files reproduce the edge cases fixed in commit `f52bb46` (`compareMS2/compareMS2.go`). Each spectrum has 25 peaks within `150–1290` *m/z*, matching `compareMS2`'s default filters (`-I 20 -L 109 -U 2000`).

## Files & What They Reproduce

| File | Scenario | Before Fix (`8a4f9cf`) | After Fix (`f52bb46`) |
| :--- | :--- | :--- | :--- |
| `tied_intensities.mgf` | 3 spectra with identical summed intensity (`40000.0`), run with `-two-pass -N 2` | Keeps `0` spectra (`40000 > 40000` is false); `WriteJSON` fails with `json: unsupported value: NaN` | Keeps `2` spectra; `setDistance = 0`, `nrGtCutoff = 4` |
| `boundary_tie.mgf` | 4 spectra with summed intensities `[40000, 65000, 65000, 90000]`, run with `-two-pass -N 2` | Keeps only `1` spectrum (scan 4), dropping both boundary-tied `65000` spectra | Keeps `2` spectra (scans 2 and 4); `setDistance = 0`, `nrGtCutoff = 4` |
| `blank_lines.mgf` | 3 spectra (`[40000, 65000, 90000]`) with blank lines before/inside/between `BEGIN IONS` blocks, run with `-two-pass -N 2` | Panics in `PreCheckMGF` (`index out of range [0] with length 0`) | Keeps `2` spectra (scans 2 and 3); `setDistance = 0`, `nrGtCutoff = 4` |
| `tied_intensities.mgf` (`-N 0`) | 0 spectra selected (`-two-pass -N 0`) under metric 2 (`-d 2`) | `WriteJSON` fails with `json: unsupported value: NaN` (`0/0`) | Exits `0`; outputs `"setDistance": "INF"` |

## How to Run (from repository root `compareMS2-cmd`)

**Windows (PowerShell):**

```powershell
.\compareMS2.exe -A examples\two_pass_edge_case\tied_intensities.mgf -B examples\two_pass_edge_case\tied_intensities.mgf -two-pass -N 2 -o $env:TEMP\tied_n2.txt -J $env:TEMP\tied_n2.json
.\compareMS2.exe -A examples\two_pass_edge_case\boundary_tie.mgf -B examples\two_pass_edge_case\boundary_tie.mgf -two-pass -N 2 -o $env:TEMP\boundary_n2.txt -J $env:TEMP\boundary_n2.json
.\compareMS2.exe -A examples\two_pass_edge_case\blank_lines.mgf -B examples\two_pass_edge_case\blank_lines.mgf -two-pass -N 2 -o $env:TEMP\blank_n2.txt -J $env:TEMP\blank_n2.json
.\compareMS2.exe -A examples\two_pass_edge_case\tied_intensities.mgf -B examples\two_pass_edge_case\tied_intensities.mgf -two-pass -N 0 -o $env:TEMP\tied_n0.txt -J $env:TEMP\tied_n0.json
```

**Linux / macOS (Bash):**

```bash
./compareMS2 -A examples/two_pass_edge_case/tied_intensities.mgf -B examples/two_pass_edge_case/tied_intensities.mgf -two-pass -N 2 -o /tmp/tied_n2.txt -J /tmp/tied_n2.json
./compareMS2 -A examples/two_pass_edge_case/boundary_tie.mgf -B examples/two_pass_edge_case/boundary_tie.mgf -two-pass -N 2 -o /tmp/boundary_n2.txt -J /tmp/boundary_n2.json
./compareMS2 -A examples/two_pass_edge_case/blank_lines.mgf -B examples/two_pass_edge_case/blank_lines.mgf -two-pass -N 2 -o /tmp/blank_n2.txt -J /tmp/blank_n2.json
./compareMS2 -A examples/two_pass_edge_case/tied_intensities.mgf -B examples/two_pass_edge_case/tied_intensities.mgf -two-pass -N 0 -o /tmp/tied_n0.txt -J /tmp/tied_n0.json
```
