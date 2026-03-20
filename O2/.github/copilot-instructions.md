# ALICE O2 Cascade Correlation Analysis

## Overview
This workspace contains analysis code for studying cascade (Ξ and Ω) particle correlations in ALICE Run 3 data. The workflow involves running O2Physics analysis tasks on AO2D data files stored on the ALICE Grid, followed by offline ROOT-based postprocessing.

## Core Architecture

### Three-Component System
1. **O2Physics Task** (`/home/rik/alice/O2Physics/PWGLF/Tasks/Strangeness/cascadecorrelations.cxx`)
   - Upstream C++ analysis task in the O2Physics framework
   - Processes AO2D files through DPL (Data Processing Layer) workflows
   - Produces `AnalysisResults.root` containing THnSparse histograms
   - Modified/maintained separately from this workspace

2. **Steering Scripts** (`runtest.sh`, `runmctest.sh`, `runclosuretest.sh`)
   - Execute O2 workflows from within the O2Physics environment
   - Chain multiple analysis services (event selection, PID, propagation, etc.)
   - Access ALICE Grid data via `alien://` protocol

3. **Postprocessing Scripts** (`.cxx` ROOT macros in workspace root)
   - `postprocessingResults.cxx`: Extract correlation functions and invariant mass spectra
   - `MCEfficiency.cxx`: Calculate reconstruction efficiencies from MC
   - `MCClosure.cxx`: Validate analysis with MC closure tests
   - `systematics.cxx`: Evaluate systematic uncertainties

### Data Flow
```
AO2D (Grid) → O2 Workflow → AnalysisResults.root → Postprocessing Macros → plots.root / PDF outputs
```

## Essential Workflows

### Postprocessing with ROOT (Primary Workflow)
The main work happens in postprocessing macros that analyze `AnalysisResults.root` files:

**Core postprocessing scripts:**
```bash
root 'postprocessingResults.cxx("trainnr")' -q
root 'MCEfficiency.cxx("trainnr", false, "correlations")' -q
root 'MCClosure.cxx("trainnr")' -q
root 'systematics.cxx("trainnr")' -q
```

**File structure expectations:**
- Input: `results/TRAINNR/AnalysisResults.root` (or use "test" for `./AnalysisResults.root`)
- Output: `plots/` directory (ROOT files) and `figures/` directory (PDFs for thesis)
- Efficiency maps: `CascadeEfficiencies*.root` files

**Common workflows:**
1. Extract correlations and mass spectra: `postprocessingResults.cxx`
2. Generate efficiency maps from MC: `MCEfficiency.cxx` (maps manually uploaded for Grid production)
3. Validate analysis with MC closure: `MCClosure.cxx`
4. Evaluate systematic uncertainties: `systematics.cxx`

### Running O2 Analysis (Testing Only)
O2 workflows are run **locally for testing compilation and runtime only**. Production runs happen on the Grid (user responsibility).

To test the O2 task locally:
```bash
alienv setenv O2Physics/latest -c bash runtest.sh
```

**Important limitations:**
- Grid token initialization (`alien-token-init`) requires user credentials - cannot be automated
- Grid issues (CCDB access, alien:// paths) are user responsibility
- Configuration files: `cascadeconfig.json`, `ssbarconfig.json` (cuts and binning)
- Check `log_o2.txt` for O2 workflow output/errors

## Project-Specific Patterns

### THnSparse Analysis Convention
All correlation analysis uses ROOT THnSparse objects with standard axis ordering defined in code enums:
```cpp
// postprocessingResults.cxx
struct corr {
  enum { dPhi, dY, signTrigg, signAssoc, ptTrigg, ptAssoc, 
         invMassTrigg, invMassAssoc, V_z, multiplicity };
};
```
Use `postprocessingTools.h::project()` function for THnSparse projections with axis range cuts.

### Grid Data Access
- AO2D paths: `alien:///alice/data/YEAR/LHCPERIOD/RUNNR/PASS/...`
- MC paths: `alien:///alice/sim/YEAR/LHCPERIOD/...`
- File lists: `Run3Samplefiles.txt`, `run3MCfiles.txt`
- Scripts auto-detect environment (local vs. Nikhef STBC cluster)

### Signal Extraction (postprocessingResults.cxx)
Two-step process for background subtraction:
1. **Fit invariant mass**: Gaussian + pol2 background to Ξ/Ω mass peaks
   - Background fit functions: `pol2bkgXi()`, `pol2bkgom()` in `postprocessingTools.h`
   - Excludes signal region during fit (hardcoded mass windows)
2. **Define regions**: Signal/sideband windows in σ units
   - `signalWindow = 3.0` (±3σ around peak)
   - `sidebandWindow = {4.0, 10.0}` (4-10σ)
   - Arrays `sigXi[][]`, `SBlowXi[][]`, `sigOm[][]` store boundaries per pT bin
3. **Extract correlations**: Project THnSparse separately for signal/sideband regions
4. **Subtract background**: Signal = (Signal region) - (Sideband region)

## Key Files & Directories

- `postprocessingTools.h`: Shared utility functions (THnSparse projection, background fit functions)
- `cascadeconfig.json`: **Critical** - defines all analysis cuts and binning
- `/home/rik/alice/O2Physics/PWGLF/`: Upstream O2Physics repository (read-only workflow dependency)
- `results/`: Expected location for train outputs
- `figures/`: Output directory for plots (referenced in thesis)

## Critical Conventions

1. **Never run O2 workflows outside `alienv setenv O2Physics/latest`** - dependencies won't resolve
2. **pT binning must match** between O2 task configuration and postprocessing macros
3. **Axis indices in enums** are hardcoded - changing THnSparse structure breaks all macros
4. **ROOT macro execution**: Use `root 'macro.cxx(args)' -q` syntax, not compiled binaries
5. **ALICE Grid access**: Requires valid GRID certificate and `alien-token-init`

## Testing & Validation

- **MC Closure**: Compare reconstructed to generated distributions (should be flat ratio)
- **Efficiency QA**: Check for pt/eta dependencies in efficiency corrections
- **Systematic variations**: Run with different cut sets (Default/Loose/Tight configs)

## Common Issues

**Postprocessing errors:**
- **"File not found"**: Check `results/TRAINNR/` directory structure or use "test" for local file
- **Wrong signal/background regions**: Recalculate boundaries after changing mass cuts or pT bins
- **Empty histograms**: Verify pT binning matches between O2 config and macro arrays
- **Segmentation fault**: Usually axis index mismatch - check THnSparse enum vs actual axis order

**O2 workflow issues (testing only):**
- **"Cannot find CCDB object"**: Grid token issue - user must run `alien-token-init` manually
- **O2 crashes**: Check `log_o2.txt`, often memory limits (`--shm-segment-size`)
- **Grid access failures**: User responsibility - requires valid GRID certificate and credentials
