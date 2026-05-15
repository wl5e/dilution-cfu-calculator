# Dilution Plate CFU Calculator

A GMP-compliant command-line tool for calculating colony-forming units (CFU) from serial dilution plate count data with uncertainty propagation.

## Overview

This tool implements USP/EP harmonized methodology for microbial enumeration assays commonly performed in pharmaceutical quality control laboratories. It handles spread plate and pour plate methods, flags TNTC (Too Numerous To Count) and TFTC (Too Few To Count) results, and generates formatted reports suitable for GMP documentation.

## Features

- **CFU/mL & CFU/g calculation** — Supports both liquid and solid sample enumeration
- **Countable range filtering** — Defaults: 30–300 CFU/plate (bacteria) or 15–150 (yeast/mold), per USP <61>/<62>
- **Uncertainty propagation** — Reports coefficient of variation (CV%) using multi-plate statistics or Poisson estimation
- **TNTC / TFTC handling** — Automatic exclusion with fallback when no plates are in countable range
- **Multi-dilution weighted averaging** — Combines data from multiple dilution levels
- **CSV & JSON input** — Flexible data entry formats
- **Quick CLI mode** — Fast calculations from command-line arguments
- **Machine-readable JSON output** — For integration with LIMS or data pipelines
- **Example template generator** — Creates ready-to-fill CSV files

## Installation

```bash
git clone https://github.com/collins-amatu-gorgerat/dilution-cfu-calculator.git
cd dilution-cfu-calculator
pip install -r requirements.txt
```

No external dependencies required for core functionality — only `pytest` for running tests.

## Usage

### Quick Calculation

```bash
# 2 dilutions × 2 replicates = 4 plates
python dilution_cfu.py quick \
  --counts "245,231,42,38" \
  --dilutions "0.01,0.001" \
  --replicates 2
```

### Calculate from CSV File

```bash
python dilution_cfu.py calculate plate_data.csv
```

CSV format (columns are auto-detected by header):

| dilution_factor | volume_plated_ml | colony_count | replicate_id |
|-----------------|------------------|--------------|---------------|
| 0.1             | 0.1              | TNTC         | A             |
| 0.01            | 0.1              | 245          | A             |
| 0.01            | 0.1              | 231          | B             |
| 0.001           | 0.1              | 42           | A             |
| 0.001           | 0.1              | 38           | B             |

### Calculate from JSON File

```bash
python dilution_cfu.py calculate plate_data.json --json-output
```

### Generate Example Template

```bash
python dilution_cfu.py example my_experiment.csv
```

### With Sample Mass (CFU/g)

```bash
python dilution_cfu.py calculate data.csv --sample-mass-g 1.0 --diluent-ml 9.0
```

### Yeast / Mold Enumeration

```bash
python dilution_cfu.py calculate data.csv --organism yeast_mold
```

## Output Example

```
==============================================================
        DILUTION PLATE CFU ANALYSIS REPORT
==============================================================
  Method:                  Spread Plate
  Countable Range:         30–300 CFU/plate
  Plates in Analysis:      4
  Dilutions Used:          1.0e-02, 1.0e-03
--------------------------------------------------------------
  CFU/mL:                  2.38e+05
  Uncertainty (CV%):       4.12%
--------------------------------------------------------------
==============================================================
```

## Requirements

- Python 3.9+
- No external runtime dependencies (standard library only)

## Testing

```bash
pytest tests/ -v
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Author

Collins Amatu Gorgerat — Microbiologist, Pharmaceutical GMP / Bioinformatics
