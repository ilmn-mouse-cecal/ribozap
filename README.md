# ribozap

**Create Additional Probes on Unknown Mixed Samples**

## Overview

**ribozap** is designed to generate custom RNaseH probes for metagenomic RNA samples without prior knowledge of the genomes or species present. This boosts the performance of the RiboZero Plus kit by targeting additional ribosomal RNA (rRNA) regions for depletion.

Depleting rRNA, which typically constitutes 80–90% of total RNA, before sequencing significantly reduces sequencing costs and increases the coverage of biologically relevant transcripts.


## Installation

**Prerequisites**

- Docker
- Python3

**Clone the repository**:

```bash
git clone https://github.com/yourusername/ribozap.git
cd ribozap
pip install -e .
docker build -t ribozap:latest .

ribozap \
  --sample-sheet samplesheet.csv \
  --output-dir ./my_out/ \
  --analysis-name my_analysis \
  --cpus 4 \
  --memory 16 \
  --resume
```


## Code Style

The codebase follows [PEP 8](https://peps.python.org/pep-0008/) conventions and uses [Black](https://github.com/psf/black) for formatting.

## Citation

If you use `ribozap` in your research, please cite our publication:

Roos M, Bunga S, Tan A, Maissy E, Skola D, Richter A, Whittaker DS, Desplats P, Zarrinpar A, Conrad R, Kuersten S. (2025). *Optimizing mouse metatranscriptome profiling by selective removal of redundant nucleic acid sequences*. **bioRxiv**. https://pubmed.ncbi.nlm.nih.gov/39868335/