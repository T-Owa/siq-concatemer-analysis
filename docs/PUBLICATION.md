# Citation and licensing

Repository: [T-Owa/siq-concatemer-analysis](https://github.com/T-Owa/siq-concatemer-analysis). The source archive can be inspected and run without a GitHub account.

## Citing this analysis

Software attribution is provided in [CITATION.cff](../CITATION.cff). Include the repository URL, version and commit used for the analysis. If that version has an archived software DOI, cite the version-specific DOI. No DOI is recorded in the bundled citation metadata.

This repository provides custom post-processing of SIQ event tables. Cite the upstream software separately:

- [SIQ repository](https://github.com/RobinVanSchendel/SIQ)
- [SIQ paper](https://doi.org/10.1093/nargab/lqac063)

The source workbooks report SIQ version 2.1. This repository does not distribute SIQ or perform upstream read merging and SIQ analysis.

## License and data

The [MIT license](../LICENSE) covers the custom scripts and documentation. Experimental data and third-party software are outside that code license. No separate experimental-data license is specified in this repository.

The bundled inputs are processed SIQ event tables, with sample metadata, reference configuration and source hashes. Sequencing FASTQ files are not included. See [Data and outputs](DATA.md) for input fields and provenance.

## Development and validation

The MH algorithm is retained from the earlier analysis. The aggregation and portable execution code were consolidated with AI assistance and checked against archived outputs and independent sequence enumeration. The checks and their scope are documented in [Validation](VALIDATION.md).
