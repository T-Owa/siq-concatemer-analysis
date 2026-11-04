# Publication and GitHub release notes

This is a prepared repository, not an already published GitHub release. No repository URL or DOI has been invented.

## Before the paper's final release

1. Replace or complete author/contributor information in `CITATION.cff` as appropriate. Add the real GitHub repository URL and, after archiving, the version-specific DOI. A manuscript citation can be added once its title and identifiers are final.
2. The supplied `LICENSE` is an MIT license for the custom code. Confirm that this choice and attribution match the contributors' and institution's intended terms. The MIT code license does not grant a separate license over the experimental data; state the chosen data terms or archive policy when releasing the study data.
3. Link the exact input datasets, configurations and final figure files used in the paper to a release. For selected-pattern figures, report which rows were displayed and retain the original denominator. Removing low-frequency internal patterns from a display does not justify deleting them from the underlying data.
4. Add the raw-read archive accession when available. SIQ processing and upstream read merging are separate from this downstream repository.
5. Make the intended release public and archive it with Zenodo or another DOI-issuing repository. Cite that version in the paper's Code availability statement.

These are publication metadata items; the bundled calculations can already be executed and checked without a GitHub account.

## Code availability template

Replace the bracketed fields only after the repository and archive exist:

> Custom code for post-processing SIQ output, including junction-compatible microhomology annotation, deletion-pattern aggregation and figure generation, is available at [GITHUB_URL] and archived at [VERSION_SPECIFIC_DOI] (version [RELEASE]). The repository includes the reference configuration, sample metadata, SIQ event tables, execution instructions and numerical validation examples. SIQ version 2.1 was used for the upstream sequence analysis.

## Methods wording, Japanese draft

SIQ出力の後処理には独自のPythonスクリプトを用いた。各イベントはcountEventsで重み付けし、各独立実験における割合を算出した後、条件ごとの平均および標本標準偏差を求めた。マイクロホモロジーは、観測された接合配列と完全に一致する単一欠失の等価配置を列挙し、プライマー配列を保持してドナー接合部に接する、またはまたがる配置に限定して、許容される最大の相同性長として評価した。評価不能なイベントはMH＝0と区別した。欠失パターンの集計にはSIQの欠失座標から再構成したアレルを用い、観測Rawとの不一致はQC情報として保持した。各解析の対象イベント、分母およびコードは公開リポジトリに記載した。

## Related software and instructions

- SIQ: https://github.com/RobinVanSchendel/SIQ
- SIQ paper: https://doi.org/10.1093/nargab/lqac063
- GitHub file upload: https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository
- GitHub/Zenodo archiving: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content

The original MH algorithm was recovered from the earlier analysis archive. The aggregation and portable execution code were consolidated with AI assistance and checked against archived outputs and independent sequence enumeration. This records how the package was prepared; responsibility for scientific interpretation remains with the study authors.
