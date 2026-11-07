# SIQコンカテマー接合部解析

SIQ出力から、Event割合、Deletion長・位置、MH候補、主要欠失パターンとTop5を再集計するPythonコードです。修飾60サンプル、阻害剤54サンプルの実データを含みます。各条件・時点の独立実験は3回です。

編集者・査読者向けの実行手順と解析・出力の対応は、英語の[Reviewer guide](docs/REVIEWER_GUIDE.md)にまとめています。

## まず実行する

リポジトリを取得し、`README.md`があるフォルダでターミナルを開いて実行します。Ubuntu / WSLでは`python3`を使います。Windows PowerShellでは、環境に応じて`python3`を`python`に置き換えてください。

```bash
python3 scripts/run_analysis.py --outdir results
python3 scripts/verify_reference_results.py --results results
python3 -m unittest discover -s tests -v
```

同梱CSVからの集計だけならPython 3.10以上の標準機能で動きます。`results`にファイルが既にある場合は、`results_2`など別の名前を指定してください。元データや既存結果は上書きしません。

図も出す場合は、Python 3.11以上を使って以下を実行します（検証環境はPython 3.12.14です）。集計のみはPython 3.10以上で動きます。

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 scripts/plot_figures.py --results results --outdir results/figures
```

集計表からSVG（編集可能）、PDF、300 dpi PNGを出力します。作図によって解析値は変わりません。

## 解析の定義

解析ごとの対象イベントと分母、MH・欠失パターンの定義は[Methods](docs/METHODS.md)、入力と出力の説明は[Data](docs/DATA.md)を参照してください。

---

[コードのライセンス](LICENSE) · [引用情報](CITATION.cff) · [利用・引用の詳細](docs/PUBLICATION.md)
