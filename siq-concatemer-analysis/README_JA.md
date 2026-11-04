# SIQコンカテマー接合部解析

SIQ出力から、Event割合、Deletion長・位置、MH候補、主要欠失パターンとTop5を再集計する公開準備用コードです。修飾60サンプル、阻害剤54サンプルの実データを含みます。各条件・時点の独立実験は3回です。

## まず実行する

ZIPを展開し、`README.md`があるフォルダでターミナルを開いて実行します。Ubuntu / WSLでは`python3`を使います。

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

SVG（編集可能）、PDF、300 dpi PNGを出力します。以前のFigure Studioと同じ解析を確認するための静的作図コードであり、Web操作画面や以前の図の見た目を完全に複製するものではありません。

## 今回の解析で区別すること

- 全Eventの分母は、`rawData`に報告された全Typeの`countEvents`合計です。MergedCorrectや入力FASTQの総リード数とは区別します。
- 全サンプルのリードをプールせず、サンプルごとの割合を求めてから独立実験の平均・標本SDを計算します。
- 厳密なMH集計は、Raw全長が欠失座標からの再構成配列と一致し、プライマーを保持して接合部に接する／またがる配置がある純DELETIONを対象にします。
- 6/12 bpのTop5は、その欠失長の純DELETION全体が分母です。内部欠失やRaw不一致のある記録も分母に保持し、QC列で区別します。
- 内部欠失などのMH評価不能はNAです。MH＝0と同じものとして数えません。
- MHは「配列から許容される最大のMH候補」です。MMEJ/TMEJの使用を証明する値ではありません。
- 接合部寄りに表示する座標は、同じ配列を生む等価配置の中から選びます。実際の切断位置を一意に決めたものではありません。

詳しい定義は[Methods](docs/METHODS.md)、元データと出力の説明は[Data](docs/DATA.md)を参照してください。

## GitHubにアップロードする

1. GitHubで新しいリポジトリを作成します。名前の例は`siq-concatemer-analysis`です。
2. ZIPを解凍し、内側の`siq-concatemer-analysis`フォルダを開きます。
3. GitHubの「Add file → Upload files」から、その中のファイルとフォルダをアップロードします。ZIPだけを置くのではなく、READMEやコードが直接見える構造にします。
4. リポジトリの最上位に`README.md`、`scripts/`、`data/`などが並ぶことを確認します。
5. 論文に使う版が確定したらReleaseを作り、Zenodoでその版を保存してDOIを付けます。

GitHubの画面からは1回100ファイル、1ファイル25 MiBまでアップロードできます。詳細は[GitHub公式手順](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository)を参照してください。

`data/`には研究の実データが入っています。コードのみ先に公開する場合でも、論文用の最終版では再現に用いた入力データの公開先を明記してください。

このパッケージを作成した時点では、GitHubへのアップロード・外部公開・DOI取得は行っていません。投稿前に決める項目は[公開用メモ](docs/PUBLICATION.md)にまとめています。
