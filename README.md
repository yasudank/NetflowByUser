# PFS Netflow Pipeline for Classical Observation - run_netflow.py

`run_netflow.py` は、PFS (Prime Focus Spectrograph) のファイバー割り当て最適化計算 (netflow) と、それに続く `pfsDesign` FITS ファイルおよび OPE 観測用制御ファイルの生成・検証を自動で実行するエンドツーエンドのパイプラインスクリプトです。以下の天文台で作成されたプログラムを参考に作成しました。コードのほぼすべては Google Gemini で生成されました。
* pfs_obsproc_planning_tools: https://github.com/Subaru-PFS/pfs_obsproc_planning_tools
* ets_pointing: https://github.com/Subaru-PFS/ets_pointing
* ets_fiberalloc: https://github.com/Subaru-PFS/ets_fiberalloc

---

## 0. 準備: 

### 0.1. Gurobi Optimizer のライセンス準備

ファイバー配置最適化（Netflow）の計算を実行するためには、商用数理計画ソルバーである **Gurobi Optimizer** が必要です。

* **ライセンスについて**: Gurobi Optimizer の利用にはライセンスキーの取得とアクティベーションが必要です。
* **アカデミック無償ライセンス**: 大学などの学術機関に所属している場合は、以下のページから**無償のアカデミック版ライセンス**を取得できます：
  * [Gurobi アカデミックライセンス](https://www.gurobi.com/jp/license#academic-license)
* **アクティベーション**: ライセンス取得後、配布される `grbgetkey` コマンドを実行して環境にライセンスを登録（アクティベーション）してください。

---

### 0.2. 仮想環境の構築

PFS 関連パッケージや solver などの依存関係が整理された Python 仮想環境は、同一ディレクトリの `create_pfs_env.py` を用いて自動生成・管理することができます。

```bash
# 設定ファイル (例: config.yaml) のパッケージ依存定義を解決して仮想環境を作成
python3 create_pfs_env.py path/to/config.yaml --venv-dir .venv

# 作成した仮想環境を有効化
source .venv/bin/activate
```

### create_pfs_env.py の主要な引数とオプション
* `config_file` (必須): `pfs:` および `gurobi:` セクションを含む YAML ファイル（例: `config.yaml`）へのパス。通常この config.yaml は観測ランごとに 観測所から SSP チームに提供されます。
* `-o`, `--venv-dir` : 仮想環境を作成する場所を指定します（デフォルト: `.venv`）。
* `-p`, `--python` : 作成する Python のバージョンを指定します（デフォルト: `3.12`）。
* `--use-local` : ディレクトリ内に `ets_fiberalloc` や `ics_cobraOps` などのソースコードがある場合、それらをローカル編集モード (`pip install -e`) でインストールします。

---

### 0.3. Calibration catalog (sky, fluxstd) と Guide Start catalog (gaia) の準備

pfsa サーバーを経由して、targetdbから取得する。config_tagetdb.toml の input.fn_ppcList にリストした pointing それぞれについて、その周辺のデータを取得する。それらを `merge_target_csv.py` で一つのファイルに結合する。

```bash
# リモートのカタログサーバーからデータを取得する
python get_targetdb_cat.py config_targetdb_cosmos.toml

# 取得したデータをマージする
python merge_target_csv.py ./cosmos/sky
python merge_target_csv.py ./cosmos/fluxstd
python merge_target_csv.py ./cosmos/gaia
```

LSST の Deep Drilling Field である、COSMOS と XMM_LSS の領域については、PFS ７視野でカバーする際に必要になるデータはすでにダウンロード(./cosmos、./xmm_lss ディレクトリ)してあり、そのまま使用できます。

もしデータが存在しない環境で動かす場合は、GitHub Release からアーカイブされたデータを一括ダウンロードして解凍するスクリプト `download_data.sh` を利用できます。

```bash
# GitHubのリリースからデータをダウンロードして展開する
bash download_data.sh
```

---

### 0.4. データベース接続設定 (`db_config.toml`)

pfsDesign 生成や一部のデータベース連携時に使用する、認証情報や接続情報（SSH・PostgreSQL）は、ソースコード上にハードコードせず db_config.toml から動的に読み込まれます。

* この設定ファイルにはパスワードやSSH鍵情報などの機密情報が含まれるため、`.gitignore` によって Git へのコミットから自動的に除外されます。
* 環境変数 `DB_CONFIG_FILE` がセットされている場合はそのパスの設定ファイルが優先され、指定がない場合は実行ディレクトリ内の `db_config.toml` をデフォルトで読み込みます。

---

### 0.5. 視野の最適化

準備したターゲットリストから、指定した priority (default 2) までの天体を、指定した視野数で最大数カバーできるように、かつ各ガイドカメラに規定の数以上のガイド星が入るように、視野の中心座標と回転角 (PA) を最適化します。PFSの視野は六角形を仮定し、ガイド星の判定には Gaia カタログと座標変換を使用します。結果は `optimized_pointings.ecsv` に保存されます。

```bash
python optimize_hex_fov_with_guidestars.py --input ./cosmos/targets_all_20260514.csv --gaia-catalog ./cosmos/gaia.ecsv --max-priority 2 --num-fovs 4
```

以下の`run_netflow.py`で、観測視野のリストを指定しなかった場合は、このスクリプトが呼び出されて、最適化された視野リストが生成され使用される。

---

### 0.6. 視野の局所最適化（ローカルサーチ）

既存のポインティングリスト（例: `hexagons_cosmos_flat_centers.ecsv`）に対して、各座標でガイド星制約（必要なガイド星の数や、明るすぎる星の除外など）が満たされているかを確認し、満たしていない場合はその座標の周辺（RA, Dec, PA）を探索して制約をクリアできる最も近い座標を探し出します。最適化されたリストは新しいファイルに保存されます。

```bash
python optimize_hex_fov_local_search.py --input hexagons_cosmos_flat_centers.ecsv --output hexagons_cosmos_flat_centers_opt.ecsv
```

探索範囲やステップサイズは `--search_radius`, `--search_step`, `--pa_radius`, `--pa_step` 引数で調整可能です。
また、`--avoid-gaps` オプションを指定することで、隣接する視野同士に隙間ができないように（理想的な配置間隔から離れすぎないように）ペナルティを課して調整座標を選ぶことができます。

## 1. 実行方法

アクティベートした仮想環境の Python を使用して、あるいは作成された `.venv/bin/python` を明示的に指定して、スクリプトを実行します。

```bash
# 仮想環境を有効化した状態で実行
python run_netflow.py --config netflow_pipeline_config.yaml

# または、作成した仮想環境の python インタプリタを直接指定して実行
.venv/bin/python run_netflow.py --config netflow_pipeline_config.yaml
```

### コマンドライン引数
* `-c`, `--config` : パイプライン全体を制御する YAML 設定ファイルのパスを指定します（デフォルト: `netflow_pipeline_config.yaml`）。
* `--config-yaml` : 観測ランごとに提供される Gurobi / PFS ソフトウェアパラメータの設定ファイル（例: `config.yaml`）のパスを指定します（デフォルト: `config.yaml`）。

### 注意事項: config.yaml からのパラメータ自動マージについて
`run_netflow.py` を実行する際、指定された `config.yaml` の内容が動的に読み込まれ、`netflow_pipeline_config.yaml` の設定を上書き・補完します。

* **Gurobi パラメータ**: `config.yaml` 内の `gurobi.param` セクションに記述されたパラメータ（`seed`, `method`, `mipgap`, `threads` 等）が、ソルバー実行時に自動適用されます。
* **PFS ソフトウェアパラメータ**: `config.yaml` 内の `pfs` セクションに記述された以下のパラメータが、ソルバーやバリデーション実行時に自動反映されます。
  * `black_dot_radius_margin` (ベンチ/黒点回避マージン。Ffi/Designのバリデーション時にも連動して反映されます)
  * `brokenCobrasMargin` (故障Cobra回避マージン)
  * `fiducialsAvoidDistance` (Fiducial回避距離)
  * `dot_penalty` (黒点近傍配置ペナルティ)
  * `cobraSafetyMargin` (Cobra安全マージン)
  * `numReservedFibers` (予約ファイバー数)
  * `fiberNonAllocationCost` (ファイバー非割り当てペナルティ)

---

## 2. 主要な機能と動作フロー

`run_netflow.py` を実行すると、以下のフローに従って自動的に処理が実行されます。

### ① 設定ファイルと天体カタログの読み込み
指定された設定ファイル（YAML/TOML）を読み込み、観測プログラム情報、入力天体パス、出力ディレクトリなどを特定し、対象となる `science`, `fluxstd`, `sky` の各天体カタログをメモリ上にロードします。

### ② 視野 (Pointing) の自動最適化 (FoV Optimization)
設定ファイルの `inputs.pointing_file` が `null` の場合、**視野自動最適化機能**が自動的にトリガーされます。
* `optimize_hex_fov_with_guidestars.py` 内のアルゴリズム（ガイド星制約付き Greedy Maximum Coverage）を使用し、ロードされた `science_targets` のうち `priority <= max_priority` を満たすターゲット天体を最も効率よくカバーしつつ、6つのガイドカメラそれぞれに十分な数（デフォルト: 各2つ以上）のガイド星が入るように視野の中心座標と回転角 (PA) を探索・最適化します。
* 指定された視野数 (`num_fields`) のポインティング情報を算出し、自動的に `optimized_pointings.ecsv` として保存します。このファイルが後続の処理の入力ポインティングとなります。

### ③ 露出時間の一時上書きとアサイン計算
netflow ソルバー（Gurobi）による割り当て計算を行う直前に、以下の処理を安全に実行します。
* 観測対象の `science` 天体の本来の必要露出時間 (`_obs_time`) を一時的に退避します。
* 各天体の露出時間を、設定ファイルに定義された1回あたりの観測露出時間 `t_obs` (デフォルト: `900.0`s) に一時的に上書きします。これにより、カタログ上の設定と切り離した露出条件での最適配置が計算されます。
* この最適化では、すべての天体を一回の観測で終了するとして最適化を行います。観測の途中でファイバー配置を変えることは、このステップでは考慮されていません。
* **SFAに基づくスカイの空間均等化制約の適用**: 焦点面上を20個の幾何学的なセクター（Laszlo領域）に分割し、それぞれの領域に対して「最低12本のスカイファイバーを必ず割り当てる」という強いペナルティ制約（`locationGroupPenalty=1e11`）をソルバーに課します。これにより、サイエンスターゲットの割り当て数を実質的に一切犠牲にすることなく（十分な余剰ファイバーを活用して）、空間的に均等で高品質なスカイ観測が保証されます。
* 割り当て計算（`solve_assignment`）が正常に終了、または万が一例外で中断された場合でも、`try...finally` ブロックによって天体データの `_obs_time` を**確実に元の露出時間へ復元**します。

### ④ アサイン結果の出力とプロット
アサインされた天体の組み合わせ結果を `targets/` ディレクトリ配下に ECSV 形式（`science`, `fluxstd`, `sky`）で出力し、ターゲット全体のスカイ分布図（プロット画像）を生成します。

### ⑤ pfsDesign (FITS) / OPE 生成と Validation
* 生成されたアサイン結果を基に make_pfs_design.pyを呼び出し、`pfsDesign` FITS ファイルおよび観測実施用の `.ope` ファイルを生成します。
* 設定の `pointing_file` が `null` の場合、自動的に ② で生成された `optimized_pointings.ecsv` を参照するようにフォールバックします。
* バリデーション実行時には、外部の PostgreSQL 接続（HiloのGaia DB）を使用する代わりに、ローカルの ECSV カタログ（`gaia_catalog`）を直接読み込んで座標と等級で高速にフィルタリング・モッククエリすることで、明るすぎる星のチェックやガイド星の数のチェックなどを行うローカル Validation を実行します。

---

## 3. 設定ファイルの設定項目 (`netflow_pipeline_config.yaml`)

netflow_pipeline_config.yaml の主要な設定パラメータは以下の通りです。

```yaml
# 共通設定
proposal_id: "S26A-104"
obstime: "2026-05-09T06:00:00Z"              # 観測日時(UTC)

# 入力設定
inputs:
  pointing_file: null                        # 既存のポインティングファイル、null指定で自動最適化を実行
  catalog_dir: "./cosmos/"
  science_targets: "targets_all_20260514.csv"
  fluxstd_targets: "fluxstd.ecsv"
  sky_targets: "sky.ecsv"
  gaia_catalog: "./cosmos/gaia.ecsv"            # Validation用ローカルGaiaカタログ

# 出力設定
outputs:
  targets_dir: "targets"                     # アサイン結果ECSVの出力先フォルダ
  plot_file: "sky_distribution.png"          # 最終アサイン結果のプロットファイル
  fov_plot_file: "fov_coverage.png"           # 視野自動最適化時のカバー分布図プロットファイル
  text_output: "output.txt"

# netflow計算パラメータ
netflow:
  t_obs: 900.0                              # 露出上書き用の観測時間 (秒)
  num_fields: 4                              # ポインティング自動最適化時の視野数 (FoVの配置数)
  max_priority: 2                            # ポインティング自動最適化で考慮する最大優先度 (小さい値ほど高優先)
  min_stars_per_cam: 2                       # 1つのガイドカメラ内に必要なガイド星の最小数
  min_cams_with_stars: 6                     # min_stars_per_cam の条件を満たすべきガイドカメラの最小数（最大6）
```

---
