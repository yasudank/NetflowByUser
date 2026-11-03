# P0 天体を最大化するポインティング最適化と複数回観測の割り当て
(`optimize_pointings_p0.py` と `run_netflow.py` の `prefer_visits` オプション)

本ドキュメントでは、[optimize_pointings_p0.py](../src/netflow_by_user/optimize_pointings_p0.py) の目的・アルゴリズム・使い方と、複数パターン（複数回）観測のための netflow の割り当てオプション `prefer_visits` について解説します。

---

## 1. 概要

`optimize_pointings_p0` は、敷き詰められたポインティング群（例: `hexagons_xmm_lss_flat_centers.ecsv` の 7 視野）を基準として、**実際にファイバーを割り当てられる最優先天体 (P0) の数**が最大になるように、各視野の中心座標と位置角 (PA) を最適化します。

* 視野を六角形として数えるのではなく、netflow と同じ装置モデル（黒点、故障 Cobra、フィデューシャル）で **各 Cobra が届く天体** を計算して評価します。
* ガイド星の制約（`min_stars_per_cam`, `min_cams_with_stars`）と、故障ファイバー近傍の輝星の制約を満たさない視野は天体を観測できないものとして扱い、その数を最優先で 0 にします。
* **複数パターン**（同じ領域を覆う敷き詰めを K 組、順に観測する）にも対応します。総観測時間を K 等分して、「1 回以上観測される P0 の数」を最大化したうえで、観測回数ができるだけ揃うように配置します。

---

## 2. 評価方法

### 2.1. 1 視野の評価
候補となる視野（RA, Dec, PA）ごとに、次を計算してキャッシュします。

1. ガイドカメラ AG0〜AG5 ごとのガイド星数と、制約を満たすかどうか
2. 故障ファイバー近傍の輝星の有無
3. 対象天体（中心から `--target-radius` 以内、優先度 `--max-priority` 以下）の焦点面座標と、それぞれの天体に届く Cobra の一覧

この計算は `--processes` 個のプロセスで並列に実行されます。

### 2.2. 割り当て可能数（最大流）
全視野の「天体 → Cobra」の辺から次のネットワークを作り、最大流で P0 の割り当て可能数を求めます。

```
source → 天体 (容量 j) → (天体, パターン) (容量 1) → Cobra (容量 1) → sink
```

* **1 パターンにつき 1 天体 1 本**: netflow はパターンごとに別々に実行されるため、隣り合う視野の重なりで同じ天体に 2 本割り当てても、観測回数は 1 回にしかなりません。この制約を入れることで、予測値が netflow の結果と一致します。
* ファイバー衝突と較正天体（sky, fluxstd）は考慮しないので、厳密には上限値です。XMM-LSS での検証では、P0 については netflow の結果と 1 天体単位で一致しました。

### 2.3. スコア（辞書式比較）
レイアウトは次のタプルを辞書式に比較して評価します（左ほど優先）。

1. ガイド星などの制約を満たさない視野の数（少ないほど良い）
2. 1 回以上観測される P0 の数、2 回以上の数、…、K 回以上の数（K はパターン数）
3. P1〜P`max_priority` で届く天体の数（タイブレーク）
4. 基準の敷き詰めに対する隙間の割合（小さいほど良い）

「j 回以上観測される P0 の数」は、各天体の容量を j にしたときの最大流 T_j を使って T_j − T_{j−1} として求めます。割り当て可能なファイバー数ベクトルはポリマトロイドをなすので、1 つの割り当てがすべての j について Σ min(n_t, j) を同時に最大化します。このため、これらの値は与えられた視野に対して厳密です。

---

## 3. 探索アルゴリズム

### ステップ 1: 大域探索
基準の敷き詰め全体を剛体として回転（`--global-rot-step` 刻みで 0〜60°）および平行移動（±`--global-shift`, `--global-shift-step` 刻み）させます。さらに各視野を 60° の倍数だけ回転させます（六角形なので敷き詰めは崩れません）。60° 倍数の選び方は座標降下法で決めます。

### ステップ 2: 局所探索
大域探索の上位 `--n-refine` 個から出発し、各視野を個別に小さく動かして（±`--local-shift`, ±`--local-rot`）座標降下法でスコアを改善します。隙間の割合が `--max-gap` を超える配置は採用しません。制約を満たしていない視野から先に動かします。

### ステップ 3: 複数パターンの場合
1. パターンを 1 つずつ追加します。既存のパターンは固定し、新しいパターンを上の 2 段階で最適化します（`--init` を与えたパターンは探索せず、そのレイアウトから始めます）。
2. その後、各パターンを他のパターンを固定したまま順に再最適化します（`--cycles` 回まで、改善がなくなれば終了）。このときは現在のレイアウトと大域探索の上位 `--cycle-n-refine` 個から局所探索します。現在のレイアウトより良い場合のみ置き換えるので、スコアは悪くなりません。

---

## 4. 使い方

`pip install -e .` でパッケージをインストールすると `optimize_pointings_p0` コマンドが使えます（`python src/netflow_by_user/optimize_pointings_p0.py` でも実行できます）。対象天体・Gaia カタログ・ガイド星の制約などは `--config` の YAML（`run_netflow` と同じ設定ファイル）から読み込みます。

### 4.1. 1 パターン
```bash
optimize_pointings_p0 -c netflow_config_2026Nov.yaml \
    --grid xmm_lss/hexagons_xmm_lss_flat_centers.ecsv \
    -o xmm_lss/hexagons_xmm_lss_flat_centers_p0opt.ecsv \
    --global-shift 0.6 --global-rot-step 10 --processes 20 \
    --add-columns --report-file xmm_lss/p0opt_report.txt --plot xmm_lss/p0opt.png
```

### 4.2. 既存のパターンに 1 パターン追加する（2 パターン観測）
既存のパターンを `--fixed` で固定し、新しいパターンのコードには `--code-suffix` で接尾辞を付けます。
```bash
optimize_pointings_p0 -c netflow_config_2026Nov.yaml \
    --grid xmm_lss/hexagons_xmm_lss_flat_centers.ecsv \
    --fixed xmm_lss/hexagons_xmm_lss_flat_centers_p0opt.ecsv --code-suffix _B \
    -o xmm_lss/hexagons_xmm_lss_flat_centers_p0opt_B.ecsv \
    --global-shift 0.6 --global-rot-step 10 --n-refine 3 --processes 20 \
    --add-columns --report-file xmm_lss/p0opt_B_report.txt --plot xmm_lss/p0opt_B.png
```
以前の結果を `--init` に与えると、それを出発点の 1 つとして含めます（結果がそれより悪くなることはありません）。

### 4.3. K パターンをまとめて最適化する
```bash
optimize_pointings_p0 -c netflow_config_2026Nov.yaml \
    --grid xmm_lss/hexagons_xmm_lss_flat_centers.ecsv \
    --n-patterns 6 --init xmm_lss/hexagons_xmm_lss_flat_centers_p0opt.ecsv \
    --global-shift 0.6 --global-rot-step 10 --n-refine 3 --cycles 4 --processes 20 \
    -o xmm_lss/hexagons_xmm_lss_flat_centers_p0opt_6pat.ecsv \
    --add-columns --report-file xmm_lss/p0opt_6pat_report.txt --plot xmm_lss/p0opt_6pat.png
```
出力の `ppc_code` には `_A`, `_B`, … が付き、`pattern` カラムにパターン名が入ります。

### 4.4. 実行時間の目安
20 プロセスで、1 パターンの追加が数分〜十数分、4 パターンで約 20 分、6 パターンで 30 分〜1 時間程度です。評価結果はパターン間・サイクル間で共有されます。計算は Python から呼ぶ装置モデルが中心なので、GPU では速くなりません。

---

## 5. コマンドライン引数

| 引数 | デフォルト | 説明 |
|---|---|---|
| `--config`, `-c` | (必須) | パイプライン設定 YAML（天体カタログ、Gaia カタログ、ガイド星制約など） |
| `--config-yaml` | `config.yaml` | PFS パラメータを含む config.yaml |
| `--grid`, `-i` | (必須) | 基準の敷き詰め ECSV（`ppc_code`, `ppc_ra`, `ppc_dec`[, `ppc_pa`]） |
| `--output`, `-o` | (必須) | 出力するポインティング ECSV |
| `--targets` | 設定から | 科学天体 CSV |
| `--gaia` | 設定から | Gaia カタログ ECSV |
| `--center RA DEC` | 基準の平均に最も近い視野 | 対象領域の中心 |
| `--target-radius` | 1.8 | 中心からこの半径（度）以内の天体のみを使う |
| `--max-priority` | 2 | 読み込む天体の最大優先度（P0 を最大化し、他はタイブレーク） |
| `--n-patterns` | 1 | 最適化するパターン数 |
| `--fixed` | なし | 固定して追加で観測するパターンの ECSV（複数可） |
| `--init` | なし | 最初のパターンの初期レイアウトとして使う ECSV（複数可） |
| `--cycles` | 2 | 各パターンを順に再最適化する最大回数（複数パターン時） |
| `--cycle-n-refine` | 1 | 再最適化で局所探索する大域探索の上位解の数 |
| `--global-shift` | 0.2 | 敷き詰め全体の最大平行移動（度） |
| `--global-shift-step` | 0.1 | 平行移動の刻み（度） |
| `--global-rot-step` | 5.0 | 全体回転の刻み（度、0〜60° を探索） |
| `--n-refine` | 5 | 局所探索する大域探索の上位解の数 |
| `--local-shift` | 0.08 | 各視野の最大移動（度） |
| `--local-shift-step` | 0.02 | 局所移動の刻み（度） |
| `--local-rot` | 5.0 | 各視野の追加回転の最大値（度） |
| `--local-rot-step` | 2.5 | 局所回転の刻み（度） |
| `--max-gap` | 0.03 | 基準の敷き詰めに対して許す隙間の割合 |
| `--max-passes` | 3 | 局所探索の座標降下の最大パス数 |
| `--baseline` | なし | 比較用に評価するだけのポインティング ECSV（最適化には使わない） |
| `--code-suffix` | なし | 出力 `ppc_code` に付ける接尾辞（例: `_B`） |
| `--processes` | min(16, CPU 数) | ワーカープロセス数 |
| `--add-columns` | オフ | 出力に `ag0`〜`ag5`, `n_p0_access`, `n_p0_assigned` を追加 |
| `--report-file` | なし | レポートをテキストファイルに保存 |
| `--plot` | なし | 最適化したレイアウトの図を保存 |

### 出力
* **ECSV**: `ppc_code`, `ppc_ra`, `ppc_dec`, `ppc_pa`, `ppc_resolution`, `ppc_priority`（複数パターン時は `pattern` も）。そのまま `run_netflow` の `pointing_file` に使えます。
* **レポート**: スコア、各視野のずれ・ガイド星数・届く P0 数・割り当てる P0 数、P0 の観測回数の分布（複数パターン時）、割り当てられない P0 の一覧。
* **図**: 視野の六角形と `ppc_code`、P0 を観測回数で色分け、P1〜P2 を灰色で表示します。半径 `--target-radius` の外の P0 は白抜きで示します。

---

## 6. netflow での複数回観測の割り当て（`prefer_visits`）

最適化した K パターンを netflow で割り当てる方法は 2 通りあります。

### 6.1. パターンごとに別々に実行する
各パターンの視野だけを `pointing_file` にした設定で、`run_netflow` を K 回実行します。各パターンの割り当ては独立に決まるので、2 パターン目で同じ天体が選ばれる保証はありません。ファイバーが足りない下位優先度ほど、1 回だけ観測される天体が多くなります。

### 6.2. 全パターンをまとめて実行する（推奨）
全パターンの視野を 1 つの `pointing_file` にまとめ、設定に `prefer_visits` を加えて 1 回の `run_netflow` で解きます。
```yaml
netflow:
  t_obs: 450.0              # 1 視野の露出時間（総時間を 1 パターンと同じにするなら t_obs / K）
  prefer_visits: 2          # 2 視野以上から届く天体は 2 回観測するか、観測しない
  single_visit_weight: 0.5  # 1 視野からしか届かない天体の 1 回観測の価値（2 回観測 = 1 に対して）
```

* **2 視野以上から届く天体**: 要求露出を `prefer_visits × t_obs` にし、netflow 本来の all-or-nothing（全部観測するか、まったく観測しないか）で扱います。ファイバーが足りないときに 1 回だけ観測されることがなくなり、観測回数は「2 回 : 1 回 : 0 回 ≈ x : 0 : (1−x)」になります。
* **1 視野からしか届かない天体**（六角形の縁、Cobra の隙間、故障 Cobra など）: 1 回の観測を許し、その価値を `single_visit_weight × N` とします（N はその優先度のコスト）。既定値 0.5 では、ファイバー 1 本あたりの価値が 2 回観測と同じになります。
* 優先度の順序は保たれます（優先度ごとのコストは 3 倍ずつ異なります）。1 天体あたりの観測回数は `prefer_visits` 回までに制限されます。
* 現在の実装は `prefer_visits: 2` で検証しています。

ログには、優先度ごとに「届く視野の数（1 / 2+）× 割り当てた回数（0x / 1x / 2x）」の表が出力されます。
```
Science targets by class, accessible visits and assigned visits:
    sci_P0 accessible 1 : 0x:      1, 1x:     54, 2x:      0
    sci_P0 accessible 2+: 0x:      0, 1x:      0, 2x:    129
    ...
```
`prefer_visits` を設定しなければ、従来どおり各天体を 1 回観測として割り当てます。

### 6.3. 実例（XMM-LSS、中心から 1.8° 以内）

| 方法 | P0 (180) 2回 / 1回 / 0回 | P3 (2052) 2回 / 1回 / 0回 | P9 (23489) 2回 / 1回 / 0回 |
|---|---|---|---|
| 1 パターン | – / 153 / 27 | – / 1556 / 496 | – / 8330 / 15159 |
| 2 パターン、別々に実行 | 129 / 47 / 4 | 1161 / 725 / 166 | 3782 / 9080 / 10627 |
| 2 パターン、`prefer_visits: 2` | 129 / 47 / 4 | 1230 / 560 / 262 | 5834 / 4484 / 13171 |

1 パターンでは全露出時間、2 パターンでは 1 回の観測が半分の露出時間です。P0 の 2 回観測の 129 は、最適化で求めた「両パターンから届く数」と一致します。
