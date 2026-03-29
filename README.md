# TCGA Lung Multi-Omics

TCGA-LUAD + TCGA-LUSC 多體學整合分析 pipeline，基於 [pymaftools](./pymaftools)。

## Pipeline 概覽

```
generate_manifests  →  align_manifests  →  download  →  build_tables  →  filter_tables  →  build_cohort
      (GDC API)           (取交集)         (gdc-client)   (HDF5 tables)    (Primary Tumor)     (Cohort)
```

資料類型（見 `config.toml`）：Expression · Mutation · CNV (ASCAT3) · Methylation

## 快速開始

```bash
# 安裝依賴
uv sync

# 查看目前進度
uv run python pipeline.py status

# 跑完整 pipeline
uv run python pipeline.py run

# 從某一步繼續
uv run python pipeline.py run --from download

# 只跑特定步驟
uv run python pipeline.py run --steps build_tables filter_tables build_cohort
```

## 設定

編輯 `config.toml`：

```toml
projects = ["TCGA-LUAD", "TCGA-LUSC"]

[download]
gdc_client = "tools/gdc-client"      # gdc-client 路徑
token      = "/path/to/token.txt"    # GDC token（受控資料需要）
threads    = 20
retries    = 5

[data_types.expression]
data_type     = "Gene Expression Quantification"
workflow_type = "STAR - Counts"
# ... 其他 data types
```

## Portal 模式（替代方案）

若從 [GDC Data Portal](https://portal.gdc.cancer.gov/) 手動下載 manifest：

```bash
# 必須包含所有 data type（expression、mutation、cnv、methylation）
uv run python scripts/align_manifests.py --mode portal --input data/gdc_manifest.txt
# 缺少任何 data type 時會立即報錯提示
```

## 各步驟說明

| 步驟 | 腳本 | 輸出 |
|------|------|------|
| generate_manifests | `scripts/generate_manifests.py` | `data/manifests/full/` · `data/file_to_case.tsv` |
| align_manifests | `scripts/align_manifests.py` | `data/manifests/aligned/` · `data/aligned_cases.tsv` |
| download | `scripts/download.py` | `data/raw/` · `data/download_log.tsv` |
| build_tables | `scripts/build_tables.py` | `data/processed/tables/` |
| filter_tables | `scripts/filter_tables.py` | `data/processed/filtered/` |
| build_cohort | `scripts/build_cohort.py` | `data/processed/tcga_lung.h5` |

各腳本均可獨立執行，支援 `--types` 等參數，詳見各腳本 `--help`。

## 下載 gdc-client

```bash
wget https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip
unzip gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip
unzip gdc-client_2.3_Ubuntu_x64.zip -d tools/
chmod +x tools/gdc-client
```
