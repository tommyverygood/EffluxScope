# Frozen-score benchmark

本目录只运行回顾性计算，不训练模型，不修改旧数据、旧模型或论文。使用同版本整张输入表，v2主分析、v1敏感性分别运行。

在项目根目录运行：

```sh
.venv-experiments/bin/python outputs/20261009_adc_pair_experiments/benchmark/run_benchmark.py
.venv-experiments/bin/python outputs/20261009_adc_pair_experiments/benchmark/plot_benchmark.py
```

`run_benchmark.py` 只需要 Python 与 NumPy；`plot_benchmark.py` 另外需要 Matplotlib。脚本根据自身位置定位项目输入，所有输出写在脚本同目录。准确依赖版本与输入 SHA256 分别记录于 `protocol.json` 和 `figure_manifest.json`。已有 runtime 可运行；迁移给其他人时，先安装记录的 NumPy/Matplotlib 版本并提供两个原始 analysis CSV，不能只拷贝PNG当作可复算数据。

验收条件：`input_audit.json` 中 `passed=true`；全v0.4必须重现453个accepted、TP221/TN119/FP84/FN29、BAcc 0.7351034483；同coverage contrasts的两侧接受数相同；2000次共同wholecluster抽样，差值在同一draw内计算；full pool不足时目标unavailable；旧输入哈希未变。

阅读 `实验报告.md` 及 `paired_contrasts.csv`。同coverage下full v0.4相对辅助v0.3 own margin的负结果必须保留。no-local全部551个与同coverage重新选出的453个名单不同，不能混写结果。v1与v2相同不构成独立重复验证；本实验也不是ADC released-species配对湿实验。
