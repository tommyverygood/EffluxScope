# Current results figures and TOC

This component redraws the approved 2026-11 manuscript Figures 4, 5 and 6, and the TOC graphic, from frozen plotting inputs. Legacy generator/output names are preserved to keep the original scientific lineage traceable.

| Manuscript item | Entry point | Output |
|---|---|---|
| Figure 4, source diagnostics | `论文主稿/latex_jcim/tools/make_figure3_jcim.py` | `论文主稿/latex_jcim/figures/Figure3_validation_gates.pdf` |
| Figure 5, structural splits | `论文主稿/latex_jcim/tools/make_figure2_jcim.py` | `论文主稿/latex_jcim/figures/Figure2_split_separation.pdf` |
| Figure 6, errors and abstention | `论文主稿/latex_jcim/tools/make_figure4_jcim.py` (delegates to `论文主稿/figures_src/make_figure4.py`) | `论文主稿/latex_jcim/figures/Figure4_error_distance.pdf` |
| TOC | `toc/tools/make_toc_graphic_wheel.py` | `toc/toc_graphic/TOC_graphic_wheel.pdf` |

Run `python run_results_figures.py`, then `python verify_results_figures.py`. Both commands resolve their paths relative to this directory and can be invoked from another working directory. The plotting command writes PDFs and additional print/raster formats; the verification command checks PDF text, page dimensions and 150 dpi Poppler rendering against the included approved reference figures.

## Dependencies and verification

Verified in a fresh Python 3.12 environment with NumPy 2.5.1, pandas 2.2.3, Matplotlib 3.11.2, Pillow 12.3.0, RDKit 2026.3.4 and pypdf 7.0.0. The host supplies Ghostscript for EPS export, Poppler `pdftoppm` for pixel comparison, and Arial for the approved typography. No font files are redistributed. Exact figure comparison requires the separately installed Arial family; another font may preserve the data but change text geometry. The entry point fails clearly if Arial is absent.

All four outputs passed identical extracted text, page dimensions and rendered pixels (zero changed pixels) in `verification/redraw_comparison.json`. PDF hashes can differ because creation timestamps and other metadata differ. No model, molecular embedding, structural partition or confidence interval was refitted by this redraw.

## Frozen checks

Figure 4 checks the ASAP association and criterion, identity-audited Roth GI50 resistance diagnostic, all eight regulatory criteria, five regulatory metrics, and full-coverage values against the bundled frozen result files. The regulatory pool is pre-score-locked; this is not a claim of outcome blinding.

Figure 5 checks seven numerical/coordinate input hashes, 2,624 external structures, 443 reference coordinates and the original split memberships (62 scaffold-removed and 544 KMeans-OOD compounds). The two maps share the same frozen density field. Its `render.json` is generated output and is not a new statistical analysis.

Figure 6 checks aggregate confusion-count arithmetic for ten inside/outside model comparisons and five full-set balanced accuracies. In this component it loads the archived intervals. It does not need or redistribute third-party individual prediction scores. Interval recomputation is a separate workflow in the complete release and must be evaluated separately from this redraw.

The TOC asserts the archived matched-count balanced-accuracy ordering before drawing its qualitative rank labels.

## Export changes and provenance

The original scientific inputs and approved graphical values are retained. The only data-document change is replacing a machine-specific `protocol_lock_path` in `数据/processed/abcb1_regulatory_v4_blind_summary_v1.json` with the equivalent project-relative path; this string is historical provenance, not a required redraw dependency. The old Figure 7 wording in the Figure 6 wrapper docstring was corrected. Original file hashes are retained in `original_input_manifest.json`; exported hashes and these changes are listed in `export_provenance.json`.

The source/aggregate files were already represented in the public 2026-10 study and correction releases. The latest plotting code, shared palette, density display and paired-bar presentation are the updates supplied here. Apply the complete release's source-specific licence notices: author-written software uses MIT; source-derived data retain their original attribution and licence terms. Approved reference graphics are author-produced outputs, not third-party article figures. This component includes no manuscript full text, credentials, private correspondence, third-party article full texts, or model checkpoints.
