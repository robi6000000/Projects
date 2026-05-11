---
name: Protocol and notebook conventions
description: Writing conventions for protocols and results.ipynb preferences
type: feedback
originSessionId: 3823c645-9b3e-436d-8129-ba9c964f2fd8
---
Protocol writing conventions observed across sessions:

- Date-stamped markdown files in `protocols/` folder
- Concise bullet-point style, not prose
- TODOs listed at the top or in a dedicated section
- SLURM commands included verbatim with job names and array specs

Notebook conventions:
- No bold/emphasis in cell output text
- RENAME dict applied just before plotting, not globally
- Helper functions (load_roc_data, plot_roc, plot_auc_bars) in the imports cell
- File-existence guards (`os.path.exists`) for cells that depend on incomplete jobs
- Save figures to `scripts/analysis/` with numbered filenames

**Why:** thesis-ready outputs, clean separation of data cells and plot cells.
**How to apply:** always read the current notebook state before editing; don't assume cell order from memory.
