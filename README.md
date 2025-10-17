# Reconstruct_the_sample_surface_based_on_APT_data_and_quantify_the_surface_enrichment
Reconstruct the Sample Surface from APT Data & Quantify Surface Enrichment (MATLAB)

For the specific principle, please refer to the document below.

[Readme.docx](https://github.com/user-attachments/files/22966786/Readme.docx)

---

## ✨ What this does

* Builds a **HEA-only surface envelope** from APT points (Ir/Ru/Rh/Pt/Pd).
* Performs **grid-wise rank peeling**: in each XY cell, atoms are ranked by signed normal distance to the surface; layer-1 = topmost, layer-2 = next, etc.
* Lets you **extract any layer on demand** and compute per-layer stats (mean/median Z, counts, composition).
* Quantifies **Ir-centered local enrichment** of neighbors via KDE PDFs, posterior composition (P(s\mid r)), and enrichment (E_s(r)=P(s\mid r)/\pi_s).
* Includes a **persistent interactive app**: choose “surface = top N layers”, center element, neighbor scope, radius, baseline — and export results.

---

## 🧩 Inputs

1. **APT point cloud** (from AP Suite 6 export):
   CSV with at least these columns

   * `X`, `Y`, `Z` — coordinates (nm)
   * `m/z` (or similar; script auto-detects)

2. **m/z windows** (from AP Suite 6 export):
   Excel or CSV with columns

   * `Species`, `mz_min`, `mz_max`
     …used to map each ion to a species.

> The surface is constructed from the **HEA set**: `Ir, Ru, Rh, Pt, Pd` (configurable in code).

---

## ⚙️ Requirements

* MATLAB **R2021a or later**
* **Statistics and Machine Learning Toolbox** (uses `createns`, `rangesearch`, `ksdensity`, `groupsummary`)

Tested on R2021b–R2024a.

---

## 🚀 Quick start

1. Clone the repo and open MATLAB.
2. Put your APT CSV and m/z window file in the same folder (or note their paths).
3. Run the main script:

   * `ExtractSurf_HEA_Layers_App.m`
     *(also compatible with `ExtractSurf_EnglishVersion.m` if present; the app/UI is in the “createSurfaceApp” section).*
4. Choose parameters when prompted:

   * **Grid step** (h), **quantile** (q) for the upper envelope (e.g., 0.98), **min HEA points/cell** (N_{min}), and **Layer-1 definition** (`topmost` or `closest`).
5. Use the **interactive app** to define *surface = top N layers*, pick the **center element**, **neighbor scope** (HEA-only vs all), **radius** (R), and **baseline** (global vs surface).
   The app computes enrichment and plots PDFs/Posterior/Enrichment; tables/figures can be exported.

---

## 📊 Output

### Figures

* **Surface mesh** (HEA-only envelope)
* **Layer-wise HEA scatters**: layers 1…5 (or up to available), element-colored (Ir/Ru/Rh/Pt/Pd)
* **Combined top-N** layers scatter
* **Ir-centered pair-distance analysis**: KDE PDFs, posterior (P(s\mid r)), enrichment (E_s(r))

### Excel workbook (default: `HEA_surface_output.xlsx`)

* `SurfaceGrid` — sampled surface ((X,Y,Z))
* `SurfaceLayer1` … `SurfaceLayer5` — HEA atoms per peeled layer
* `LayerZ_Stats` — per-layer counts & Z stats (Mean/Median/Std/Min/Max)
* `LayerD_Stats` — per-layer normal-distance stats
* `Layer*_Summary` — per-layer composition (counts & fractions)
* `PDF_KDE_HEA`, `Enrichment_HEA` — Ir-centered KDE & enrichment tables
* (If enabled) `ShellProb`, `ShellCounts`, `ShellMeta` — FCC shell summaries

---

## 🧠 Method (very briefly)

* **Surface envelope:** Partition XY into a grid (step (h)). In each occupied cell, take the **(q)-quantile** of Z (e.g., (q=0.98)) as the local surface height, and fit a local plane for the **surface normal**.
* **Signed normal distance:** For each atom, compute (d = \mathbf{n}\cdot(\mathbf{x}-\mathbf{x}_\text{surf})). Positive (d) points into the material.
* **Grid-wise rank peeling:** In each cell, **rank HEA atoms** by (d) (descending for `topmost` or ascending for `closest`). The 1st ranked atom is **Layer-1**, the 2nd is **Layer-2**, etc. This guarantees layers beyond 1 are populated wherever enough atoms exist.
* **Enrichment:** For center species (A) in the selected surface (top (N) layers), estimate neighbor PDFs (f_s(r)) for species (s) with KDE. With baseline fractions (\pi_s) (global or surface), posterior (P(s\mid r) \propto f_s(r)\pi_s), enrichment (E_s(r)=P(s\mid r)/\pi_s).

---

## 🔧 Key parameters

* `h` — XY grid step (nm)
* `q` — surface quantile (0–1)
* `Nmin` — min HEA points per cell
* `layer_order` — `topmost` (by largest (d)) or `closest` (by smallest (d))
* In the app: **N layers**, **center**, **neighbor scope**, **R**, **baseline**

---

## 📎 Citation

If you use or refer to this code, please cite:

> Qu, M. (2025). *Reconstruct the sample surface from APT data and quantify the surface enrichment* [GitHub repository].
> [https://github.com/QuMingrong/Reconstruct_the_sample_surface_based_on_APT_data_and_quantify_the_surface_enrichment](https://github.com/QuMingrong/Reconstruct_the_sample_surface_based_on_APT_data_and_quantify_the_surface_enrichment)

You may also include the link directly in figure captions or methods sections, e.g.:

> “MATLAB code for APT surface reconstruction and grid-wise layer peeling is available at
> [https://github.com/QuMingrong/Reconstruct_the_sample_surface_based_on_APT_data_and_quantify_the_surface_enrichment.”](https://github.com/QuMingrong/Reconstruct_the_sample_surface_based_on_APT_data_and_quantify_the_surface_enrichment.”)

---

## 📨 Contact

Questions or collaborations: **[qmr@mail.ustc.edu.cn](mailto:qmr@mail.ustc.edu.cn)**
