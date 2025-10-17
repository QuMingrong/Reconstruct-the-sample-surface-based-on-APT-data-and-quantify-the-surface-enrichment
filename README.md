# Reconstruct_the_sample_surface_based_on_APT_data_and_quantify_the_surface_enrichment
MATLAB script for Atom Probe Tomography (APT) data of High-Entropy Alloys (HEA: Ir/Ru/Rh/Pt/Pd). It reconstructs the HEA surface from HEA atoms only, performs per-cell ranked layer peeling so you can view any layer on demand, and computes 1-nm enrichment with KDE curves. Everything is exportable to Excel and controllable from an interactive UI.

Features
	Surface from HEA atoms only (upper envelope + local normals)
	Per-cell ranked peeling: in each XY grid cell, HEA atoms are sorted by signed normal distance d; the 1st/2nd/… picks form Layer-1/Layer-2/…
	On-demand layers: view/export any layer; per-layer Z statistics (mean/median/std/min/max)
	Interactive app: choose “surface = top N layers”, center element, neighbor universe, baseline (global/surface), and radius R; compute instantly
	1-nm enrichment (or any R): composition P, baseline π, enrichment E=P/π
	KDE curves: pair-distance PDF f_s(r), posterior P(s|r), enrichment E_s(r)
	(Optional) FCC shell reference lines
	Clear visualizations (surface mesh, element-colored scatters, curves) and Excel exports
________________________________________
Input Data (exported from AP Suite 6)
You need two files:
	APT atoms table (*.csv or *.xlsx) with at least:
	X, Y, Z (nm)
	m/z (or mz, MassToCharge, etc.)
Example:
X,Y,Z,m/z
5.312,-2.441,-6.023,101.905
...
	m/z window table (APT_mz_windows.xlsx or CSV) with exact columns:
	Species (e.g., Ir, Ru, Rh, Pt, Pd)
	mz_min
	mz_max
Example:
Species,mz_min,mz_max
Ir,101.5,102.4
Ru,99.2,100.1
Rh,102.7,103.6
Pt,194.5,195.6
Pd,105.5,106.4
Multiple windows may map to the same Species (isotopes/molecular ions).
________________________________________
Requirements
	MATLAB R2021b+ (recommended)
	Statistics and Machine Learning Toolbox (ksdensity, createns, rangesearch)
	Windows/macOS/Linux
	Excel export via writetable (no extra dependency)
________________________________________
Quickstart
	Clone this repo.
	In MATLAB, run: ExtractSurf_HEA_Layers_App.m.
	Choose:
	the APT atoms file (x,y,z,m/z),
	the m/z windows file,
	an output Excel path.
	Parameter dialog (defaults are good starting points):
	h (grid step, nm): 0.3
	q (upper-envelope quantile 0–1): 0.98
	Nmin (min HEA points per cell): 20
	Layer-1 definition: topmost (largest d) or closest (smallest d)
	The script builds the surface, ranks HEA atoms per cell, exports first layers and stats, then opens the interactive app for further analysis.
________________________________________
What the App Does
	Surface = top N layers (choose N at any time)
	Center element (default from available HEA species)
	Neighbor mode
	HEA (exclude center) — only HEA neighbors except the center element
	All species (exclude center) — all species except the center element
	Baseline
	Global baseline — global atomic shares within the selected universe
	Surface baseline — shares computed only from atoms in “surface = top N layers”
	Radius r (nm) — default 1.0
On “Compute 1 nm enrichment”, the app produces:
	table with PairsWithinR, P_withinR, π_baseline, E=P/π, and a tendency label (Aggregate/Repel/Neutral),
	three curves: PDF, Posterior, Enrichment.
The Per-layer Z stats table updates accordingly.
________________________________________
Excel Outputs (sheets)
	SurfaceGrid: surface mesh points X,Y,Z
	SurfaceLayer1, SurfaceLayer2, …: HEA atoms per layer (X,Y,Z,Species,ColorHex)
	LayerZ_Stats: Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ
	LayerD_Stats: per-layer signed normal distance d (MeanD, MedianD, StdD)
	Layerk_Summary: species counts and fractions for Layer k
	PDF_KDE_HEA, Enrichment_HEA: example Ir-center HEA-neighbor KDE and enrichment long tables
	From the app: SurfN{N}_Enrich_{Center}, SurfN{N}_PDF_{Center}, SurfN{N}_Post_{Center}, SurfN{N}_Enrich_{Center}
If the Excel file is open, writing may fail—close it and retry.
________________________________________
Method (Short)
	Surface from HEA only
In XY cell C_ij, take upper-envelope height
z_s (i,j)=Q_q ({z_k∈C_ij})

fit local plane z≈ax+by+c; unit normal
n ̂=((-a,-b,1))/√(a^2+b^2+1).
	Signed normal distance
For point p=(x,y,z)and surface point s=(x,y,Z_s):
d=(p-s)⋅n ̂(" " d≥0:inside" ")

	Per-cell ranked peeling (no fixed K)
Use HEA & d≥0only; within each cell, sort by d.
Cell’s 1st goes to Layer-1, 2nd to Layer-2, etc.
Each layer has ≤1 atom per cell.
	Enrichment
For radius R, composition P_s^((R))vs baseline π_s:
E_s=(P_s^((R)))/π_s 

	KDE & Posterior
PDFs f_s (r)via KDE; mixture g(r)=∑_s▒〖π_s f_s (r)〗;
P(s∣r)=(π_s f_s (r))/(g(r)),E_s (r)=(P(s∣r))/π_s .

________________________________________
Tuning Tips
	h: larger → more points per cell (deeper layers more stable); smaller → more cells but fewer ranks per cell.
	q: larger → “higher” surface, fewer candidates; typical 0.95–0.995.
	Nmin: raise to ignore sparse cells; lower to cover more area.
	Layer-1: topmost (descending d) vs closest (ascending d) changes spatial assignment.


[Readme.docx](https://github.com/user-attachments/files/22966786/Readme.docx)
