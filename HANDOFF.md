# ISAC / AFDM project: handoff for the next session

Snapshot written 2026-10-08. Read this first. Owner: Nguyen Quang Khai (PTIT). Working language with the user: Vietnamese.

## 1. What this project is

A university report on **AFDM for ISAC** (Integrated Sensing and Communication), with MATLAB simulations of AFDM / OTFS / OFDM (/ OCDM) in doubly-dispersive (multipath + high-Doppler) channels.

This folder started from the CC0 GitHub repo "doubly-dispersive-channel-simulation" (root `*.m` files, see `LICENSE`) and was extended with:

| Folder | Content | Status |
|---|---|---|
| root `*.m` | Original waveform/channel/sensing code (`waveformMod`, `genDDChannel`, `simBERSER`, `simSensing`, `simDiversityML`, ...) plus own additions (`estimateDiversitySlope`, `applyChannelDamage`, `simDamageRobustness`, velocity sweeps) | Older experiments, **removed from the current report**, kept for reference |
| `rm_isac_demo/` | Older ISAC range/velocity NMSE and CAF experiments | Reference only |
| `lisat_check/` | Python waveform sanity check | Reference only |
| `afdm_paper_check/` | **Current work**: own MATLAB implementation of the simulation scenarios of Bemani et al. 2023 (7 figures) | Active |
| `report_build/` | docx-js build scripts of the report (`build_theory.js`) and figure copies | Active |

## 2. Source paper (not included in the repo: copyright)

[2] A. Bemani, N. Ksairi, M. Kountouris, "Affine Frequency Division Multiplexing for Next Generation Wireless Communications", IEEE Trans. Wireless Commun., 2023. arXiv:2204.12798 (https://arxiv.org/abs/2204.12798).
[1] Same authors, "AFDM: A Full Diversity Next Generation Waveform for High Mobility Communications", arXiv:2104.11331.

To re-check numbers, download the PDF and (optionally) put it at `afdm_paper_check/afdm_paper.pdf` (git-ignored).

`afdm_paper_check/Gen_channel_mtx.m` is **third-party code** (R. Zhang, UCAS, BSD-style licence, header must be kept; the reference named in its header must be cited if used in a publication).

## 3. The 7 simulation experiments (`afdm_paper_check/`)

Own code written from the formulas, not copied from the paper. Parameters follow [2] Section VII (fc = 4 GHz everywhere; alpha_max = 1 is about 405 km/h, alpha_max = 2 is about 540 km/h).

| Report figure | Script | Results `.mat` | Replot script | What it shows |
|---|---|---|---|---|
| Hinh 2 | `fig8a_afdm_diversity.m` | `fig8a_results.mat` | `replot_fig8a.m` | AFDM ML, N=16, P=2/3/4, integer Doppler: BER slope grows with P (diversity) |
| Hinh 3 | `fig8b_four_waveforms.m` | `fig8b_results.mat` | `replot_fig8b.m` | ML, OFDM/OCDM/OTFS/AFDM, N=16, P=3, Jakes Doppler |
| Hinh 4 | `fig9_mmse_comparison.m` | `fig9a_results.mat` | `replot_fig9a.m` | MMSE, 4 waveforms, N=256 |
| Hinh 5 | `fig10_detectors.m` (+ `fig10a_fix_tail.m`, `fig10b_fix_lowc.m` patches) | `fig10_results.mat` | `replot_fig10.m` | MMSE vs MRC-DFE vs banded low-complexity MMSE; integer (left) and fractional xi_nu=1 (right) Doppler |
| Hinh 6 | `fig11_channel_estimation.m` | `fig11_results.mat` | `replot_fig11.m` | Embedded-pilot estimation, SNRp=25/30/35 dB vs Ideal, integer Doppler |
| Hinh 7 | `fig12a_channel_estimation_frac.m` | `fig12a_results.mat` | `replot_fig12a.m` | Same, fractional Doppler xi_nu=1, SNRp=30/35/40 dB |
| Hinh 8 | `fig12b_xinu_effect.m` | `fig12b_results.mat` | `replot_fig12b.m` | xi_nu = 1/4/8 at SNRp=40 dB |

Support: `AFDM_mod.m`, `AFDM_demod.m`, `afdm_pilot_estimate.m` (successive-cancellation matching pursuit pilot estimator, own implementation), `Gen_channel_mtx.m` (third party), `diag_*.m` (diagnostics used while debugging), `*_run*.log` (run logs).

Workflow: run the `figXX_*.m` script (slow, writes `.mat`), then `replot_*.m` (fast, loads `.mat`, writes `figXX_reproduction.png` **without titles**, LineWidth 2). The replot scripts contain **hard-coded absolute paths** `D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\...`: edit them on another machine. Copy the PNGs to `report_build/afdm_repro/` for the report.

### Decisions and pitfalls found (do not undo)

1. **Doppler model (fig10)**: continuous Jakes sampling `doppler_norm = alpha_max*cos(theta)`, then `round()` for the integer case. Forcing distinct integer Dopplers (`randperm`) was wrong. Natural Doppler collisions are exactly what produces the OFDM error floor (about 1.25e-4 at 30 dB, close to the paper's value) while AFDM stays robust.
2. **Banded low-complexity MMSE**: `Lband = 3*(2*xi_nu+1)*P` for the fractional case (empirical factor 3, because the simulated Heff is the full circulant matrix with no null-guard truncation). Without it AFDM-lowC wrongly looks worse than OFDM-lowC. The factor lives in `fig10_detectors.m`; the fractional lowC curve was patched via `fig10b_fix_lowc.m`.
3. **Monte Carlo noise**: numMonte is 300 to 30000 (paper: 1e6). High-SNR tail of Fig 10a used numMonte=30000 (`fig10a_fix_tail.m`) to avoid fake "upturns" caused by 4 to 8 counted errors. Never run two scripts that load/modify/save the same `.mat` in parallel.
4. Zero BER points are floored (`floorZeros` in replot scripts) so the log plot works.

## 4. The report

- Source: `report_build/build_theory.js` (docx-js, `npm install` in `report_build/` first; `package.json` pins `docx`). Despite its name it generates the current **ISAC/AFDM report** `LyThuyet_AFDM_OTFS_OFDM.docx` (17 pages, delivered to the user as `E:\Baocao_AFDM_OTFS_OFDM_v18.docx`).
- Structure: Part I theory (sections 1-7: M/A matrices, diversity theorems, etc.); Part II simulation: 8 intro, 9 parameters (Table 1 incl. velocity column) / detector and estimation formulas / complexity, 10 effective channel matrix + 10.1 banded structure (formulas L and Q), 11-17 the 7 figures (Hinh 2-8), 18 discussion, 19 limitations, 20 conclusion, Appendix A (file list).
- Build: `node build_theory.js`; before rebuilding kill stale Word (`Get-Process WINWORD | Stop-Process -Force`) or you get `EBUSY`. PDF via Word COM (`SaveAs`, format 17); page images via PyMuPDF (`render_pages.py`).
- Verified visually page by page in v18 (no overflow, figures titleless).
- The user's **full merged report** lives outside this repo: `D:\BAO_CAO_ISAC_HOAN_CHINH.docx` (chapters 1-2 from `BAOCAO_ISAC.docx`, chapter 3 = this v18, PTIT cover from `E:\bia.docx`). Built by a Python raw-XML pipeline in an older session scratchpad (not saved here). Cover placeholders (giang vien, nhom, lop, thanh vien/MSV) are still yellow-highlighted and must be filled by the user. Known corrections vs v18 found during that merge: M for AFDM is IDAFT (A^H); sign of W^nu is minus; c2 in the code is rational (1/N^2, 1/(2N^2)); table of lmax / xi_nu corrected against the MATLAB; `Gen_channel_mtx.m` must cite Bemani TWC. **`build_theory.js` has NOT yet received these corrections**: apply them if v18 is reused.

### User's writing rules (mandatory, from explicit feedback)

- It is the **user's own report**; the paper is only a reference for method and parameters. Never write "khop voi bai bao goc" / "xac nhan bai bao" / "kiem chung doc lap". Phrase as own investigation, cite [2] for formulas and parameters.
- No arbitrary ALL-CAPS emphasis; professional tone.
- Figures: **no title above the chart**; explanations must be mechanistic and not sparse (Quan sat / Giai thich / Ket luan, lead with the conclusion, ISAC-oriented framing: sensing/communication trade-off, multi-target, pilot cost).
- Vietnamese academic style; keep each figure analysis short but explained.

## 5. Open items / what to do next

1. **Add the user's own experiment** (the user asked whether results only reproduce the paper; the honest answer is yes, parameters follow [2]). Proposed options, not yet chosen: velocity sweep from low to above 540 km/h (when does AFDM degrade as alpha_max exceeds the design value; Theorem 1 premise), larger P / lmax, an ISAC-like multi-target scenario with close Dopplers. Ask the user which one.
2. Optionally fold the three explanations already given in chat into the report for readability: why BER drops with P (diversity, probability that all P independent paths fade ~ SNR^-P), what MMSE / MRC-DFE / banded MMSE do and what SNRp is (pilot power over noise, `SNRp = |x_pilot|^2/N0`, controls matching-pursuit accuracy), why the three detectors coincide for integer Doppler (Heff exactly banded with L=P non-zero entries, so truncation loses nothing; for fractional Doppler L=(2 xi_nu+1)P only approximates the leakage).
3. Propagate the corrections listed in section 4 into `build_theory.js` and the merged report; update the user's file and redeliver.
4. Fill the cover placeholders (user task).

## 6. Environment notes

Windows 10, MATLAB, Node (docx-js), Python with PyMuPDF (`fitz`), MS Word (COM automation). In Git Bash `zip` is missing: use Python `zipfile`. Pandoc reports 0 headings for the docx-js output; that is cosmetic (Word renders the headings correctly), verify through the PDF instead.
