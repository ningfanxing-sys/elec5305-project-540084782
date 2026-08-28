# Sparse-to-Dense HRTF Reconstruction

**Cross-database effects of time alignment and regularised spherical-harmonic interpolation**

FAXNXING NING  
SID: 540084782  

[Project proposal](ELEC5305_Project_Proposal.pdf) | [GitHub repository](https://github.com/ningfanxing-sys/elec5305-project-540084782) | [Project site](https://ningfanxing-sys.github.io/elec5305-project-540084782/)

## Project summary

Head-related transfer functions (HRTFs) describe how the head and ears filter sound arriving from different directions. They are central to binaural audio, but measuring a dense set of directions for each listener is slow and requires specialist equipment.

This project asks whether a sparse set of measured directions can be reconstructed accurately enough to preserve spectral and binaural cues. It will compare nearest-neighbour interpolation with unregularised and ridge-regularised spherical-harmonic reconstruction, both with and without per-ear time alignment [1-4]. Parameters selected on HUTUBS [5] will be fixed before they are evaluated on SADIE II [6].

The project does not propose a new interpolation algorithm. It applies established signal-processing methods to test whether settings chosen in one HRTF database remain reliable in another, and to identify the directions, subjects or frequency bands where they fail. This focus is motivated by documented systematic differences among HRTF databases [7].

## Research question

How do sparse directional sampling, spherical-harmonic order, regularisation and time alignment affect HRTF reconstruction, and do settings selected on HUTUBS transfer to SADIE II without retuning?

## Study design

| Component | Plan |
|---|---|
| Data | Public HUTUBS and SADIE II HRIRs in SOFA format |
| Sampling | Sparse context sets of 25, 50, 100 and 200 directions |
| Baselines | Nearest neighbour and unregularised spherical harmonics |
| Main comparison | Ridge-regularised spherical harmonics, with and without time alignment |
| Evaluation | Log-spectral distortion, ITD error, ILD error, fitting time, query time and memory |
| Validation | Select settings on held-out HUTUBS subjects, freeze them, then test SADIE II |

The analysis will report subject-level uncertainty and directional failure maps rather than only an overall average. A short binaural A/B example will illustrate reconstructed output, but it will not be treated as evidence of perceptual transparency because no listening experiment is planned.

## Course relevance and feasibility

The work applies Fourier-domain filtering, spectral analysis, spatial sampling, correlation-based delay estimation, regularisation and binaural synthesis from ELEC5305. All core experiments will run in Python on a Mac CPU using public data. No new measurements, GPU or iPhone are required.

The implementation is deliberately bounded. It begins with a nearest-neighbour baseline and a testable spherical-harmonic pipeline before adding time alignment and cross-database evaluation. This keeps a working result available even if the more complex reconstruction route fails.

## Planned outputs

- Reproducible Python code and fixed experiment configurations
- Machine-readable evaluation results
- Scientific plots for spectral, binaural-cue and directional errors
- A compute-versus-fidelity comparison
- Written report, GitHub documentation, audio example and final video report

## Schedule

| Course weeks | Work |
|---|---|
| 1-2 | Select the topic and define the research question |
| 3-5 | Review literature and collect datasets; submit the proposal and initial project site in Week 4 |
| 6-9 | Implement preprocessing, baseline reconstruction, time alignment and initial tests |
| 10-11 | Optimise the pipeline, run the cross-database evaluation and produce figures |
| 12-13 | Finalise the working code, written report, GitHub documentation and video report |

## Selected references

Selected references are listed below. The complete numbered bibliography is provided in the project proposal.

1. Romigh, G. D. *et al.* Efficient real spherical harmonic representation of head-related transfer functions. *IEEE Journal of Selected Topics in Signal Processing* **9**, 921-930 (2015). https://doi.org/10.1109/JSTSP.2015.2421876
2. Ben-Hur, Z. *et al.* Efficient representation and sparse sampling of head-related transfer functions using phase correction based on ear alignment. *IEEE/ACM Transactions on Audio, Speech, and Language Processing* **27**, 2249-2262 (2019). https://doi.org/10.1109/TASLP.2019.2945479
3. Arend, J. M. *et al.* Assessing spherical harmonics interpolation of time-aligned head-related transfer functions. *Journal of the Audio Engineering Society* **69**, 104-117 (2021). https://doi.org/10.17743/jaes.2020.0070
4. Bau, D. *et al.* Estimation of the optimal spherical harmonics order for the interpolation of head-related transfer functions sampled on sparse irregular grids. *Frontiers in Signal Processing* **2**, 884541 (2022). https://doi.org/10.3389/frsip.2022.884541
5. Brinkmann, F. *et al.* A cross-evaluated database of measured and simulated HRTFs including 3D head meshes, anthropometric features, and headphone impulse responses. *Journal of the Audio Engineering Society* **67**, 705-718 (2019). https://doi.org/10.17743/jaes.2019.0024
6. AudioLab, University of York & Kearney, G. C. *SADIE II Database*, version 2-2 (2024). https://doi.org/10.5281/zenodo.12092466
7. Wen, Y., Zhang, Y. & Duan, Z. Mitigating cross-database differences for learning unified HRTF representation. *2023 IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA)*, 1-5 (2023). https://doi.org/10.1109/WASPAA58266.2023.10248178

## AI use acknowledgement

AI was used in this assignment, primarily for translation and polishing. AI was also used for image generation.
