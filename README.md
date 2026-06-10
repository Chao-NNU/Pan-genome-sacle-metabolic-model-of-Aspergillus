# Pan-genome-scale Metabolic Models of *Aspergillus*

[![GitHub license](https://img.shields.io/github/license/Chao-NNU/Pan-genome-sacle-metabolic-model-of-Aspergillus)](https://github.com/Chao-NNU/Pan-genome-sacle-metabolic-model-of-Aspergillus)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)

📌 Introduction
This repository hosts high-quality, genome-scale metabolic models (GEMs) reconstructed for a wide range of species and strains within the genus ***Aspergillus***. By integrating pan-genomic data, this project aims to provide a comprehensive computational framework to explore the metabolic diversity, substrate utilization profiles, and secondary metabolic potentials of *Aspergillus* species from a systems biology perspective.

Currently, the repository contains **over 400 species/strain-specific metabolic models**, encompassing industrial cell factories, opportunistic pathogens, and environmentally significant strains (e.g., *Aspergillus niger*, *Aspergillus terreus*, *Aspergillus tubingensis*).

---

📂 Repository Structure
```text
├── models/                         # Core directory containing all GEMs
│   ├── Aspergillus_rambellii_SRRC1468.xml
│   ├── Aspergillus_terreus_ATCC_20542.xml
│   └── ... (400+ SBML-formatted model files)
└── README.md                       # Repository documentation

📊 Model Specifications & Validation
Scale & Diversity:
* **Scale & Diversity:** An integrated pan-GSMM framework spanning 441 *Aspergillus* strains across 77 distinct species.
* **Network Size:** The reconstructed global pan-metabolic network encompasses 4,202 metabolic reactions and 2,726 metabolites, built upon a pangenome of 24,182 gene families (where 75.9% represent variable genes driving high metabolic flexibility).
* **Format:** SBML Level 3 Version 1 (.xml), fully compatible with COBRApy and The COBRA Toolbox.
* **Predictive Reliability:** Achieved 73.28% accuracy in essential gene predictions and a 0.07 Mean Absolute Error (MAE) in growth prediction, ensuring highly robust computational simulation.

🚀 Quick Start
These models can be easily loaded and analyzed using Python (COBRApy) or MATLAB (The COBRA Toolbox) for Flux Balance Analysis (FBA).

Python Example
Ensure you have the cobra package installed (pip install cobra):
import cobra

# Load a specific Aspergillus model
model_path = "models/Aspergillus_terreus_ATCC_20542.xml"
model = cobra.io.read_sbml_model(model_path)

# Inspect model properties
print(f"Metabolites: {len(model.metabolites)}")
print(f"Reactions: {len(model.reactions)}")
print(f"Genes: {len(model.genes)}")

# Run a baseline Flux Balance Analysis (FBA)
solution = model.optimize()
print(f"Growth Rate (Objective Value): {solution.objective_value:.4f}")

📈 Key Applications & Insights
Based on constraint-based flux simulations, topology analysis, and machine learning clustering, these models enable several critical applications as demonstrated in our study:

1. **Genus-Level Metabolic Diversity & Clustering**
   Phenotypic growth simulations across 638 nutrient conditions successfully classify the 441 strains into five metabolically distinct clusters, bypassing the traditional reliance on a few limited model strains.

2. **Chassis Strain Selection via Metabolic Profiling**
   - **Cluster 2:** Characterized by robust anaerobic fermentative metabolism; computationally proven to favor the biosynthesis of fermentative organic acids like short-chain fatty acids.
   - **Cluster 4:** Characterized by enhanced aerobic and TCA cycle-associated pathways; exhibiting clear biosynthetic preferences for TCA cycle-related organic acids including citrate, malate, and pyruvate.

3. **FSEOF/FVSEOF-Driven Metabolic Engineering Target Discovery**
   Identification of 1,310 feasible metabolic engineering targets capable of enhancing organic acid production while maintaining cellular biomass growth. This includes both experimentally validated reactions and novel transport- and cofactor-associated engineering targets, bridging the gap between genus-level diversity and rational cell-factory design.

🤝 Citation & Contact
If you use these models or the dataset in your research, please cite our work:

Chao, et al. (2026). Machine learning–guided pan-genome-scale metabolic modeling reveals genus-wide metabolic diversity and biosynthetic potential in Aspergillus. (Manuscript in preparation).

Author: Chao (Chao-NNU)

Institution: Nanjing Normal University (NNU)

Email: chaoye09@njnu.edu.cn
