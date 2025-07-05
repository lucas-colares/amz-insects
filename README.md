# 🦟 Life History Drives Insect Responses to Forest Loss

**Repository for:** *Life history induces markedly divergent insect responses to habitat loss*
**Published in:** *Journal of Animal Ecology (JAE-2024-00872)*
**Lead author:** Lucas Ferreira Colares
**DOI:** [10.5281/zenodo.15238078](https://doi.org/10.5281/zenodo.15238078)

---

## 🌎 Study Area

This research was conducted in the **Balbina Hydroelectric Reservoir**, Central Amazonia (Brazil), which forms the **world’s largest man-made forest archipelago**. In 1986, the dam submerged 251,216 ha of continuous old-growth forest, creating >3,500 forest islands of varying size and isolation. This setting provides a natural experiment to assess the effects of habitat loss and fragmentation on tropical biodiversity.

* 🌍 **Coordinates:** \~1°55′S, 59°30′W
* 🏝️ **Reservoir size:** \~360,350 ha
* 🌳 **Vegetation:** Dense terra firme rainforest
* 📏 **Island size range:** 0.2 ha to 4,878 ha
* 🐜 **Protection status:** Surrounded by the Uatumã Biological Reserve (940,000 ha)

---

## 🧪 Sampling Design

* **Sampling period:** October–November 2021
* **Traps used:** 236 yellow sticky traps (double-sided, 20 × 15 cm = 600 cm²)
* **Trap duration:** 24 hours
* **Deployment locations:**

  * 17 forest islands (edge-to-core transects)
  * 3 continuous mainland forest sites
  * 8 floating aquatic traps per island in the reservoir matrix (spaced logarithmically: 10–4000 m)

**Each trap was photographed on both sides**, resulting in 472 total high-resolution images.

---

## 🖼️ Image Processing Workflow

1. **Camera used:** Panasonic DMC-FZ60 (fixed position and settings)
2. **Image segmentation (R):**

   * Grayscale conversion
   * Adaptive thresholding (to handle uneven lighting)
   * Morphological dilation
   * Cluster analysis (Gower’s distance) to isolate individual insects
3. **Output:** 14,090 sub-images generated for object detection
4. **Manual annotation:** 999 sub-images annotated for model training

Packages used: `imager`, `EBImage`, `stats`, `cluster`, `tidyverse`

---

## 🤖 Deep Learning for Detection

* **Model:** YOLOv8 (Ultralytics implementation)

* **Training dataset:** 999 manually annotated insect sub-images

* **Training strategy:**

  * 5-fold cross-validation
  * 250 epochs, batch size = 50
  * Image size = 640×640

* **Inference:**

  * Full-size images = 4608×3465 pixels
  * Detection performed with **SAHI (Sliced Aided Hyper Inference)** to improve performance on small insects

* **Metrics assessed per fold:**

  * Precision, Recall, F1, mAP\@0.5, optimal threshold via simulation

* **Post-processing:**

  * Inferences filtered by IoU > 0.5 and max confidence score
  * Bounding box area used as a **proxy for body size**

Packages used: `ultralytics`, `sahi`, `opencv`, `numpy`, `pandas`, `matplotlib`

🖥️ Training performed on **NVIDIA A100 GPU**

---

## 🧬 Functional Classification of Insects

Insects were grouped by **life history** strategy:

* **Aquatic insects** (with at least one aquatic larval stage):

  * Ephemeroptera (mayflies)
  * Nematocera (mosquitoes)
  * Trichoptera (caddisflies)
* **Terrestrial insects**:

  * Brachycera (flies)
  * Coleoptera (beetles)
  * Hemiptera (cicadas)
  * Hymenoptera (bees and wasps)

Insects with <30 annotated individuals or poor classification metrics (e.g., Orthoptera, Lepidoptera) were excluded.

---

## 🧮 Ecological & Statistical Analysis

All statistical analyses were performed in **R**:

### 📊 Diversity and Composition

* **Rarefaction analysis (individual-based):** Partitioned into:

  * N-component (total abundance)
  * SAD-component (species pool size)
  * Using methods from Engel et al. (2022)

* **Community composition:**

  * **PCoA** (Principal Coordinates Analysis) on:

    * Jaccard distance (presence/absence)
    * Bray–Curtis distance (abundance)
  * **βC diversity index** (coverage-based rarefaction with null models)

### 📈 Modeling

* **Generalized Additive Models (GAMs):**

  * Modeled the effect of **forest cover (%)** on:

    * Diversity (N, SAD)
    * Composition (PCoA scores)
    * Taxon-specific abundance and body size
  * **Spatial smoothing** included: `s(x, y)`
  * **Distribution types**:

    * Gaussian (diversity, PCoA)
    * Poisson (abundance)
    * Quasi-Poisson (body size)
  * **Model selection:** Based on lowest AIC across 20 `k` values

Packages used: `mgcv`, `vegan`, `betapart`, `iNEXT`, `ggplot2`, `picante`

---

## 📦 External Data & Model Resources

| Resource                              | DOI / URL                                                |
| ------------------------------------- | -------------------------------------------------------- |
| 📸 Full-resolution trap images (472)  | [Figshare](https://doi.org/10.6084/m9.figshare.23823591) |
| 🧪 5-fold YOLOv8 training dataset     | [Figshare](https://doi.org/10.6084/m9.figshare.28688198) |
| 🤖 YOLOv8 trained models (5 folds)    | [Figshare](https://doi.org/10.6084/m9.figshare.28820993) |
| 💾 Code + processed data for analysis | [Zenodo](https://doi.org/10.5281/zenodo.15238078)        |

---

## 📚 Citation

If using this repository, please cite:

> **Colares, L.F.**, et al. (2025). *Life history induces markedly divergent insect responses to habitat loss*. Journal of Animal Ecology. [https://doi.org/10.5281/zenodo.15238078](https://doi.org/10.5281/zenodo.15238078)

Also cite the **Figshare DOIs** if using the models or datasets directly.

---

## 📬 Contact

**Lucas Ferreira Colares**
Postdoctoral Researcher – UFPA | Project IARAA
📧 [lucasfcolares@gmail.com](mailto:lucasfcolares@gmail.com)
[ResearchGate](https://www.researchgate.net/profile/Lucas-Colares)
