# InterroGATOR
Created: October/19/2021

Last update: December/22/2022 (v1.2)         InterroGATOR release!

Author: Leandro Balzano-Nogueira

Diabetes Institute, University of Florida (Gainesville)

## What is InterroGATOR?

InterroGATOR is a Shiny application for conducting differential expression and Spatial analyses on different cell types in different human tissues. This application allows the user to analyze cell expression profiles from CITEseq and CODEX data previously integrated (Manuscript in preparation).

InterroGATOR is the visualization and data analysis after two integrative pipelines developed in Brusko Lab, University of Florida Diabetes Institute (https://bruskolab.diabetes.ufl.edu/).

* Step 1: Single-cell CITE-seq processing and Integration.
* Step 2: CODEX/CITEseq Integration strategy.

Step 1: We used Cell Ranger package (Bryan, 2016) to demultiplex, map to the human reference genome (GRCh38), and count UMIs in the mRNA libraries. Then, we filtered out cells with more than 10% UMIs from mitochondrially encoded genes or less than 1200 mRNA UMIs in total. A second filtering step was performed following a modified procedure from Seurat v4 pipeline (Butler et al., 2018). Briefly, we calculated a statistical threshold for droplet filtering in which the idea was to retain all droplets with size larger than three times the minimum median absolute deviation (mad) and smaller than three times the maximum mad. With this we guarantee the inclusion of donor particularities in terms of cell size. Next, we followed Seurat v4 package pipeline up until the CITE-seq protein data normalization. Then, we used dsb method (Mulè et al., 2020) to to remove droplets with no cells and with more than one cell.
Once the data was pre-processed the integration of all cell donors was performed using Seurat’s scRNA-seq integration pipeline (https://satijalab.org/seurat/articles/integration_introduction.html). The idea was to create a larger cell dataset to evaluate its performance and integrate larger number of cells with CODEX datasets of these and other donors, so that we can use this data as reference for the InterroGATOR analyses.
Finally, the cell type identification was performed using Azimuth v0.4.3 (Hao, et al., 2021) using peripheral blood mononuclear cell (PBMC) as reference.  

Step 2: The strategy to integrate CODEX and CITEseq data is based on (Govek et al., 2021). Briefly, the algorithm maps the CODEX protein space to the CITE-seq protein space using a modified version of the anchor correction proposed by (Stuart et al., 2019) and implemented in Seurat v4 (Hao et al., 2021). Specifically, a set of anchors is identified using a mutual nearest neighbors approach with kanchor = 20. The algorithm finds the nearest neighbors using Euclidean distance in a common 29-dimensional space obtained by canonical correlation analysis (CCA). Then, the anchors that do not preserve the structure of the original protein space were filtered out and the cells in the CODEX dataset were aligned into the CITE-seq protein space. To be able to transfer quantities between the CITE-seq and CODEX dataset, a transfer matrix was built. 

Spatial relationship among cell populations. To assess the spatial relationship between two features we followed STvEA algorithm (Govek et al., 2021). Briefly, a k nearest neighbor using Euclidean distance of the CODEX spatial dimensions was computed. Then an adjacency matrix was calculated and used to calculate an adjacency score for each marker. The larger the score value, the larger the values of the features in adjacent cells.  The false discovery rate for multiple hypothesis testing is controlled by using the Benjamini-Hochberg procedure.

Once the data is integrated, InterroGATOR can be used to identify different features associated with a particular cell type as well as the location of these particular cells in a tissue. InterroGATOR allows to infer many aspects related to a tissue, such as:

* The identification of certain tissue structures and evaluate its behavior to determine efficiency or disease causes.
* Identify features as markers of disease progression in a particular tissue.
* Compare Differential expression profiles of different cell types to evaluate a particular tissue or even identify cell sub-types.
Amongst other utilities!


## How to install and run InterroGATOR

There are two main options to run InterroGATOR:
* Run InterroGATOR as Docker image
* Directly download and run InterroGATOR as shiny app.

## Run InterroGATOR as Docker Image
install the Docker engine.

Run InterroGATOR with the following command in terminal (Mac/Linux) or PowerShell (Win):

```
docker run --rm -p 8787:8787 leobalzano/InterroGATOR:1.0
```

Open the URL shows in the terminal (typically http://[::]:8787) in any web browser.

## Run InterroGATOR as Shiny app
* Download or clone InterroGATOR repository

```
git clone https://github.com/leobalzano/InterroGATOR.git
```

* Unzip the file in data folder
* Install all the packages listed in the script
* Run the app

## References:
Bryan, J. (2016). cellranger: Translate Spreadsheet Cell Ranges to Rows and Columns. R package version 1.1.0. 

Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411–420. https://doi.org/10.1038/nbt.4096.

Govek, K. W., Troisi, E. C., Miao, Z., Aubin, R. G., Woodhouse, S., & Camara, P. G. (2021). Single-cell transcriptomic analysis of mIHC images via antigen mapping. Science Advances, 7(10). https://doi.org/10.1126/sciadv.abc5464

Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler, A., Lee, M. J., Wilk, A. J., Darby, C., Zager, M., Hoffman, P., Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A., Stuart, T., Fleming, L. M., Yeung, B., … Satija, R. (2021). Integrated analysis of multimodal single-cell data. Cell, 184(13), 3573-3587.e29. https://doi.org/10.1016/j.cell.2021.04.048

Mulè, M. P., Martins, A. J., & Tsang, J. S. (2020). Normalizing and denoising protein expression data from droplet-based single cell profiling. BioRxiv. https://doi.org/10.1101/2020.02.24.963603

Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck, W. M., Hao, Y., Stoeckius, M., Smibert, P., & Satija, R. (2019). Comprehensive Integration of Single-Cell Data. Cell, 177(7), 1888-1902.e21. https://doi.org/10.1016/j.cell.2019.05.031.
