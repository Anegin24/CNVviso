# CNVviso: Visualisation solution for copy number variation data from ultra-low coverage whole-genome sequencing from preimplantation genetic testing 
![Image Alt](https://github.com/Anegin24/CNVviso/blob/4eb49cf8bfb98f97899488b7dded1e1a0a96750d/CNVviso.png)
#### Installation
You can download this git file to test it. Note: keep the folder structure intact to avoid path errors.

cd /home/user

git clone https://github.com/Anegin24/CNVviso.git

#### Shiny-server script

          CNVviso.R
          CNVdat.R
#### R packages dependencies
          library(ggplot2)
          library(DT)
          library(ggtext)
          library(tidyverse)
          library(dplyr)
          library(plotly)
          library(smoother)
          library(shiny)
          library(magick)
          library(patchwork)
          library(RSQLite)
          library(purrr)
          library(here)
#### Tutorial
          1. Compress cnvkit result file to .zip
          2. Run CNVviso.R
          3. Import zip files in "cnvkit_testadata"
          4. Click to file name and watch result
          5. Save to database
          6. You can watch imported sample without input again .zip file using CNVdat.R.
#### Citation:
An Nguyen Huu*, Nhung Trinh Thi Hong, Huong Bui Bich, Tien Vuong Quang, Dung Tien Pham and Hoi Le Thi. CNVviso: Visualisation Solution for Copy Number Variation Data from Ultra-Low Coverage Whole-Genome Sequencing in Preimplantation Genetic Testing. _Current Bioinformatics_. DOI: 10.2174/0115748936340523241113051828
          
          
          
