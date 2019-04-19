# ShinyDEGsAnalysis: Shiny App for Interactive Visualization of Differentially Expressed Gene Results Generated by systemPipeR
Development repository  

ShinyDEGsAnalysis provides interactive visualizations to enable wet-lab scientists easier access to the DEG results generated by systemPipeR. This enables quicker interpretation of the results needed to further their research.

## Installation

Clone from github.  
At your terminal prompt type: git clone https://github.com/gdmosher/ShinyDEGsAnalysis.git  

Requires Shiny, ShinyDashboard  
Requires systemPipeR, DESeq2, ape, edgeR, pheatmap  


## Usage

Easiest to run it from inside RStudio.  Use the project file ShinyDEGsAnalysis.Rproj to open the project.

In RStudio open the file runApp_ShinyDEGsAnalysis.R for editing, then click [Run App] button to run it.  
Demo data is included.  
Click on two tabs in Shiny Dashboard App to see DEG visualizations.  Move sliders to adjust parameters for plots.

Click to upload a new counts file. edgeR and biomaRt will run automatically.  

Images can be downloaded by right clicking on them. They are also saved to disk as Png files in the folder sda/results.  