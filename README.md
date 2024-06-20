# VMAC Depends on explicit awareness project

This project is the repository for all data, scripts, and materials used for the manuscript entitled: "***Value-modulated attentional capture depends on explicit awareness***".

### Project Structure ###

The project is organized into the following main sections:

- **Input**: Where the raw data are stored.
- **Output**: Involves the use of R objects that are time-consuming to compute and the plot generated from **scripts**.
- **Scripts**: Includes scripts used to load and preprocess data, and execute the analysis.
- **Materials**: Contains the files with the task employed in experiment 1 and experiment 2.
- **Re_rusz_2020**: R project to re-analyze the meta-analysis by Rusz et al. (2020).


##### Scripts Folder Structure #####

The **scripts** folder is straightforward:

- *e1_analysis.R*: The script used in the **Result** section of *Experiment 1*.
- *e2_analysis.R*: The script used in the **Result and discussion** section of *Experiment 1*.
- *multiverse_analysis.R*: The script where the multiverse analysis presented as supplementary material is executed.

All the above .R files implicitly call `source(...)` to run *load_data.R* and *functions.R*. In the former, data are loaded and prepossessed; in the latter, helper functions are loaded into the workspace.

