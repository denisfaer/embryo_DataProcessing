## Scripts used to analyze transcription factor dynamics in early mouse embryogenesis
### DF Faerberg | Posfai Lab | Princeton University

## Description
This is a custom Posfai lab repository used to develop codes to process and visualize data showing dynamics of transcription factors in pre-implantation mouse embryogenesis.
It was originally designed for use with the outputs of the computational pipeline described in (Nunley H. et al., *Development* 2024; doi: https://doi.org/10.1242/dev.202817), modified to use with the pipeline in (Kim-Yip R. & Denberg D. et al., *Current Biology* 2025; https://doi.org/10.1016/j.cub.2025.07.031), and now being adapted to processing the embryo data more broadly.

### Installation instructions
Create a new environment with  the following versions of Python and libraries:
``Python 3.10.4``
``numpy 1.21.5``
``matplotlib 3.5.3``

### Installation steps
* Install hatch
  ```
  pip install hatch
  ```

* Clone the repository
  ```
  git clone https://github.com/denisfaer/embryo_DataProcessing.git
  ```

* Navigate into the GIT folder & install the package
  ```
  pip install .
  ```

* Test the package by running StackLoad_test.ipynb

## ProcessEmbryo functions

### Load_Stack

**load_stack.py** contains core stack processing, saving and loading functions
* **frame_IDs** returns a list of cells in a given track at the given frame
* **lineage_reconstruct** recursively reconstructs a cell's lineage as a string spaced with '<'
* **lineage_transform** transforms a lineage_reconstruct string to an array
* **find_cells** finds all channel intensities for a given cell
* **lineage_express** reconstructs lineage_transform into datasets with centroid and channel(s) data
* **prune_stack** filters out short lineages in a processed stack
* **process_stack** main function that loads, reconstructs and processes the stack
* **save_stack** saves a processed stack
* **get_stack** loads a processed stack
