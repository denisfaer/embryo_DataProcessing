## Scripts used to analyze and visualize transcription factor dynamics in early mouse embryogenesis
### DF Faerberg | Posfai Lab | Princeton University

## Description
This is a custom Posfai Lab repository used to process and visualize data on dynamics of transcription factor expression in pre-implantation mouse embryogenesis.
It is explicitly focused on processing the outputs of the computational pipeline described in (Nunley H. et al., Development 2024; doi: https://doi.org/10.1242/dev.202817)

## Installation instructions

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

### Required libraries
The following libraries are used:
* numpy
* csv
* os
* pathlib
* pickle
* matplotib

### Details on src/ProcessEmbryo
**load_stack.py** contains core stack processing, saving and loading functions
