# ddi-ps3n
PS3N: Leveraging Protein Sequence-Structure Similarity for Novel Drug-Drug Interaction Discovery

# PS3N: Leveraging Protein Sequence-Structure Similarity for Novel Drug-Drug Interaction Discovery

## Project Overview
This repository contains the code and resources for our study on predicting drug-drug interactions (DDIs) using the PS3N model. Our approach leverages protein sequence and structure information to build a similarity network for identifying novel DDIs. Due to space constraints, data files are not included in the repository but will be provided upon request or in the reviewers' comment section.

The project is organized into the following directories:

---

## Repository Structure

### 1. `Input Parsing`
This folder contains scripts for parsing and extracting relevant information from the DrugBank dataset.
- **`parse-drugbank.py`**: Reads the DrugBank XML file and extracts information such as:
  - Protein pathways
  - Drug active ingredient information
  - Drug-protein targets
  - Drug interactions
  - Drug ATC classification and pathways

Ensure the DrugBank XML data is available and adjust the file paths accordingly before running the script.

---

### 2. `Similarity Metrics Calculation`
This folder contains scripts for calculating similarity metrics for protein sequences and structures:
- **Protein Sequence and Structure Similarity:**
  - Levenshtein distance
  - Cosine similarity
  - KL divergence

- **File Information:**
  - Some scripts use generalized file names as we calculated multiple variations of these metrics.
  - Ensure data source files are available and directory paths are configured correctly to run the scripts.

---

### 3. `Model`
This folder contains the implementation and evaluation of the PS3N model:
- **`snf_implementation/`**: Implements Similarity Network Fusion (SNF) to combine multiple similarity matrices into a single fused matrix.
- **`similarity_matrix_vector.py`**: Converts the fused similarity matrix into vectors suitable for model training.
- **`k_fold_cross_validation.py`**: Performs k-fold cross-validation to evaluate the model.
- **`compare_models.py`**: Compares the PS3N model's performance with traditional machine learning models such as KNN and SVM.

---

## Running the Code

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-folder>
   ```

2. Configure paths for data files:
   - Download the required data files from the provided Google Drive link.
   - Adjust the directory paths in the scripts to match the location of the data files.

3. Run individual scripts:
   - **Input Parsing:**
     ```bash
     python InputParsing/parse-drugbank.py
     ```
   - **Similarity Metrics Calculation:**
     ```bash
     python SimilarityMetrics/<script-name>.py
     ```
   - **Model Implementation:**
     ```bash
     python Model/k_fold_cross_validation.py
     ```

---

## Data Access
Due to space limitations, the dataset is not included in this repository. A link to the data files will be provided on request. The data includes:
- DrugBank dataset (XML format)
- Processed similarity matrices for protein sequences and structures
- Feature Space
- Network for new found DDI

TO run the code you need to make sure to update the paths to the datasets in the respective scripts after downloading the files.

---

## Dependencies
Below is a list of core libraries (compatible with Python 3.7) needed to run the code in this repository:
- **numpy** (e.g., >=1.18)
- **pandas** (e.g., >=1.0)
- **scipy** (e.g., >=1.4)
- **sklearn** (e.g., >=0.22)
- **tensorflow** (e.g., >=2.2)
- **keras** (often bundled with TensorFlow 2.x)
- **networkx** (e.g., >=2.4)
- **biopython** (e.g., >=1.76)

Additionally, import statements in the code indicate dependencies such as `pairwise2` from Biopython, `model_from_json` from Keras, etc. Make sure to install them accordingly:
```bash
pip install numpy pandas scipy scikit-learn tensorflow networkx biopython
```

---

## License
This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). See the [LICENSE](LICENSE) file for more details.

## Citation
If you use this code, please cite our paper:

```
@inproceedings{islam2021detecting,
  title={Detecting drug-drug interactions using protein sequence-structure similarity networks},
  author={Islam, Saminur and Abbasi, Ahmed and Agarwal, Nitin and Zheng, Wanhong and Doretto, Gianfranco and Adjeroh, Donald A},
  booktitle={2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)},
  pages={3472--3477},
  year={2021},
  organization={IEEE}
}
```

---

## Contact
For questions or clarifications, feel free to contact:
- **Saminur Islam**: sislam8@ncsu.edu



