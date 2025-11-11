## Classification Datasets for Various Domains

## Real-World Datasets

Here are 6 datasets for classification tasks from the fields of medical diagnostics, rare-event detection, and industrial quality control, along with direct links.

### Medical Diagnostics

1.  **Heart Failure Prediction Dataset**
    *   **Description:** This dataset contains 11 clinical features that can be used to predict whether a patient has heart disease. It is a classic binary classification problem.
    *   **Target Variable:** `HeartDisease` (1 for heart disease, 0 for normal).
    *   **Link:** [https://www.kaggle.com/datasets/fedesoriano/heart-failure-prediction](https://www.kaggle.com/datasets/fedesoriano/heart-failure-prediction)

2.  **Breast Cancer Wisconsin (Diagnostic) Data Set**
    *   **Description:** The features in this dataset are computed from a digitized image of a fine needle aspirate (FNA) of a breast mass. They describe characteristics of the cell nuclei present in the image. The task is to classify tumors as malignant or benign.
    *   **Target Variable:** `diagnosis` (M = malignant, B = benign).
    *   **Link:** [https://archive.ics.uci.edu/dataset/17/breast-cancer-wisconsin-diagnostic](https://archive.ics.uci.edu/dataset/17/breast-cancer-wisconsin-diagnostic)

3.  **EEG Brainwave Dataset: Feeling Emotions**
    *   **Description:** This dataset contains EEG brainwave data that has been processed using statistical feature extraction. The data was collected from two individuals (one male and one female) for three minutes per emotional state: positive, neutral, and negative. The dataset is intended for use in tasks such as sentiment classification based on EEG data.
    *   **Target Variable:** `emotion` (Positive, Neutral, Negative).
    *   **Link:** [https://www.kaggle.com/datasets/birdy654/eeg-brainwave-dataset-feeling-emotions](https://www.kaggle.com/datasets/birdy654/eeg-brainwave-dataset-feeling-emotions)

### Rare-Event Detection

4.  **Credit Card Fraud Detection**
    *   **Description:** This dataset contains credit card transactions made over two days by European cardholders. The challenge is to identify fraudulent transactions, which constitute a very small fraction of all transactions (a minority class).
    *   **Target Variable:** `Class` (1 for fraud, 0 for a legitimate transaction).
    *   **Link:** [https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud)

### Industrial Quality Control

5.  **Steel Plates Faults Dataset**
    *   **Description:** This dataset is from experiments on steel plates to identify different types of faults. The goal is to classify the fault type based on 27 different attributes. This is a multi-class classification problem.
    *   **Target Variable:** 7 different types of faults (e.g., `Pastry`, `Z_Scratch`, `K_Scatch`).
    *   **Link:** [https://archive.ics.uci.edu/dataset/198/steel-plates-faults](https://archive.ics.uci.edu/dataset/198/steel-plates-faults)

6.  **Sensorless Drive Diagnosis Data Set**
    *   **Description:** Data was collected from a sensorless drive to diagnose its condition. The task is to determine if a fault is present, based on 48 features from the current signal.
    *   **Target Variable:** The drive's condition (11 different fault classes).
    *   **Link:** [https://archive.ics.uci.edu/dataset/325/sensorless+drive+diagnosis](https://archive.ics.uci.edu/dataset/325/sensorless+drive+diagnosis)

## Synthetic Datasets

### Synthetic Data Generation for Classifier Evaluation

*   **Description:** This is not a fixed dataset, but a detailed process for generating synthetic data to test and compare classification models, particularly discriminant analysis methods. The generation process allows for fine-grained control over the data's statistical properties.

    **The data generation process consists of the following steps:**
    1.  **Sampling the S parameter for the Wishart distribution:** A positive-definite matrix `S` is generated to serve as the scale matrix for the Wishart distribution.
    2.  **Sampling covariance matrices:** Covariance matrices are drawn from the Wishart distribution `W_p(S, df)`.
    3.  **Permutations of covariance matrices:** The sampled covariance matrices are projected onto specific permutations, with the projection method depending on the chosen scenario.
    4.  **Sampling means:** The mean vectors for each class are sampled as random points from a p-dimensional unit cube.
    5.  **Modifying class separability (psi parameter):** A parameter `psi` is used to multiply the covariance matrix, controlling the degree of linear separability between classes (a larger `psi` increases the overlap).
    6.  **Data generation:** The final data for each class is sampled from a multivariate normal distribution using the means and covariance matrices from the previous steps.

    **Five scenarios for data generation are considered:**
    *   **Scenario 1: gipsLDA model:** Same permutations and covariance matrices, different means. Analogous to classical LDA.
    *   **Scenario 2: gipsMultQDA model:** Same permutations, different covariance matrices and means.
    *   **Scenario 3: gipsQDA model:** Different permutations, covariance matrices, and means. Analogous to QDA.
    *   **Scenario 4: Classic LDA model:** Same covariance matrix, different means. Serves as a baseline.
    *   **Scenario 5: Classic QDA model:** Different covariance matrices, different means. Serves as a baseline.
*   **Target Variable:** A class label generated according to the specified scenario parameters.
*   This is a data generation methodology, not a downloadable dataset. The process allows for creating custom datasets for research and benchmarking purposes.