# Project Summary:

This study introduces a novel ensemble classification method that leverages the predictions of existing wear-time detection algorithms to improve the accuracy of classifying wear time versus physical activity, sleep, and sedentary time. The Choi [1], Nhanes [2], and hidden Markov model [3] algorithms were utilized, and their predictions were combined using two unsupervised ensemble methods [4,5]. Detecting intervals of accelerometer wear time in field-collected data is a challenging problem and may lead to data that does not reflect actual physical activity. By accurately identifying non-wear periods, this method can reduce biases in data analysis and lead to more reliable health outcome assessments.

# Dataset:

Activity signals from three datasets were analyzed: the Sleep Study with 15 children (ground truth labels available), the Stanford GOALS trial with 268 children (no ground truth labels) [6], and a simulation dataset. 

# Repository Structure:

The repo consists of five files listed below:

1. WeartimeDetection: This R file contains the necessary code for weartime detection using the Choi, Nhanes and HMM algorithms on the specified accelerometer data.
2. RunWeartimeDetection: This R file contains the necessary code for running wear time detection algorithms across all subjects.
3. SUMMA_BI: This R file contains the necessary code for wear time detection using two unsupervised ensemble methods, SUMMA and BI, utilizing the predictions of existing algorithms such as the Choi, Nhanes, and HMM algorithms.
4. SimulateDataset: This script generates data that simulates 1-day of accelerometer data, where the weartime is from the GOALS study and the nonwear is from the Sleep study.
   Each simulated sample represents 1-day accelerometry readings for a simulated subject.
5. confusion_matrix: This script generates confusion matrix of different predictors

# Citations:

 1. Choi, L., Liu, Z., Matthews, C.E., & Buchowski, M.S. (2011). Validation of accelerometer wear and non-wear time classification algorithm. Medicine & Science in Sports & Exercise, 43(2), 357–364.  
 2. Troiano RP, Berrigan D, Dodd KW, Mâsse LC, Tilert T, McDowell M. Physical activity in the United States measured by accelerometer. Med Sci Sports Exerc 2008;40(1):181-8. PubMed. 
 3. Randhawa, S., Sharma, M., Fiterau, M., Banda, J. A., Haydel, F., Kapphahn, K., Matheson, D., Moore, H., IV, Ball, R. L., Kushida, C., Delp, S., Wall, D. P., Robinson, T., & Desai, M. (2023). Statistical Learning Methods to Identify Nonwear Periods From Accelerometer Data. Journal for the Measurement of Physical Behaviour, 6(2), 124-133. Retrieved May 4, 2024, from https://doi.org/10.1123/jmpb.2022-0034.
 5. Emmanouil Antonios Platanios, Avrim Blum, and Tom M Mitchell. Estimating Accuracy from Unlabeled Data. In Conference on Uncertainty in Artificial Intelligence, pages 1-10, 2014. 
 6. Mehmet Eren Ahsen, Robert M Vogel, Gustavo A Stolovitzky, Unsupervised Evaluation and Weighted Aggregation of Ranked Classification Predictions. Journal of Machine Learning Research 20 (2019) 1-40.
 7. Robinson TN, Matheson DM, Desai M, Wilson DM, Weintraub DL, Haskell WL, et al. Family, community and clinic collaboration to treat overweight and obese children: Stanford GOALS—a randomized controlled trial of a three-year, multi-component, multi-level, multi-setting intervention. Contemp Clin Trials. 2013;36:421–435. 10.1016/j.cct.2013.09.001

![Performance_metrics.png)
