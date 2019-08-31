# Recent-Computational-Methods-for-White-Blood-Cell-Nuclei-Segmentation-A-Comparative-Study
This repository contains implementations from methods evaluated in Recent Computational Methods for White Blood Cell Nuclei Segmentation A Comparative Study.

Abstract
Background and objective: Leukaemia is a disease found worldwide; it is a type of cancer that originates in
the bone marrow and is characterised by an abnormal proliferation of white blood cells (leukocytes). In order
to correctly identify this abnormality, haematologists examine blood smears from patients. A diagnosis obtained
by this method may be influenced by factors such as the experience and level of fatigue of the haematologist, 
resulting in non-standard reports and even errors. In the literature, several methods have been proposed that 
involve algorithms to diagnose this disease. However, no reviews or surveys have been conducted. This paper 
therefore presents an empirical investigation of computational methods focusing on the segmentation of leukocytes.

Methods: In our study, 15 segmentation methods were evaluated using five public image databases: ALL-IDB2, BloodSeg
, Leukocytes, JTSC Database and CellaVision. Following the standard methodology for literature evaluation, we conducted
a pixel-level segmentation evaluation by comparing the segmented image with its corresponding ground truth. In order
to identify the strengths and weaknesses of these methods, we performed an evaluation using six evaluation metrics:
accuracy, specificity, precision, recall, kappa, Dice, and true positive rate.

Results: The segmentation algorithms performed significantly differently for different image databases, and for each
database, a different algorithm achieved the best results. Moreover, the two best methods achieved average accuracy
values higher than 97%, with an excellent kappa index. Also, the average Dice index indicated that the similarity between
the segmented leukocyte and its ground truth was higher than 0.85 for these two methods This result confirms the high
level of similarity between these images but does not guarantee that a method has segmented all leukocyte nuclei. 
We also found that the method that performed best segmented only 58.44% of all leukocytes.

Conclusions: Of the techniques used to segment leukocytes, we note that clustering algorithms, the Otsu threshold,
simple arithmetic operations and region growing are the approaches most widely used for this purpose. However, these
computational methods have not yet overcome all the challenges posed by this problem.

Full paper at https://www.sciencedirect.com/science/article/pii/S0169260718311064

Cite with bibtex.

@article{ANDRADE20191,
title = "Recent computational methods for white blood cell nuclei segmentation: A comparative study",
journal = "Computer Methods and Programs in Biomedicine",
volume = "173",
pages = "1 - 14",
year = "2019",
issn = "0169-2607",
doi = "https://doi.org/10.1016/j.cmpb.2019.03.001",
url = "http://www.sciencedirect.com/science/article/pii/S0169260718311064",
author = "Alan R. Andrade and Luis H.S. Vogado and Rodrigo de M.S. Veras and Romuere R.V. Silva and Flávio H.D. Araujo and Fátima N.S. Medeiros",
}
