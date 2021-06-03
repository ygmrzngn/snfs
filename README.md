# snfs
Biomarker Detection with Social Network Feature Selection in R

Dimension reduction is a crucial part of a classification problem when the data is high dimensional. Particularly in cancer research, discovering potential biomarker genes is important due to the fact that using reduced data could improve diagnostic accuracy. Therefore, researchers continuously try to explore more efficient ways for reducing the large number of features/genes to a small but informative subset before the classification task.

In this study we apply the Social Network Feature Selection method as a hybrid feature selection approach in single environment, R to improve Support Vector Machine classifier’s performance. In addition, we evaluate and compare the performances of several combinations used in the different steps of the method with a simulation experiment. As a result, the method slightly improves the classifier’s performance compare to using whole feature set in all the cases we investigated.   

The function apply the SNFS method and summarize the results. Users can apply different combinations of the algorithms in SNFS to determine which combination performs better in binary classification problem. In addition, there are two functions which can be used for converting exprsSet object to ExpressionSet object and spliting ExpressionSet object to test/train with given train ratio before the analysis. 
