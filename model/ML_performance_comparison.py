# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:49:21 2020

@author: samin
"""


from pandas import read_excel
import numpy as np
import pandas as pd
# import load_features as lf
import sys

import os
import mahotas
import matplotlib.pyplot as plt
import itertools
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier 
from sklearn.ensemble import RandomForestClassifier 
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from matplotlib import pyplot
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.model_selection import KFold
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC

def plot_confusion_matrix(cm,target_names,title='Confusion matrix',cmap=None,normalize=True):
    accuracy = np.trace(cm) / float(np.sum(cm))
    misclass = 1 - accuracy

    if cmap is None:
        cmap = plt.get_cmap('Blues')

    plt.figure(figsize=(8, 6))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()

    if target_names is not None:
        tick_marks = np.arange(len(target_names))
        plt.xticks(tick_marks, target_names, rotation=45)
        plt.yticks(tick_marks, target_names)

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]


    thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        if normalize:
            plt.text(j, i, "{:0.4f}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
        else:
            plt.text(j, i, "{:,}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")


    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label\naccuracy={:0.4f}; misclass={:0.4f}'.format(accuracy, misclass))
    plt.show()
    
target = 'I:\\Study\\Research WVU\\PDB\drug_similarity\\Generated Outputs\\feature space\\feature_space_all.csv'
df_final = pd.read_csv(target)
df_final = df_final.drop(['Unnamed: 0'],axis = 1)
df_final = lf.combine_all_features()

float_max = 3.0e+38
min_float = 1.175494351e-38
df_dd = df_final.dropna(subset=['Label'])

df_processed = pd.DataFrame()
df_processed = df_dd
df_processed = df_processed.replace([np.inf, -np.inf], float_max)
df_processed = df_processed.fillna(df_processed.mean()) 

X = df_processed.drop(['DAI1','DAI2','Label'],axis = 1)
Y = df_processed['Label']

#feature selection process neumerical input , neumerical output 
# regression feature selection 
# pearson's correlation feature selection 
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression
from sklearn.feature_selection import chi2
# define feature selection
fs = SelectKBest(score_func=f_regression, k=30)
X_selected = fs.fit_transform(X,Y)

#Univariate selection 
bestfeatures = SelectKBest(score_func=f_regression, k=30)
fifeatures = bestfeatures.fit(X,Y)
dfscores = pd.DataFrame(fifeatures.scores_)
dfcolumns = pd.DataFrame(X.columns)
#concat two dataframes for better visualization 
featureScores = pd.concat([dfcolumns,dfscores],axis=1)
featureScores.columns = ['Specs','Score']  #naming the dataframe columns
print(featureScores.nlargest(30,'Score'))  #print 10 best features

#feature importance 
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesClassifier
model = ExtraTreesClassifier()
model.fit(X,Y)
print(model.feature_importances_) #use inbuilt class feature_importances of tree based classifiers
#plot graph of feature importances for better visualization
plt.figure(figsize=(20,20))
feat_importances = pd.Series(model.feature_importances_, index=X.columns)
feat_importances.nlargest(51).plot(kind='barh',xticks =[.02,.04,.06,.08,0.1,.12,.14,.16,.18,.20,.22,.24,.26,.28,.30,.32,.34,.36,.38,.40])
plt.show()

path_corr_pear = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\Correlation\\pearson.csv'
path_corr_spear = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\Correlation\\spearson.csv'
path_corr_ken = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\Correlation\\kendall.csv'


#correlation hit map
import seaborn as sns 
from Correlation_calculation import calculate_std
from Correlation_calculation import calculate_correlation_range
from Correlation_calculation import calculate_range_counts
# pearson correlation 
corrmat = X.corr(method="pearson")
top_corr_features = corrmat.index
plt.figure(figsize=(50,50))
g=sns.heatmap(X[top_corr_features].corr(method="pearson"),annot=True,cmap="RdYlGn")
df_pear = calculate_std(corrmat)
df_pear.to_csv(path_corr_pear)
target = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\Correlation\\correlation_range.csv'
# calculate_correlation_range(corrmat,target)  
df_corr = calculate_range_counts(corrmat) 
df_corr.to_csv(target)  

 

#spearson Correlation 
corrmat = X.corr(method="spearman")
# from Correlation_calculation import calculate_correlation
# calculate_corrMat, pavl_mat = calculate_correlation(X)
top_corr_features = corrmat.index
plt.figure(figsize=(50,50))
#plot heat map
g=sns.heatmap(X[top_corr_features].corr(method="spearman"),annot=True,cmap="RdYlGn")
df_spear = calculate_std(corrmat)
df_spear.to_csv(path_corr_spear)
df_spear_std.to_csv(path_corr_spear_std)
target = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\Correlation\\correlation_range_spearman.csv'
df_corr = calculate_range_counts(corrmat) 
df_corr.to_csv(target)  

##Kendall's Tau
corrmat = X.corr(method="kendall")
top_corr_features = corrmat.index
plt.figure(figsize=(50,50))
g=sns.heatmap(X[top_corr_features].corr(method="kendall"),annot=True,cmap="RdYlGn")
#calculate avg and 
df_ken = calculate_std(corrmat)
df_ken.to_csv(path_corr_ken)
target = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\Correlation\\correlation_range_kendall.csv'
df_corr = calculate_range_counts(corrmat) 
df_corr.to_csv(target)  

#wrapper method
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from sklearn.feature_selection import RFE
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso

model = LinearRegression()
#Initializing RFE model
rfe = RFE(model, 20)
#Transforming data using RFE
X_rfe = rfe.fit_transform(X,Y)  
#Fitting the data to model
model.fit(X_rfe,Y)
print(rfe.support_)
print(rfe.ranking_)

#no of features
nof_list=np.arange(1,51)            
high_score=0
#Variable to store the optimum features
nof=0           
score_list =[]

for n in range(len(nof_list)):
    X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size = 0.3, random_state = 0)
    model = LinearRegression()
    rfe = RFE(model,nof_list[n])
    X_train_rfe = rfe.fit_transform(X_train,y_train)
    X_test_rfe = rfe.transform(X_test)
    model.fit(X_train_rfe,y_train)
    score = model.score(X_test_rfe,y_test)
    score_list.append(score)
    if(score>high_score):
        high_score = score
        nof = nof_list[n]
print("Optimum number of features: %d" %nof)
print("Score with %d features: %f" % (nof, high_score))

cols = list(X.columns)
model = LinearRegression()
#Initializing RFE model
rfe = RFE(model, 33)             
#Transforming data using RFE
X_rfe = rfe.fit_transform(X,Y)  
#Fitting the data to model
model.fit(X_rfe,Y)              
temp = pd.Series(rfe.support_,index = cols)
selected_features_rfe = temp[temp==True].index

print(selected_features_rfe)
X_new = X[selected_features_rfe]
    
from sklearn import preprocessing

le = preprocessing.LabelEncoder()
X = X.apply(le.fit_transform)
print(X.head())

test_size = 0.2
scoring = "accuracy"

trainDataGlobal, testDataGlobal,trainLabelsGlobal, testLabelsGlobal = train_test_split(X_new,Y,test_size = test_size,random_state = 42)
model = SVC(C=1000,kernel= 'rbf',random_state = 42,max_iter=10000,gamma=0.01)
model.fit(trainDataGlobal, trainLabelsGlobal)
svm_predict= model.predict(testDataGlobal)
acc_svm = model.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on svm: ",acc_svm)

cm_svm = confusion_matrix(testLabelsGlobal,svm_predict)
plot_confusion_matrix(cm_svm, normalize= False,target_names=[0,1], title = "Confusion Matrix")

    
models = []
seed=42
num_trees = 10
models.append(('SVM',SVC(random_state = seed)))
models.append(('LR',LogisticRegression(random_state = seed)))
models.append(('LDA',LinearDiscriminantAnalysis()))
models.append(('KNN',KNeighborsClassifier()))
models.append(('CART',DecisionTreeClassifier(random_state = seed)))
models.append(('RF',RandomForestClassifier(n_estimators = num_trees, random_state = seed)))
models.append(('NB', GaussianNB()))

results = []
names = []
for name, model in models: 
    kfold =KFold(n_splits = 10, random_state = seed)
    cv_results = cross_val_score(model, trainDataGlobal,trainLabelsGlobal, cv = kfold, scoring = scoring)
    results.append(cv_results)
    names.append(name)
    msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
    print(msg)
    
#boxplot algorithm comparison
fig = pyplot.figure()
fig.suptitle('machine Learning algorithm comparison')
ax = fig.add_subplot(111)    
pyplot.boxplot(results)
ax.set_xticklabels(names)
pyplot.show()



from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.preprocessing import LabelBinarizer
lb = LabelBinarizer()
vv = lb.fit_transform(Y)

labels = np.hstack((vv, 1 - vv))

trainDataGlobal, testDataGlobal,trainLabelsGlobal, testLabelsGlobal = train_test_split(X_new,Y,test_size = test_size,random_state = 42)

from sklearn.model_selection import RandomizedSearchCV
from sklearn.ensemble import RandomForestRegressor
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
max_features = ['auto', 'sqrt']
max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth.append(None)
min_samples_split = [2, 5, 10]
min_samples_leaf = [1, 2, 4]
bootstrap = [True, False]
random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}
rf = RandomForestRegressor(random_state = 42)
rf_random = RandomizedSearchCV(estimator=rf, param_distributions=random_grid,
                              n_iter = 1000, scoring='neg_mean_absolute_error', 
                              cv = 3, verbose=2, random_state=42, n_jobs=-1,
                              return_train_score=True)

rf_random.fit(trainDataGlobal, trainLabelsGlobal)
rf_random.best_params_


#random forest model for prediction
modelRF=RandomForestClassifier(n_estimators= 800, random_state=42,max_depth=80,min_samples_split=2,min_samples_leaf=1)
modelRF.fit(trainDataGlobal, trainLabelsGlobal)
rf_predict= modelRF.predict(testDataGlobal)
y_score = modelRF.predict_proba(testDataGlobal)
acc_rf = modelRF.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on RF: ",acc_rf)


#confusion metrix
cm_rf = confusion_matrix(testLabelsGlobal,rf_predict)
plot_confusion_matrix(cm_rf, normalize= False,target_names=[0,1], title = "Confusion Matrix RF")

modelKnn = KNeighborsClassifier(algorithm='ball_tree',p=2)
modelKnn.fit(trainDataGlobal, trainLabelsGlobal)
knn_predict = modelKnn.predict(testDataGlobal)
acc_knn = modelKnn.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on knn: ",acc_knn)

#confusion metrix
cm_knn = confusion_matrix(testLabelsGlobal,knn_predict)
plot_confusion_matrix(cm_knn, normalize= False,target_names=[0,1], title = "Confusion Matrix knn")


precision = dict()
recall = dict()
average_precision = dict()

y_test=[]
y_score=[]
for ii in testLabelsGlobal:
    if ii ==0:
        y_test.append(0)
    else:
        y_test.append(1)

for i in range(len(rf_predict)):
    if rf_predict[i]==0:
        y_score.append(0)
    else:
        y_score.append(1)
precision["micro"], recall["micro"], _ = precision_recall_curve(y_test,y_score)
average_precision["micro"] = average_precision_score(y_test, y_score,average="micro")
print('Average precision score, micro-averaged over all classes: {0:0.2f}'
      .format(average_precision["micro"]))

plt.step(recall['micro'], precision['micro'], where='post')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title(
    'Average precision score, micro-averaged over all classes: AP={0:0.2f}'
    .format(average_precision["micro"]))

from sklearn.metrics import roc_curve, auc

fpr, tpr, _ = roc_curve(y_test, y_score)


plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()



for i in range(2):
    precision[i], recall[i], _ = precision_recall_curve(testLabelsGlobal[:, i], rf_predict[:, i])
    plt.plot(recall[i], precision[i], lw=2, label='class {}'.format(i))

plt.xlabel("recall")
plt.ylabel("precision")
plt.legend(loc="best")
plt.title("precision vs. recall curve")
plt.show()

# roc curve
from sklearn.metrics import roc_curve
fpr = dict()
tpr = dict()

for i in range(2):
    fpr[i], tpr[i], _ = roc_curve(testLabelsGlobal[:, i],
                                  rf_predict[:, i])
    
    plt.plot(fpr[i], tpr[i], lw=2, label='class {}'.format(i))
roc_auc = auc(fpr[i], tpr[i])
plt.xlabel("false positive rate")
plt.ylabel("true positive rate")
plt.legend(loc="best")
plt.title("ROC curve")
plt.show()
# print(b.get(('DB01454','DB01571')))     

# df_cc = df_final   
# df_cc = df_cc.dropna()
# print(dd.get('DB00002')[drug_idList.index('DB00013')])
# val = pathway_dict.get('DB00001').get('DB00001')
# print(val)




