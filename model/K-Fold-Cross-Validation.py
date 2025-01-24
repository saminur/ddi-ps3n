# -*- codini: utf-8 -*-
"""
Created on Sun Nov  1 20:20:50 2020

@author: samin
"""

from pandas import read_excel
import numpy as np
import pandas as pd
import csv 
# alignment_path = 'I:\\Study\\Research WVU\\PDB\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS3\\pdb_KL_weighted_simi.csv'
# alignment_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\DeepDDI dataset\\SNF matrix\\merged_feature_sequence.csv'
# alignment_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS4_V2\\test\\'
alignment_path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\df_seq_merge.csv'
# f_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\NDD Dataset\\ds1_pdb_ndd\\merged_pdb_vector_ds1.csv'
df_alignment = pd.read_csv(alignment_path)
df_alignment = df_alignment.set_index('Unnamed: 0')
# df_alignment= df_alignment.fillna(1)
# df2 = (df_alignment - df_alignment.values.min()) / (df_alignment.values.max()-df_alignment.values.min())
# df_alignment = df2
f_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS1_V2\\merged_feature.csv'

# f_path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\merged_vector_new.csv'
f_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS3_V2\\merged_feature.csv'
drugList = list(df_alignment.keys())
with open(f_path, 'w',newline='') as f:
    writer = csv.writer(f)
    for ii in range(len(drugList)):
        for jj in range(ii+1,len(drugList)):
            firstList = df_alignment[drugList[ii]]
            scndList = df_alignment[drugList[jj]]
            output = list(firstList - scndList)
            row = [drugList[ii],drugList[jj]] + output
            writer.writerow(row)

# ddi_path = 'I:\\Study\\Research WVU\\PDB\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS6\\ddi_indi.csv'
# ddi_al = pd.read_csv(ddi_path)
# ddi_al = ddi_al.set_index('Unnamed: 0')
# path = 'I:\\Study\\Research WVU\\PDB\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS5\\ddi_kmers.csv'
# path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS1_V2\\ddi.csv'
# path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\ddi_seq_cardio.csv'
path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS3_V2\\ddi_pdb.csv'
# path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\ddi_merged_new.csv'
ddi_al = pd.read_csv(path)
ddi_al = ddi_al.set_index('Unnamed: 0')
path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\ddi_seq_cardio_vector.csv'
drugList = list(ddi_al.keys())
with open(path, 'w',newline='') as f:
    writer = csv.writer(f)
    for ii in range(len(drugList)):
        for jj in range(ii+1,len(drugList)):
            row = [drugList[ii],drugList[jj]] + [ddi_al[drugList[ii]][drugList[jj]]]
            writer.writerow(row)

# from scipy.stats import zscore
# ddi_al = pd.read_csv(alignment_path)
# ddi_al = ddi_al.set_index('Unnamed: 0')
# ddi_al= ddi_al.fillna(ddi_al.mean(axis=0))
# ddi_al  =ddi_al.T
# sd = np.std(ddi_al)
# mean= np.mean(ddi_al)
# numerator = ddi_al - mean
# z_score = numerator/sd
# z_norm = z_score.T
# f_path = 'I:\\Study\\Research WVU\\PDB\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS5\\alignment_score_median_simi_zscore.csv'
# z_norm.to_csv(f_path)

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
    

feature_space = pd.read_csv(f_path,header=None)
FS = pd.DataFrame(feature_space)

FS = FS.replace(-np.inf, np.nan)
FS = FS.replace(np.inf, np.nan)
FS = FS.fillna(FS.mean())
X = FS.drop([0,1],axis=1)
labels = pd.read_csv(path,header=None)
labels = labels.drop([0,1],axis=1)
Y = labels[2]

from sklearn.preprocessing import LabelBinarizer
lb = LabelBinarizer()
vv = lb.fit_transform(Y)

Y_labels = np.hstack((vv, 1 - vv))

trainDataGlobal, testDataGlobal,trainLabelsGlobal, testLabelsGlobal = train_test_split(X,Y_labels,test_size = 0.2,random_state = 42)

from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score

model = SVC(C =10,random_state = 42,max_iter=10000,gamma=0.05)
model.fit(trainDataGlobal, trainLabelsGlobal)
svm_predict= model.predict(testDataGlobal)
acc_svm = model.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on svm: ",acc_svm)

# cm_svm = confusion_matrix(testLabelsGlobal,svm_predict)
# plot_confusion_matrix(cm_svm, normalize= False,target_names=[0,1], title = "Confusion Matrix")


precision = precision_score(testLabelsGlobal, svm_predict)
print(precision)
recall = recall_score(testLabelsGlobal, svm_predict)
print(recall)
f1 = f1_score(testLabelsGlobal, svm_predict)
print(f1)
auc = roc_auc_score(testLabelsGlobal, svm_predict)
print(auc)

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold

# models = []
seed=42
# num_trees = 10
# models.append(('LR',LogisticRegression(random_state = seed)))
# models.append(('LDA',LinearDiscriminantAnalysis()))
# models.append(('KNN',KNeighborsClassifier()))
# models.append(('CART',DecisionTreeClassifier(random_state = seed)))
# models.append(('RF',RandomForestClassifier(n_estimators = num_trees, random_state = seed)))
# models.append(('NB', GaussianNB()))

# results = []
# names = []
# for name, model in models: 
#     kfold =KFold(n_splits = 10, random_state = seed)
#     cv_results = cross_val_score(model, trainDataGlobal,trainLabelsGlobal, cv = kfold, scoring = 'accuracy')
#     results.append(cv_results)
#     names.append(name)
#     msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
#     print(msg)
    
model = LogisticRegression(random_state = seed)
model.fit(trainDataGlobal, trainLabelsGlobal)
lr_predict= model.predict(testDataGlobal)
acc_lr = model.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on LR: ",acc_lr)
precision = precision_score(testLabelsGlobal, lr_predict)
print(precision)
recall = recall_score(testLabelsGlobal, lr_predict)
print(recall)
f1 = f1_score(testLabelsGlobal, lr_predict)
print(f1)
auc = roc_auc_score(testLabelsGlobal, lr_predict)
print(auc)

model = LinearDiscriminantAnalysis()
model.fit(trainDataGlobal, trainLabelsGlobal)
lda_predict= model.predict(testDataGlobal)
acc_lda = model.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on LDA: ",acc_lda)
precision = precision_score(testLabelsGlobal, lda_predict)
print(precision)
recall = recall_score(testLabelsGlobal, lda_predict)
print(recall)
f1 = f1_score(testLabelsGlobal, lda_predict)
print(f1)
auc = roc_auc_score(testLabelsGlobal, lda_predict)
print(auc)

model = KNeighborsClassifier()
model.fit(trainDataGlobal, trainLabelsGlobal)
knn_predict= model.predict(testDataGlobal)
acc_knn = model.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on knn: ",acc_knn)
precision = precision_score(testLabelsGlobal, knn_predict)
print(precision)
recall = recall_score(testLabelsGlobal, knn_predict)
print(recall)
f1 = f1_score(testLabelsGlobal, knn_predict)
print(f1)
auc = roc_auc_score(testLabelsGlobal, knn_predict)
print(auc)

model = DecisionTreeClassifier(random_state = seed)
model.fit(trainDataGlobal, trainLabelsGlobal)
dc_predict= model.predict(testDataGlobal)
acc_dc = model.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on DC: ",acc_dc)
precision = precision_score(testLabelsGlobal, dc_predict)
print(precision)
recall = recall_score(testLabelsGlobal, dc_predict)
print(recall)
f1 = f1_score(testLabelsGlobal, dc_predict)
print(f1)
auc = roc_auc_score(testLabelsGlobal, dc_predict)
print(auc)

# #confusion metrix
# cm_rf = confusion_matrix(testLabelsGlobal,rf_predict)
# plot_confusion_matrix(cm_rf, normalize= False,target_names=[0,1], title = "Confusion Matrix RF")

# from sklearn.metrics import average_precision_score
# precision = average_precision_score(testLabelsGlobal, rf_predict)
# print(precision)

modelRF=RandomForestClassifier(n_estimators= 800, random_state=42,max_depth=80,min_samples_split=2,min_samples_leaf=1)
modelRF.fit(trainDataGlobal, trainLabelsGlobal)
rf_predict= modelRF.predict(testDataGlobal)
y_score = modelRF.predict_proba(testDataGlobal)
acc_rf = modelRF.score(testDataGlobal,testLabelsGlobal)
print("Accuracy on RF: ",acc_rf)

precision = precision_score(testLabelsGlobal, rf_predict)
print(precision)


recall = recall_score(testLabelsGlobal, rf_predict)
print(recall)


f1 = f1_score(testLabelsGlobal, rf_predict)
print(f1)


auc = roc_auc_score(testLabelsGlobal, rf_predict)
print(auc)

def calculate_performace(test_num, pred_y,  labels):
    tp =0
    fp = 0
    tn = 0
    fn = 0
    for index in range(test_num):
        if labels[index] ==1:
            if labels[index] == pred_y[index]:
                tp = tp +1
            else:
                fn = fn + 1
        else:
            if labels[index] == pred_y[index]:
                tn = tn +1
            else:
                fp = fp + 1 
    acc = float(tp + tn)/test_num
    if tp == 0 and fp == 0:
        precision = 0
        sensitivity = float(tp)/ (tp+fn)
        specificity = float(tn)/(tn + fp)
    else:
        precision = float(tp)/(tp+ fp)
        sensitivity = float(tp)/ (tp+fn)
        specificity = float(tn)/(tn + fp)
        # MCC = float(tp*tn-fp*fn)/(np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return acc, precision, sensitivity, specificity #, MCC 



from keras.layers.core import Dropout, Activation
from keras.models import Sequential
from keras.layers import Dense
from keras import optimizers

model = Sequential()
model.add(Dense(input_dim=904, output_dim=600,init='he_normal'))
model.add(Activation('relu'))
model.add(Dropout(0.2))
# model.add(Dense(input_dim=1398, output_dim=800,init='he_normal'))
# model.add(Activation('relu'))
# model.add(Dropout(0.2))

#Initial model
model.add(Dense(input_dim=600, output_dim=300,init='he_normal'))
model.add(Activation('relu'))
model.add(Dropout(0.2))
model.add(Dense(input_dim=300, output_dim=150,init='he_normal'))
model.add(Activation('relu'))
model.add(Dropout(0.2))
model.add(Dense(input_dim=150, output_dim=50,init='he_normal'))
model.add(Activation('relu'))
model.add(Dropout(0.2))
model.add(Dense(input_dim=50, output_dim=2,init='glorot_normal'))
model.add(Activation('sigmoid'))

# model.add(Dense(input_dim=125, output_dim=100,init='he_normal'))
# model.add(Activation('relu'))
# model.add(Dropout(0.2))
# model.add(Dense(input_dim=100, output_dim=75,init='he_normal'))
# model.add(Activation('relu'))
# model.add(Dropout(0.2))
# model.add(Dense(input_dim=75, output_dim=40,init='he_normal'))
# model.add(Activation('relu'))
# model.add(Dropout(0.2))
# model.add(Dense(input_dim=40, output_dim=2,init='glorot_normal'))
# model.add(Activation('sigmoid'))

# sgd = optimizers.SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
# adam = optimizers.Adam(lr=0.5, decay=1e-6)
# model.compile(loss='binary_crossentropy', optimizer=adam)   

print(model.summary())
#new model

# x_train = np.array(trainDataGlobal)
# x_train = x_train.reshape(trainDataGlobal.shape[0],46,9).astype('float32')
# x_test = np.array(testDataGlobal).reshape(testDataGlobal.shape[0],46,9).astype('float32')

# from keras.models import Sequential
# from keras.layers import Dense
# from keras.layers import Flatten
# from keras.layers import Dropout
# from keras.layers.convolutional import Conv1D
# from keras.layers.convolutional import MaxPooling1D
# from keras.utils import to_categorical
# from keras.layers import GlobalAveragePooling1D

# model = Sequential()
# model.add(Conv1D(filters=256, kernel_size=2, activation='relu', input_shape=(46, 9)))
# model.add(Conv1D(filters=256, kernel_size=2, activation='relu'))

# model.add(MaxPooling1D(pool_size=2 ))
# model.add(Conv1D(filters=128, kernel_size=2, activation='relu'))
# model.add(Conv1D(filters=128, kernel_size=2, activation='relu'))
# # model.add(GlobalAveragePooling1D())
# model.add(Dropout(0.3))
# model.add(MaxPooling1D(pool_size=2 ))

# model.add(Flatten())
# model.add(Dense(500, activation='relu'))
# model.add(Dense(250, activation='relu'))
# model.add(Dense(100, activation='relu'))
# model.add(Dense(2, activation='softmax'))
# print(model.summary())

model.compile(loss='categorical_crossentropy',
                optimizer='adam', metrics=['accuracy'])

history = model.fit(X, Y_labels,
            batch_size=100,
            epochs=80,
            validation_split=0.1)

train_accuracy = history.history[ 'accuracy' ]
validation_accuracy = history.history[ 'val_accuracy' ]
# score = model.evaluate(testDataGlobal, testLabelsGlobal, verbose=0)

model.save('I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\model_seq_new.h5')
# print('Score: ',score)
proba = model.predict_classes(testDataGlobal,batch_size=200,verbose=True)
ae_y_pred_prob = model.predict_proba(testDataGlobal,batch_size=100,verbose=True)
Y_test = lb.inverse_transform(testLabelsGlobal)
acc, precision, sensitivity, specificity = calculate_performace(len(Y_test), proba,  Y_test)
# model.save_weights("model.h5")
# from keras.utils.vis_utils import plot_model
# plot_model(model, to_file='model_plot.png', show_shapes=True, show_layer_names=True)
# import os
# print(os.environ["PATH"])
from sklearn.metrics import roc_curve, auc
fpr, tpr, auc_thresholds = roc_curve(testLabelsGlobal[:,1], ae_y_pred_prob[:,1])
auc_score = auc(fpr, tpr)

from sklearn.metrics import precision_recall_curve
precision1, recall, pr_threshods = precision_recall_curve(testLabelsGlobal[:,1], ae_y_pred_prob[:,1])
aupr_score = auc(recall, precision1)


print("Precision: ",precision)
print("Recall: ",sensitivity)


from sklearn.metrics import f1_score
f1 = f1_score(Y_test, proba)
print("F1 Score: ",f1)

from sklearn.metrics import roc_curve, auc  
fpr, tpr, auc_thresholds = roc_curve(Y_test, ae_y_pred_prob[:,1])
auc_score = auc(fpr, tpr)
print("AUC: ",auc_score)
print("accuracy: ",acc)
print('AUPR: ',aupr_score)


#test Dataset

alignment_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS1\\merged_feature.csv'
# f_path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\NDD Dataset\\ds1_pdb_ndd\\merged_pdb_vector_ds1.csv'
df_alignment = pd.read_csv(alignment_path)
df_alignment = df_alignment.set_index('Unnamed: 0')
# df_alignment= df_alignment.fillna(1)
# df2 = (df_alignment - df_alignment.values.min()) / (df_alignment.values.max()-df_alignment.values.min())
# df_alignment = df2
# f_path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\merged_vector_new.csv'
f_path ='I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS1_V2\\merged_feature.csv'

            
# f_path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\df_seq_merge_vector.csv'
# path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\ddi_seq_cardio.csv'
path = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\DS1_V2\\ddi.csv'

#preparing the dataset without cardio vascular to train
ddi_seq = pd.read_csv(path).set_index('Unnamed: 0')
ddi = pd.read_csv(path2).set_index('Unnamed: 0')
drugList = list(ddi_seq.keys())
dr_list2 = list(pd.read_csv(path2).set_index('Unnamed: 0').keys())
with open(f_path, 'w',newline='') as f:
    writer = csv.writer(f)
    for ii in range(len(dr_list2)):
        Id1 = dr_list2[ii]
        if Id1 in drugList:
            continue
        for jj in range(ii+1,len(dr_list2)):
            Id2 = dr_list2[jj]
            if Id2 in drugList:
                continue
            firstList = df_alignment[dr_list2[ii]]
            scndList = df_alignment[dr_list2[jj]]
            output = list(firstList - scndList)
            row = [dr_list2[ii],dr_list2[jj]] + output
            writer.writerow(row)

path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\ddi_merged_new.csv'
# drugList = list(ddi_al.keys())
with open(path, 'w',newline='') as f:
    writer = csv.writer(f)
    for ii in range(len(dr_list2)):
        Id1 = dr_list2[ii]
        if Id1 in drugList:
            continue
        for jj in range(ii+1,len(dr_list2)):
            Id2 = dr_list2[jj]
            if Id2 in drugList:
                continue
            row = [dr_list2[ii],dr_list2[jj]] + [ddi[dr_list2[ii]][dr_list2[jj]]]
            writer.writerow(row)


#load models
from tensorflow import keras
model = keras.models.load_model('I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\model_seq.h5')
feature_space = pd.read_csv(f_path,header=None)
FS = pd.DataFrame(feature_space)

FS = FS.replace(-np.inf, np.nan)
FS = FS.replace(np.inf, np.nan)
FS = FS.fillna(FS.mean())


X = FS.drop([0,1],axis=1)
proba = model.predict_classes(X,batch_size=50,verbose=True)
proba= proba.tolist()
labels = pd.read_csv(path,header=None)

f_path = 'I:\\Study\\Research WVU\\PDB\\cardioVascular Disease\\matrices\\comparison_predcition_seq_cardio_new.csv'
Columns = ['Drug Id1','Drug Id2','label','predicted_label']
 
with open(f_path, 'w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(Columns)
    for ii in range(len(labels)):
        if proba[ii] == 1:
            
            row = [labels.iloc[ii][0],labels.iloc[ii][1],labels.iloc[ii][2],0] 
            writer.writerow(row)
        else:
            row = [labels.iloc[ii][0],labels.iloc[ii][1],labels.iloc[ii][2],1] 
            writer.writerow(row)
            


from sklearn.model_selection import KFold

kfold = KFold(n_splits=10,shuffle=True,random_state=42)

fold_no = 1
acc_per_fold = []
loss_per_fold = []

for train,test in kfold.split(X,Y_labels):

    model = Sequential()
    model.add(Dense(input_dim=662, output_dim=600,init='glorot_normal'))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(input_dim=600, output_dim=400,init='glorot_normal'))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(input_dim=400, output_dim=200,init='glorot_normal'))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(input_dim=200, output_dim=2,init='glorot_normal'))
    model.add(Activation('sigmoid'))
    sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    
 
    # Generate a print
    print('------------------------------------------------------------------------')
    print(f'Training for fold {fold_no} ...')
    
    model.compile(loss='binary_crossentropy', optimizer=sgd,metrics=['accuracy'])  

    history = model.fit(X.iloc[train], Y_labels[train],
            batch_size=100,
            epochs=50,
            validation_split=0,
            verbose=1)
    
    # Generate generalization metrics
    scores = model.evaluate(X.iloc[test], Y_labels[test], verbose=0)
    print(f'Score for fold {fold_no}: {model.metrics_names[0]} of {scores[0]}; {model.metrics_names[1]} of {scores[1]*100}%')
    acc_per_fold.append(scores[1] * 100)
    loss_per_fold.append(scores[0])

    # Increase fold number
    fold_no = fold_no + 1


# == Provide average scores ==
print('------------------------------------------------------------------------')
print('Score per fold')
for i in range(0, len(acc_per_fold)):
  print('------------------------------------------------------------------------')
  print(f'> Fold {i+1} - Loss: {loss_per_fold[i]} - Accuracy: {acc_per_fold[i]}%')
print('------------------------------------------------------------------------')
print('Average scores for all folds:')
print(f'> Accuracy: {np.mean(acc_per_fold)} (+- {np.std(acc_per_fold)})')
print(f'> Loss: {np.mean(loss_per_fold)}')
print('------------------------------------------------------------------------')


from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report

tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-2,1e-3, 1e-4],
                     'C': [1, 10, 100, 1000]}]
scores = ['precision', 'recall']
        
for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

    clf = GridSearchCV(
        SVC(), tuned_parameters, scoring='%s_macro' % score
    )
    clf.fit(trainDataGlobal, trainLabelsGlobal)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = testLabelsGlobal, clf.predict(testDataGlobal)
    print(classification_report(y_true, y_pred))
    print()       
        
        