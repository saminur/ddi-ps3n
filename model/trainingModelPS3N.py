# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 23:59:34 2023

@author: Saminur Islam
"""

from pandas import read_excel
import numpy as np
import pandas as pd
import csv 

f_path = 'merged_feature.csv'

path = 'ddi_pdb.csv'

from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
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

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras import optimizers
from tensorflow.keras import regularizers

from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.models import model_from_json

model = Sequential()
model.add(Dense(500, activation='relu',input_dim=904,kernel_initializer='glorot_uniform',kernel_regularizer=regularizers.l2(0.1)))
model.add(BatchNormalization())
# model.add(Dense(400, activation='relu',input_dim=600,kernel_initializer='glorot_uniform'))
model.add(Dropout(0.1))
model.add(Dense(300, input_dim=500,activation='relu',kernel_initializer='glorot_uniform',kernel_regularizer=regularizers.l2(0.01)))
model.add(BatchNormalization())
model.add(Dropout(0.1))
model.add(Dense(100, input_dim=300,activation='relu',kernel_initializer='glorot_uniform'))
# model.add(Dropout(0.2))
# model.add(Dense(75,input_dim=300,activation='relu',kernel_initializer='glorot_normal',kernel_regularizer=regularizers.l2(0.1)))
model.add(BatchNormalization())
# model.add(Dropout(0.2))
model.add(Dense(2, input_dim=100))
model.add(Activation('sigmoid'))

print(model.summary())
optimizer = optimizers.Adam(learning_rate=0.0001)

model.compile(loss='categorical_crossentropy',
                optimizer= optimizer, metrics=['accuracy'])

history = model.fit(X, Y_labels,
            batch_size=512,
            epochs=50,
            validation_split=0.1)

train_accuracy = history.history[ 'accuracy' ]
validation_accuracy = history.history[ 'val_accuracy' ]

#save model
# serialize model to JSON
pathJson = 'trainedModel\\model.json'
pathweights = 'trainedModel\\ddi_pdb_ps.h5'
model_json = model.to_json()
with open(pathJson, "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
model.save_weights(pathweights)
print("Saved model to disk")

# Create the plot
epochs = [i for i in range(1, 50 + 1)]
plt.figure()
plt.plot(epochs, train_accuracy, label='Training Accuracy')
plt.plot(epochs, validation_accuracy, label='Validation Accuracy')

# Add labels and a legend
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.title('Training and Validation Accuracy')
plt.legend()

# Display or save the plot
plt.show()

# score = model.evaluate(testDataGlobal, testLabelsGlobal, verbose=0)

# model.save('H:\\Research Work\\TrinetXDataSet\\New Results for Jounral\\trainedModel\\ddi_pdb_ps.h5')

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

proba = model.predict(testDataGlobal,batch_size=200,verbose=True)
proba = np.argmax(proba, axis=1)
ae_y_pred_prob = model.predict(testDataGlobal,batch_size=100,verbose=True)
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

#check the model performance on the

#load model 
json_file = open(pathJson, 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
# load weights into new model
loaded_model.load_weights(pathweights)

optimizer = optimizers.Adam(learning_rate=0.0001)

loaded_model.compile(loss='categorical_crossentropy',
                optimizer= optimizer, metrics=['accuracy'])



ndd_ds1_path = 'merged_protein_pdb_vector_ds1.csv'
path = 'pdb_ddi_vector_ds1.csv'
feature_space = pd.read_csv(ndd_ds1_path,header=None)
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

#reshaped array into Sample size X 904 
reshaped_array = np.concatenate((X.values, np.empty((X.values.shape[0], 490))), axis=1)

proba = loaded_model.predict(pd.DataFrame(reshaped_array),batch_size=200,verbose=True)

path_similarity = 'merged_protein_pdb_ds1.csv'

path_simi_cosine_pdb = 'ds1_pdb_ndd\\pdb matrices\\pdb_cosine_avg.csv'
path_simi_cosine_protein = 'ds1_protein_ndd\\protein_cosine_avg.csv'

path_found_ddi = 'New Results for Jounral\\new_found_ddi_ndd.csv'
path_drugs = 'drug_products.txt' 
df_drug_info = pd.read_csv(path_drugs, sep="\t",header=0).drop("id", axis=1)
df_path_found_ddi = pd.read_csv(path_found_ddi)
df_path_smilarity = pd.read_csv(path_similarity).set_index('Unnamed: 0')
df_protein_cosine = pd.read_csv(path_simi_cosine_protein).set_index('Unnamed: 0')
df_pdb_cosine = pd.read_csv(path_simi_cosine_pdb).set_index('Unnamed: 0')

dict_drug_band = {}
for index,row in df_drug_info.iterrows():
    # print(row['drug_id'])
    if row['drug_id'] not in dict_drug_band:
        dict_drug_band[row['drug_id']] = row['brand_name']
    # print(df_drug_info[i]['drug id'])
    # break;
    
dict_similarity = {}

path = 'network_ndd_all_similarity.csv'
columns = ['Drug ID1', 'Drug Name1', 'Drug ID2', 'Drug Name2', 'Interaction','fushion similarity','protein similarity','pdb similarity','Predicted Prob']
df_path_smilarity['DB00115']['DB00136']

with open(path, 'w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(columns)
    for index,row in df_path_found_ddi.iterrows():
        DrugName1 = dict_drug_band.get(row['drug1'])
        DrugName2 = dict_drug_band.get(row['Drug2'])
        similarity = df_path_smilarity[row['drug1']][row['Drug2']]
        protein_similarity = df_protein_cosine[row['drug1']][row['Drug2']]
        pdb_similarity = df_pdb_cosine[row['drug1']][row['Drug2']]
        predicted_prob = max (proba[index][0],proba[index][1])
        write_row = [row['drug1'],DrugName1,row['Drug2'],DrugName2,row['Interactions'],similarity,protein_similarity,pdb_similarity,predicted_prob]
        writer.writerow(write_row)
        

df_network_ndd = pd.read_csv(path)    

proba = np.argmax(proba, axis=1)
ae_y_pred_prob = loaded_model.predict(pd.DataFrame(reshaped_array),batch_size=100,verbose=True)
Y_test = lb.inverse_transform(Y_labels)
acc, precision, sensitivity, specificity = calculate_performace(len(Y_test), proba,  Y_test)
# model.save_weights("model.h5")
# from keras.utils.vis_utils import plot_model
# plot_model(model, to_file='model_plot.png', show_shapes=True, show_layer_names=True)
# import os
# print(os.environ["PATH"])
from sklearn.metrics import roc_curve, auc
fpr, tpr, auc_thresholds = roc_curve(Y_labels[:,1], ae_y_pred_prob[:,1])
auc_score = auc(fpr, tpr)

from sklearn.metrics import precision_recall_curve
precision1, recall, pr_threshods = precision_recall_curve(Y_labels[:,1], ae_y_pred_prob[:,1])
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

#found new interactions from NDD
labels = pd.read_csv(path,header=None)
filename = 'new_found_ddi_ndd.csv'
columns = ["drug1","Drug2","Interactions"]
proba = np.array(proba)
proba[0]
with open(filename, 'w',newline='') as csvfile: 
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(columns) 
    for index, row in labels.iterrows():
        if row[2] == 0 and proba[index] == 1:
            print(row[0])
            csvwriter.writerow([row[0],row[1],1])


    