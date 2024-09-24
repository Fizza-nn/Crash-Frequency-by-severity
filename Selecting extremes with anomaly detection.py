import pandas as pd
from pycaret.anomaly import *
import numpy as np
from scipy.stats import genextreme, genpareto
import math
import statsmodels.api as sm
import pylab

FILE = 'data.xlsx'

df = pd.read_excel(FILE)          
df2 = df[['mttc','deltav']]
dataset['mttc'] = -1*dataset['mttc']
dataset = df2.loc[(df['deltav'])]
dataset = df2.loc[(df2['mttc'])]

details = setup(dataset)  #setup() function must be called before using create_model()
model = create_model('iforest',fraction=0.05) 

print(model)
#This function returns a trained model object.
model_predictions = predict_model(model, data = dataset)
X = np.array(model_predictions)

anomaly = []

for i in range (len(X)):
    if (X[i][2]==1):
        anomaly.append(X[i,0:2])

n = 0
anomaly = np.array(anomaly)

