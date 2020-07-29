#!/usr/bin/python
import numpy as np
import pandas as pd
import sklearn
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
import random
import pickle
import re
import sys

indir=sys.argv[1].strip()
muts=sys.argv[2].strip()
model_select=sys.argv[3].strip()

with open(model_select, 'rb') as f:
    model = pickle.load(f)
    for x in muts:
        outdir=indir+"/"+x
        x_read=pd.read_csv(outdir+"/DMBS_feature.csv")
        predict_pro=model.predict_proba(x_read)[:,1]
        predict=model.predict(x_read)
        file_for_write=open(outdir+"/pred.txt","w+")
        text = "Score:{sc}".format(sc=float(predict_pro))
        text += "\n"
        text += "Pred:{pr}".format(pr=int(predict))
        file_for_write.write(text)
        file_for_write.close()

a=""
for x in muts:
    a+=x
    a+="\n"
    outdir=indir+"/"+x
    f=open(outdir+"/pred.txt","r")
    part_a=f.read()
    f.close()
    a+=part_a
    a+="\n"

w=open(indir+"/pred_all.txt","w+")
w.write(a)
w.close()
