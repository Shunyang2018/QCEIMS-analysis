#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 15:31:58 2019

@author: shunyang
"""
#for NIST17
import pandas as pd
import numpy as np
import os
import sys
import Daphnis.methods

from jdxtool import cos_sim, dot_sim, readjdx, peak, match,  trans

"""
I/O section.

indexfile is the reference spectra database in xls file, which using a preset index as a marker to compare in-silico and reference spectra.
for example : GFn539.jdx, index is 539

output is the result file name in csv format

path should contain all in-silico spectra in jdx format
"""


#debug mode switch
debug = False
replace = 'result3TMS'#replace the name label
output = './result_clean_GFn.csv'
indexfile = './TMS_under700.xlsx'


#mast have[0],to remove OrderedDict in header
path = "./TMS/spectra" #file path that contains jdx

def compare(filename,df):
    '''
    compare computational and experimental mass spectrums
    Input
        filemane: str
            path and name of computational ms
        df: pandas.dataframe
            databases of experimental ms
    Output
        cos:float
            cos similarity
        dot:float
            dot product similarity like NIST search
    '''
    #get mass index in data base
    massindex = filename.split('/')[-1].split('.')[0]
    massindex = massindex.replace(replace,'')
    massindex = int(massindex)
    x = readjdx(filename)#read computational results

    if int(massindex) not in Index['No.'].tolist():
        print(str(massindex)+'is not in database')
        cos='NA'
        dot='NA'
        m='NA'
    else:
        spec=Index[Index['No.'] == massindex]
        if 'e' in spec['m/z'].iloc[0]:
#            print('eeeee')
            y = spec['m/z'].iloc[0].split(' ')[1:-1]
            y.insert(0,0)
        else:
#            print("no ee")
            y =  spec['m/z'].iloc[0].split()[1:-1]

        y = list(map(float, y))
        y = np.copy(y/np.amax(y)) #normalization to 100

        if debug:
            print('y norm', y)
            print(y[119])
            # x[207]=y[207]
            # x[92]=y[92]
            # x[73]=y[73]
            x[145]=y[145]
            x[137]=y[137]

        cos = cos_sim(x,y)
        dot = dot_sim(x,y)


        x_trans = trans(x)
        y_trans = trans(y)
        entropy = 1 - Daphnis.methods.distance(x_trans , y_trans , method="weighted_entropy",
                                            normalize_result=True, spectrum_refined=False,
                                            ms2_da=0.05)
        m = match(x,y)
        confusion_matrix = peak(x, y)
        # print(entropy)
        # print('CM',confusion_matrix)
    return cos, dot, m , confusion_matrix,entropy

#%%




print("reading databases")
Index = pd.read_excel(indexfile,sheet_name='3TMS',index_col=None, usecols = "A:P",header=0)
print("starting analysing")


#%%
# entropy
if debug:
    test = trans(readjdx('./spectra/GFn200.jdx'))
#%%
#main body
# def main():
os.chdir(path)
files= os.listdir(path)
Dot = []
Cos = []
Match = []
check = []
entropy = []
CM= np.empty(shape=[0,4])
# print('empty',CM.shape)
for file in files:

    # if os.path.splitext(file)[1] == ".jdx":
    if '3TMS' in os.path.splitext(file)[0]:
        print(file)
        cos, dot, m, cm,en= compare(file, Index)
        Cos.append(cos)
        Dot.append(dot)
        Match.append(m)
        check.append(file)
        CM = np.append(CM, cm, axis=0)
        entropy.append(en)
        if debug == True:

            break
        if len(check) != len(Cos):
            print(file)
            sys.exit(1)
#
#
dataframe = pd.DataFrame({'Index':check, 'Dot':Dot, 'Cos':Cos,
                          'Match':Match, 'entropy':entropy, 'TP': CM[:,0], 'FP':CM[:,1], 'FN':CM[:,2], 'TN':CM[:,3]})
dataframe.to_csv(output,index=False,sep=',')
if debug== True:

    print(dataframe)
dataframe.to_csv(output,index=False,sep=',')


# if __name__=="__main__":
#     main()


#%%
import matplotlib.pyplot as plt
fig,ax = plt.subplots()
dataframe['Cos'] = dataframe['Cos']/1000
dataframe['TPR'] = dataframe['TP']/(dataframe['TP'] + dataframe['FN'])
ax.hist(dataframe['TPR'],rwidth=0.8,label='True Positive Rate')
plt.legend()
plt.show()
plt.savefig('./TPRhistogram.png')
#%%
fig,ax = plt.subplots()

ax = dataframe.plot.scatter('TPR', 'Dot' , s=10)
df2 = dataframe[ (dataframe['TPR']>0.6) & (dataframe['Dot']>650) ]
df2.plot.scatter('TPR', 'Dot' , s=10,color='g', ax=ax )
df3 = dataframe[ (dataframe['TPR']>0.6) & (dataframe['Dot']<650) ]
df3.plot.scatter('TPR', 'Dot' , s=10,color='r', ax=ax )

ax.set_xlabel('TPR',fontsize=16,fontweight='bold', fontname="Arial")
ax.set_ylabel('Dot score', fontsize=16,fontweight='bold', fontname="Arial")
ax.spines['bottom'].set_linewidth(2)
ax.spines['top'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
plt.show()
plt.savefig('./TPRscatterdot.png', transparent=True)









