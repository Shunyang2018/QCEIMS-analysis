#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:36:46 2019

@author: shunyang
"""
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('/Users/shunyang/Downloads/QCEIMS_conclusion_v7.xlsx', header = 11, sheet_name = 'graphs')


#boxplot = df.boxplot(column=['COS'], by=['RBN'])
boxplot = df.boxplot(column=['Cos'], by=['RBN'])
boxplot.get_figure().gca().set_title("")
boxplot.set_ylabel('Cos similarity score',fontsize=15)
boxplot.set_xlabel('RBN',fontsize=17)
boxplot.set_xlim(0.5,10.5)

plt.suptitle('')
plt.savefig('cos_rbn.png', dpi=600)
#%%



df = pd.read_csv('/Users/shunyang/Box/qceims_input/conformer-2-nonene/qceims_gmmx_results/gmmx.csv')

ax = df.hist(column='Cos', grid=False, color='b',normed=True)
ser = df.Cos.tolist()
# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]  
xmin, xmax = min(xt), max(xt)  
lnspc = np.linspace(xmin, xmax, len(df.Cos.tolist()))
# lets try the normal distribution first
m, s = stats.norm.fit(ser) # get mean and standard deviation  
pdf_g = stats.norm.pdf(lnspc, m, s) # now get theoretical values in our interval  
plt.plot(lnspc, pdf_g, label="Norm", color='g') # plot it

ax[0].item().set_ylabel('Frequency ',fontsize=15)
ax[0].item().set_xlabel('Cos similarity score',fontsize=15)
ax = ax[0].item()
ax.get_figure().gca().set_title("")
plt.show()
plt.savefig('conformers_histogram.png', dpi=600)
#%%

#energy plot
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import spline

df = pd.read_excel('/Users/shunyang/Files/qceims reference/QCEIMS_SI.xlsx',sheet_name='2-nonene conformers')
energy = df['om2 energy (kcal/mol)'].tolist()
x = df['Index'].tolist()
xnew = np.linspace(1,119,300)
energy_smooth = spline(x, energy ,xnew)
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(xnew, energy_smooth, c ='b')
ax.plot(x, energy, c ='b')
ax.scatter(x, energy, c ='r', linewidths=0.1,label = 'Energy of conformers')
#ax = df.hist('om2 energy (kcal/mol)', grid=False, color='b',normed=True)
ax.set_xlabel('Molecule',fontsize=20)
ax.set_ylabel('Single point energy (kCal/Mol)',fontsize=15)
ax.legend(fontsize=12)
ax.set_xlim(0,100)
plt.grid(ls = '--')
plt.savefig('conformers_energy.png', dpi=400)
plt.show()

#%%
'energy & Cos'
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_excel('/Users/shunyang/Files/qceims reference/QCEIMS_SI.xlsx',sheet_name='2-nonene conformers')
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)

ax = df.plot(x='om2 energy (kcal/mol)',y='Cos', kind = 'scatter',ax=axes, linewidths=0.2, fontsize=10)
plt.grid(ls = '--')
ax.set_xlabel ('Single point energy (kCal/Mol)', fontsize=15)
ax.set_ylabel('Cos similarity score',fontsize=15)
#ax.set_ylim(750,840)
#ax.set_xlim(0,3)
plt.show()
plt.savefig('conformers_score.png', dpi=600)


#%%
'PHI&COS'
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('/Users/shunyang/Downloads/QCEIMS_conclusion_v7.xlsx', header = 11, sheet_name = 'graphs')



ax = df.plot(x='PHI',y='Cos', kind = 'scatter', linewidths=0.2, fontsize=10)
#df.plot(x='PHI',y='Dot', kind = 'scatter', linewidths=0.2,c='r', fontsize=10)

ax.get_figure().gca().set_title("")
ax.set_ylabel('Cos similarity score',fontsize=15)
ax.set_xlabel('PHI',fontsize=17)
ax.set_xlim(0,12)
ax.set_ylim(0,1.0)
plt.grid(ls = '--')
plt.suptitle('')
plt.savefig('cos_phi.png', dpi=600)


#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
df = pd.read_excel('/Users/shunyang/Documents/qceims reference/QCEIMS_SI.xlsx', sheet_name = 'data')

df['Cos'] = df['Cos']*1000
fig, ax = plt.subplots()
bins = np.linspace(0, 1000, 100)
ax.hist(df['Cos'].tolist(), bins, color='b', rwidth=0.8, alpha=0.5, label='Cos')
ax.hist(df['Dot'].tolist(), bins, color='g', rwidth=0.8, alpha=0.5, label='Dot')
# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]  
xmin, xmax = min(xt), max(xt)  
lnspc = np.linspace(xmin, xmax, len(df.Cos.tolist()))
# lets try the normal distribution first
m, s = stats.norm.fit(ser) # get mean and standard deviation  
pdf_g = stats.norm.pdf(lnspc, m, s) # now get theoretical values in our interval  
plt.plot(lnspc, pdf_g, label="Norm", color='g') # plot it

ax.set_ylabel('Frequency ',fontsize=15)
ax.set_xlabel('Dot similarity score',fontsize=15)
ax.set_xlim(0,1000)
ax.get_figure().gca().set_title("")
plt.show()
plt.savefig('cosall.png', dpi=600)






















