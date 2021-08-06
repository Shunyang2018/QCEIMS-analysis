import pandas
from pandas.plotting import scatter_matrix
import numpy
import sys
numpy.set_printoptions(threshold=sys.maxsize)
import pickle
import math
import matplotlib.pyplot
#import jdxtool.py
from jdxtool import *

with open('./index_insilico_ref.pkl','rb') as temp_file:
    temp_as_list=pickle.load(temp_file)




#temp_as_list=[[0,numpy.array([0.1,0.5,0.9,0,0,1]),numpy.array([0.1,0.5,0.9,0,0,1])],
#[1,numpy.array([0.1,0.5,0.9,0,0,1]),numpy.array([0.0,0.5,0.9,0,0,1])]]

temp=numpy.array(temp_as_list)
print(numpy.shape(temp))
print(numpy.shape(temp[:][:]))
#for item in temp:
#    print(item)
    #print(len(item[1]))
#    hold=input('hold')

#subtemp=temp[0:2]
subtemp=temp
print(subtemp)
#hold=input('hold')
#the ordering of things


#get jaccard overprediction, jaccard underprediction, jaccard intersection, jaccard overprediction minus underprediction
def get_tanimoto_list_for_all(x_list,y_list,how):
    tanimoto_list=list()

    for nth_compound in range(0,len(x_list)):
        temp_tanimoto_score=tanimoto_sim(x_list[nth_compound],y_list[nth_compound],how)
        tanimoto_list.append(temp_tanimoto_score)

    tanimoto_array=numpy.array(tanimoto_list)

    return tanimoto_array
    


    

#histogram those.... just see whats cooking




#define function that transforms vector to 1,0 where 1 if value is non-zero, 0 otherwise
def cast_float_list_to_binary(temp_list):
    temp_list_binary=list()
    for temp_spectrum in temp_list:
        #temp_spectrum_binary=list()
        print(temp_spectrum)
        #hold=input('hold')
        temp_binary_spectrum=numpy.ceil(temp_spectrum)
        #for j in range(0,len(temp_list[i])):
            #if the value at this index is not zero, add the number 1, else add 0
        #    if (j > 0):
        #        temp_spectrum_binary.append(1)
        #    elif (j==0):
        #        temp_spectrum_binary.append(0)
        #    else:
        #        input('error')
        #typecast array to make things like what we received in pikle file
        #temp_array=numpy.array(temp_spectrum_binary)
        temp_binary_spectrum=temp_binary_spectrum.astype(int)
        temp_list_binary.append(temp_binary_spectrum)
    temp_array_binary=numpy.array(temp_list_binary)
    return temp_array_binary


temp_again=subtemp[:,1]
#print(temp_again)
#hold=input('asdf')
my_binary_in_silico=cast_float_list_to_binary(subtemp[:,1])
my_binary_experimental=cast_float_list_to_binary(subtemp[:,2])
#my_binary.astype(int)
#print(my_binary_in_silico,my_binary_experimental)

my_jaccard_score_array=get_tanimoto_list_for_all(my_binary_in_silico,my_binary_experimental,'intersection')
#print(my_tanimoto_score_list)
my_binary_overprediction_score_array=get_tanimoto_list_for_all(my_binary_in_silico,my_binary_experimental,'overprediction')
#print(my_tanimoto_score_list)
my_binary_underprediction_score_array=get_tanimoto_list_for_all(my_binary_in_silico,my_binary_experimental,'underprediction')
#take the difference of the over and underprediciton scores, see which is dominating
my_binary_over_vs_under_array=numpy.subtract(my_binary_overprediction_score_array,my_binary_underprediction_score_array)

#then we can call the tanimoto function that we define here, because the jaccard score is just a special case of the tanimoto

def histogram_some_score_list(temp_list):
    matplotlib.pyplot.hist(temp_list,bins=numpy.arange(-1,1.05,0.05))
    matplotlib.pyplot.show()


histogram_some_score_list(my_jaccard_score_array)
histogram_some_score_list(my_binary_overprediction_score_array)
histogram_some_score_list(my_binary_underprediction_score_array)
histogram_some_score_list(my_binary_over_vs_under_array)
#get the tanimoto overprediction, tanimoto underprediction, tanimoto intersection, tanimoto overprediciton minus underprediction
#histogram those.... just see whats cooking


my_in_silico=(subtemp[:,1])
my_experimental=(subtemp[:,2])
#my_binary.astype(int)
#print(my_binary_in_silico,my_binary_experimental)

my_tanimoto_score_array=get_tanimoto_list_for_all(my_in_silico,my_experimental,'intersection')
#print(my_tanimoto_score_list)
my_overprediction_score_array=get_tanimoto_list_for_all(my_in_silico,my_experimental,'overprediction')
#print(my_tanimoto_score_list)
my_underprediction_score_array=get_tanimoto_list_for_all(my_in_silico,my_experimental,'underprediction')

my_over_vs_under_array=numpy.subtract(my_overprediction_score_array,my_underprediction_score_array)
#my_under_vs_over_array=numpy.subtract(my_underprediction_score_array,my_overprediction_score_array)



histogram_some_score_list(my_tanimoto_score_array)
histogram_some_score_list(my_overprediction_score_array)
histogram_some_score_list(my_underprediction_score_array)
histogram_some_score_list(my_over_vs_under_array)
#histogram_some_score_list(my_under_vs_over_array)


def get_weighted_dot_for_all(x_list,y_list):
    weighted_dot_list=list()

    for nth_compound in range(0,len(x_list)):
        temp_weighted_dot_score=weighted_dot_sim(x_list[nth_compound],y_list[nth_compound])
        weighted_dot_list.append(temp_weighted_dot_score)

    weighted_dot_array=numpy.array(weighted_dot_list)

    return weighted_dot_array

#get two more things....
my_binary_weighted_dot=get_weighted_dot_for_all(my_binary_in_silico,my_binary_experimental)
my_weighted_dot=get_weighted_dot_for_all(my_in_silico,my_experimental)

def histogram_some_weighted_dot_list(temp_list):
    matplotlib.pyplot.hist(temp_list,bins=20)
    matplotlib.pyplot.show()

histogram_some_weighted_dot_list(my_binary_weighted_dot)
histogram_some_weighted_dot_list(my_weighted_dot)


#create panda for mega scatter plot
my_panda=pandas.DataFrame({
'bin inter':my_jaccard_score_array,
'bin over':my_binary_overprediction_score_array,
'bin under':my_binary_underprediction_score_array,
'bin over-under':my_binary_over_vs_under_array,
'cts inter':my_tanimoto_score_array,
'cts over':my_overprediction_score_array,
'cts under':my_underprediction_score_array,
'cts over-under':my_over_vs_under_array,
'bin wt dot':my_binary_weighted_dot,
'cts wt dot':my_weighted_dot
})
attributes_of_interest=['bin inter',
'bin over',
'bin under',
'bin over-under',
'cts inter',
'cts over',
'cts under',
'cts over-under',
'bin wt dot',
'cts wt dot']

scatter_matrix(my_panda[attributes_of_interest],figsize=(20,20),alpha=0.1)
matplotlib.pyplot.show()


'''
my_jaccard_score_array=get_tanimoto_list_for_all(my_binary_in_silico,my_binary_experimental,'intersection')
my_binary_overprediction_score_array=get_tanimoto_list_for_all(my_binary_in_silico,my_binary_experimental,'overprediction')
my_binary_underprediction_score_array=get_tanimoto_list_for_all(my_binary_in_silico,my_binary_experimental,'underprediction')
my_binary_over_vs_under_array=numpy.subtract(my_binary_overprediction_score_array,my_binary_underprediction_score_array)

my_tanimoto_score_array=get_tanimoto_list_for_all(my_in_silico,my_experimental,'intersection')
my_overprediction_score_array=get_tanimoto_list_for_all(my_in_silico,my_experimental,'overprediction')
my_underprediction_score_array=get_tanimoto_list_for_all(my_in_silico,my_experimental,'underprediction')
my_over_vs_under_array=numpy.subtract(my_overprediction_score_array,my_underprediction_score_array)

my_binary_weighted_dot=get_weighted_dot_for_all(my_binary_in_silico,my_binary_experimental)
my_weighted_dot=get_weighted_dot_for_all(my_in_silico,my_experimental)
'''

#print(my_binary_underprediction_score_array.sort())
print(numpy.sort(my_binary_underprediction_score_array))

#%%

import matplotlib.pyplot as plt
plt.clf()
fig,boxplot = plt.subplots()
plt.style.use('seaborn-deep')

hist = [my_tanimoto_score_array,my_overprediction_score_array,my_underprediction_score_array]
l = ['inter','over','under']
col = ['tab:green','tab:blue','tab:orange']
# bins = np.linspace(0, 100, 25)/100
# plt.hist(hist, label=l, color=col)
# plt.hist(my_tanimoto_score_array, alpha=0.5, label='inter',rwidth=0.4)
# plt.hist(my_overprediction_score_array, alpha=0.5, label='over',rwidth=0.4)
plt.hist(my_over_vs_under_array,rwidth=0.8,color='tab:blue')



boxplot.set_ylabel('Count',fontsize=20,fontweight='bold', fontname="Arial")
boxplot.set_xlabel('Over minus Under',fontsize=20,fontweight='bold', fontname="Arial")

boxplot.spines['bottom'].set_linewidth(2)
boxplot.spines['top'].set_linewidth(2)
boxplot.spines['left'].set_linewidth(2)
boxplot.spines['right'].set_linewidth(2)
#plt.legend(loc=0)
for tick in boxplot.xaxis.get_major_ticks():
    tick.label.set_fontsize(16) 
for tick in boxplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(16) 
    
# plt.tight_layout()#won't trim the figure
boxplot.grid(False)

#%%


import pickle
import matplotlib.pyplot as plt

fig = plt.figure()
boxplot = fig.add_subplot(111)
plt.style.use('seaborn-deep')

hist = [my_tanimoto_score_array,my_overprediction_score_array,my_underprediction_score_array]

import pickle

with open('hist.pkl', 'wb') as f:
        pickle.dump(hist, f)

        
#%%
l = ['inter','over','under']
col = ['tab:green','tab:blue','tab:orange']
# bins = np.linspace(0, 100, 25)/100
# plt.hist(hist, label=l, color=col)
# plt.hist(my_tanimoto_score_array, alpha=0.5, label='inter',rwidth=0.4)
# plt.hist(my_overprediction_score_array, alpha=0.5, label='over',rwidth=0.4)
#plt.hist(my_over_vs_under_array,rwidth=0.8,color='tab:blue')

boxplot.scatter(my_weighted_dot, my_jaccard_score_array,label='inter',s=5,color='tab:green',alpha=0.6)
boxplot.scatter(my_weighted_dot, my_binary_overprediction_score_array,label='over',s=5,color='tab:blue')
boxplot.scatter(my_weighted_dot, my_binary_underprediction_score_array,label='under',s=5,color='tab:orange',alpha=0.6)

boxplot.set_ylabel('Jaccard Index',fontsize=20,fontweight='bold', fontname="Arial")
boxplot.set_xlabel('Weighter dot score',fontsize=20,fontweight='bold', fontname="Arial")

boxplot.spines['bottom'].set_linewidth(2)
boxplot.spines['top'].set_linewidth(2)
boxplot.spines['left'].set_linewidth(2)
boxplot.spines['right'].set_linewidth(2)
#plt.legend(loc=0)
for tick in boxplot.xaxis.get_major_ticks():
    tick.label.set_fontsize(16) 
for tick in boxplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(16) 
    
# plt.tight_layout()#won't trim the figure
boxplot.grid(False)



box = boxplot.get_position()
boxplot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
boxplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()





#%%

import pickle
import matplotlib.pyplot as plt

with open('hist.pkl', 'rb') as f:        
        hist = pickle.load(f)
import seaborn as sns


# sns.jointplot(my_weighted_dot, my_binary_overprediction_score_array, kind='scatter')
ax = sns.jointplot(my_weighted_dot,my_binary_underprediction_score_array-my_binary_overprediction_score_array, kind='scatter',color='tab:blue')
# sns.jointplot(my_weighted_dot, my_jaccard_score_array, kind='scatter')



ax.ax_joint.set_xlabel('Weighter dot score',fontsize=20,fontweight='bold', fontname="Arial")
ax.ax_joint.set_ylabel('Under - Over prediction',fontsize=20,fontweight='bold', fontname="Arial")
plt.tight_layout()
plt.show()



#%%

index = [i[0] for i in temp_as_list]

df = pandas.DataFrame({'index':index,'dot':my_weighted_dot
                       , 'under':my_binary_underprediction_score_array,
                       'over':my_binary_overprediction_score_array,'inter':my_jaccard_score_array})


df.to_csv('./jaccard.csv')























