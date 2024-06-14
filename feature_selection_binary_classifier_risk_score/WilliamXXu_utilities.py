import numpy as np
import math
from decimal import *
import matplotlib.pyplot as plt
from scipy.stats import binomtest,t,sem
import pandas as pd
from itertools import product,chain,combinations
from collections.abc import Iterable, Sequence
from dataclasses import dataclass,field
from joblib import Parallel, delayed
import requests
import json
import pickle


def download_gene_list_harmonizome(url):
    #url = "https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/antigen+processing+and+presentation/KEGG+Pathways"
    #'https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/nf-kb+signaling+pathway/Biocarta+Pathways'

    dic = json.loads(requests.get(url).text)
    li=[]
    for k in dic['associations']:
        li.append(k['gene']['symbol'])
    print(len(li))
    return set(li)



def table_input(func,table,output=None,multiprocess=False,n_jobs=3):


    if multiprocess:
        generator=Parallel(n_jobs=n_jobs)(delayed(func)(**table.iloc[i,:].to_dict()) for i in range(table.shape[0]))
    else:
        generator=(func(**table.iloc[i,:].to_dict()) for i in range(table.shape[0]))
    li=list(generator)
    '''
    li=[]
    for index, row in table.iterrows():
        temp=func(**row.to_dict())
        li.append(temp)
    '''
    combined=dict()
    for d in li:
        if isinstance(d,dict):
            for k in d.keys():
                combined[k]=list()
    #combined={k: list() for d in li for k in d.keys()}

    for k,v in combined.items():
        for d in li:
            if isinstance(d,dict):
                v.append(d[k])
            else:
                v.append(None)
        '''
            for d in li:
        if isinstance(d,dict):
            for k,v in d.items():
                combined[k].append(v)
        '''

    table=pd.concat([table,pd.DataFrame.from_dict(combined)],axis=1)
    print(table)
    
    if output is not None:
        pickle.dump(table,open(output+'.pkl','wb'))
        table.to_csv(output+'.csv')

def input_table_creation(input_dict):
    li=list(product(*input_dict.values()))
    return pd.DataFrame(li,columns=list(input_dict.keys()))


def powerset(iterable,no_empty=False):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(int(no_empty),len(s)+1))

def stratify(outcome,column,select=None,include_or_exclude=True):
    '''
    general purpose: stratify the pandas DataFrame according to a certain column, with the potential to select a subgroup

    parameters:
    outcome - pandas DataFrame
    column - name of a column of outcome, to be strafified
    select - if not None, return selected subgroup according to this value. If None, return all subgroups as dicitionary with column values as key.

    returns:
    described in select already
    '''


    samples=outcome
    column_type=type(column)
    select_type=type(select)
    if select is None:
        di=dict()


        if column_type==str:
            all_possibilities=set(samples[column])
            for x in all_possibilities:
                di[x]=samples.loc[samples[column]==x]
        elif isinstance(column, Iterable):
            if isinstance(column,Sequence):
                pass
            else:
                column=list(column)

            all_possibilities=[]
            for x in column:
                all_possibilities.append(set(samples[column]))
            products=product(all_possibilities)
            for y in products:
                di[y]=stratify(samples,column,select=y)
        else:
            raise Exception('column data type error')
        return di
    else:
        if include_or_exclude is False:
            minus=stratify(samples,column,select=select)
            return samples.loc[samples.index.isin(minus.index)==False]
        else:
            if column_type==str:
                if isinstance(select, Iterable) and select_type!=str:
                    return samples.loc[samples[column].isin(set(select))]
                else:
                    return samples.loc[samples[column]==select]
            elif isinstance(column,Sequence):
                if isinstance(select, Sequence) and len(select)>=len(column):
                    pass
                else:
                    raise Exception('select data not compatiable for the column data')
                temp=samples
                for z in range(len(column)):
                    temp=stratify(temp,column[z],select[z])
                return temp
                #return iterative_selection(samples,column,column_type)
            else:
                raise Exception('column data type error')


def cross_entropy(prediction,actual):
    
    delta=1e-7
    return -np.sum(actual*np.log(prediction+delta))

@dataclass(frozen=False)
class progress_convergence():
    tol:float
    min_steps:int=10
    id_for_ploting:str=' '
    mode:str='numerical convergence'
    convergence_window:int=10
    plot_switch:bool=True
    plot_interval:int=10
    ci_exclude_target:float=None
    progress_list:list[float]=field(default_factory=list)
    input_list:list[float]=field(default_factory=list)
    current:int=0

#        self.convergence_window=congerence_window
    def __post_int__(self):
        self.convergence_list=list(np.linspace(1000, 1000000, num=self.convergence_window))


    def move(self,input): #random sample for numerical covergence, number of heads for binomial. Return bool of whether to continue moving
        self.current+=1
        self.input_list.append(input)
        data=self.input_list
        stopp=False

        avg=np.nanmean(data,axis=0)
  #      std_sample=np.std(data,axis=0)
        std_mean=sem(data,nan_policy='omit',axis=0)
        #print(np.sum(np.isnan(avg)))
        #print(np.sum(np.isnan(std_mean)))
        #print(data)


        if self.mode =='numerical convergence':

            self.convergence_list.append(avg)
            self.convergence_list.pop(0)
            
            #diff=max(self.convergence_list)-min(self.convergence_list)

            diff=np.divide(np.std(self.convergence_list,axis=0),np.mean(self.convergence_list,axis=0))
            progress=np.nanmax(diff)
            condition=progress<self.tol

        elif self.mode=='t' or self.mode =='t_ci_exclude':  
            
            if self.current>1:
                ci=t.interval(0.95, df=len(data)-1, loc=avg, scale=std_mean)
                #print(ci)
                diff=ci[1]-ci[0]
                #diff[~np.isnan(diff)]
                #print(type(diff))


                if self.mode == 't_ci_exclude': 
                    #condition=np.prod((self.ci_exclude_target>ci[1]) or (self.ci_exclude_target<ci[0]) or (diff<self.tol))
                    progress=np.nanmax(np.min([self.ci_exclude_target-ci[0],ci[1]-self.ci_exclude_target,diff-self.tol],axis=0))
                    condition=progress<0
                else:
                    progress=np.nanmax(diff)
                    condition=progress<self.tol
                    #print(progress)
            else:
                progress=np.nan
                condition=False

        elif self.mode=='binomial':
            if isinstance(input,pd.Series):
                raise Exception('binomial only for single input!')
            else:
                ci=binomtest(input,self.current).proportion_ci(method='wilsoncc')
                diff=ci.high-ci.low
                condition=diff<self.tol
                progress=diff
        else:
            raise Exception('mode error!')
        
        if condition and self.current>=self.min_steps:
            stopp=True

        self.progress_list.append(progress)
        '''
                print(progress)
        print(stopp)
        print(self.plot_switch)
        print(self.current%self.plot_interval==0)
        print(self.current)
        print(self.plot_interval)
        print(self.current%self.plot_interval)
        '''


        #print(self.current)
        #print(progress)
        #print(self.plot_switch)
        if stopp:
            if self.plot_switch:
                print('Convergence achieved!')
                self.ploting()
            return False
        else:
            if self.plot_switch and (self.current%self.plot_interval==0):
                #print('activated')
                #print('Current loss: '+str(self.progress_list[-1]))
                self.ploting()
            return True
        '''
        if abs(progress-self.previous_progress)<self.tol:
            self.consecutive_convergence+=1
        else:
            self.consecutive_convergence=0
        self.previous_progress=progress
        if self.consecutive_convergence>=self.consercutive_threshold:
        '''



    def ploting(self):
        #x=list(range(self.current))
        plt.plot(self.progress_list)
        print('Current loss: '+str(self.progress_list[-1]))
        plt.xlabel('Number of iterations')
        plt.ylabel('Result value')
        plt.title('Iteration Progress Plot of '+self.id_for_ploting)
        plt.show()

def pair_generator(li,option):
    #'no diagnoal', 'no rep', 'all'

    res=[]
    for x in range(len(li)):
        if option=='all':
            for y in li:
                res.append((li[x],y))
        elif option =='no diagnoal':
            for y in li:
                if li[x]!=y:
                    res.append((li[x],y))
        elif option =='no rep':
            for y in range(x+1,len(li)):
                res.append((li[x],li[y]))
    return res

def hue_pair_generator(pairs,hues,mode='hue_constant'):

    res=[]
    for x in hues:
        for y in pairs:
            if mode=='hue_constant':
                res.append(((y[0],x),(y[1],x)))
            elif mode=='x_constant':
                res.append(((x,y[0]),(x,y[1])))
            else:
                raise Exception('hue_pair_generator mode error!')
    return res


'''
        beta_m=lambda x: np.log2(x/(1-x))
        m_beta=lambda x:1/(1+1/(2**x))
'''

def beta_m(x):
    return np.log2(x/(1-x))
def m_beta(x):
    return 1/(1+1/(2**x))


def m_choose_n(m,n,integer=False):     
    temp=lambda x:Decimal(math.factorial(x))

    with localcontext() as ctx:
        ctx.prec = 42
        s=temp(m)/temp(n)/temp(m-n)
    if integer:
        s=int(s)
    return s

'''
def intersection_pvalue(total,set1,set2,intersection):
    
def intersection_pvalue(total,set1,set2,intersection):
    getcontext().prec = 42
    temp=lambda x:Decimal(math.factorial(x))
    constan=temp(set1)*temp(set2)/temp(total)*temp(total-set1)*temp(total-set2)
    denominator=0
    for k in range(intersection):
        denominator+=Decimal(1)/(temp(k)*temp(set1-k)*temp(set2-k)*temp(total-set1-set2-k))
    return constan*denominator

    return res/
'''


def intersection_pvalue(total,set1,set2,intersection):
    getcontext().prec = 42
    constan=m_choose_n(total,set2)
    numerator=0
    for k in range(intersection):
        numerator+=m_choose_n(set1,k)*m_choose_n(total-set1,set2-k)
    return 1-numerator/constan


#intersection_pvalue(628,64,78,53)



#ranking difference
#set similarity