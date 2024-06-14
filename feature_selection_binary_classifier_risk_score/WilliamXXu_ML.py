from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression,ElasticNet
from sklearn.svm import LinearSVC,SVC
from sklearn.neighbors import KNeighborsClassifier
#from lightgbm import LGBMClassifier
from sklearn.naive_bayes import GaussianNB
from catboost import CatBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
import pandas as pd
import numpy as np
import umap
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score, f1_score,recall_score,confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve
from scipy.stats import iqr,levene,fligner,ks_2samp,anderson_ksamp,mannwhitneyu, boxcox

from sklearn.preprocessing import PowerTransformer, StandardScaler

#from autogluon.tabular import TabularDataset, TabularPredictor

def choose_preprocessor(input):
    if input[0]=='standardise':
        return StandardScaler()
    elif input[0]=='power':
        return PowerTransformer()
    elif input[0][0:10]=='univariate':
        return univariate_filter(input[0][11:],None,input[1])
    else:
        raise Exception('Preprocessing method not supported or typo.')
    
def choose_feature_engineering(input):

    if input[0] =='umap':
        return umap.UMAP(n_components=input[1])
        
    elif input[0] == 'kmeans':
        return KMeans(n_clusters=input[1])
    else:
        raise Exception('Feature engineering method not supported or typo.')

def choose_model(model_name,para):
    if model_name =='xgboost':
        return XGBClassifier(n_estimators=100,scale_pos_weight =para)
    if model_name=='NaiveBayes':
        return GaussianNB()
    if model_name=='knn':
        return KNeighborsClassifier(n_neighbors=para)
    elif model_name=='svm':
        #return LinearSVC(C=1,class_weight={1:para,0:1},dual='auto')
        #{‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’} 
        return SVC(C=1,kernel=para,probability=True)
    elif model_name=='logistic':
        return LogisticRegression(C=1,class_weight={1:para,0:1},max_iter=3000)
    elif model_name=='elasticNet':
        return LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=para[0],C=para[1]) #LogisticRegression penalty='elasticnet',solver='saga'
    elif model_name=='randomForest':
        return RandomForestClassifier(n_estimators=100,class_weight={1:para,0:1})
    elif model_name=='neuroMLP':
        return MLPClassifier(random_state=100, max_iter=1000, hidden_layer_sizes=(50, 50))
    elif model_name=='autoGluon':
        pass
        #return TabularPredictor(label='prognosis_binary',eval_metric='roc_auc',verbosity=0,log_to_file=False)
    elif model_name =='lightgbm':
        pass
        #return LGBMClassifier(reg_alpha=0,reg_lambda=para)
    elif model_name =='catboost':   
        return CatBoostClassifier(verbose=0)#reg_lambda=para
    else:
        raise Exception('Model not supported or typo.')



def choose_measure(model_name,meas,model,x,y,threshold):
    '''
    print(y_test_processed)
    print(model.predict(x_test_processed))
    print(x_test_processed)        
'''
    if model_name=='autoGluon':
        
        evaluate_data=TabularDataset(pd.concat([x,y.rename('prognosis_binary')],axis=1))
        res=model.evaluate(evaluate_data,silent=True,detailed_report=True)
        #print(res)
        if meas=='roc' or meas=='specificty' or meas=='cross_entropy':
            raise Exception('Measure under development')
        else:
            return res[meas]
    else:
        y_proba=pd.DataFrame(model.predict_proba(x), columns=model.classes_)[1]
        y_predict=(y_proba>threshold).map({True:1,False:0})

        if meas=='confusion_matrix':
                return confusion_matrix(y,y_predict)

        elif meas=='f1':
                return f1_score(y,y_predict)

        elif meas =='specificity':
                return recall_score(1-y,1-y_predict)

        elif meas=='recall':
                return recall_score(y,y_predict)
        elif meas=='accuracy':
                return accuracy_score(y,y_predict)

                '''
                if meas == 'continuous':
                    def proba_measure(prediction,actual):
                        
                        #prediction=0.5*np.ones((actual.shape[0],2))
                        diff=actual_transformed-prediction
                        return np.mean(np.array([ np.dot(j,j)/2 for j in diff]))
                    score=proba_measure(model.predict_proba(x_test.loc[cd,:]),y_test['binary'].loc[cd])
                    score_training=proba_measure(model.predict_proba(x_train),y_train['binary'])

                '''
        elif meas=='roc':
            mean_fpr=np.linspace(0,1,100)
            fpr, tpr, thresholds = roc_curve(y,y_proba)

            #viz=RocCurveDisplay.from_estimator(model,x,y)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            return interp_tpr
            
            
        elif meas == 'roc_auc':
            mean_fpr=np.linspace(0,1,100)
            return roc_auc_score(y,y_proba)

        elif meas == 'cross_entropy':
                    def cross_entropy_binary(prediction,actual):

                        #prediction=0.5*np.ones((actual.shape[0],2))
                        actual=[[1-k,k] for k in actual]
                        delta=1e-7
                        return np.mean(-np.sum(actual*np.log2(prediction+delta),axis=1))

                    y_hat=model.predict_proba(x)
                    return cross_entropy_binary(y_hat,y)
        else:
            raise Exception('Measure not supported or typo.')


def threshold_determination(method,model,x,y,cutoff=None):
    proba=pd.DataFrame(model.predict_proba(x), columns=model.classes_)[1]
    fpr, tpr, thresholds_roc = roc_curve(y,proba)
    precision, recall, thresholds_pr =precision_recall_curve(y,proba)

    if method == 'specificity':
        for i in range(len(fpr)):
            if fpr[i]<cutoff:
                res=thresholds_roc[i]
        return res
    elif method == 'g_means':
        gmeans = tpr * (1-fpr)
        ix = np.argmax(gmeans)
        return thresholds_roc[ix]
    elif method =='j_index':
        J = tpr - fpr
        ix = np.argmax(J)
        return thresholds_roc[ix]
    elif method =='top_left_index':
        tl=(1 - tpr) ** 2 + fpr ** 2
        ix = np.argmin(tl)
        return thresholds_roc[ix]
    elif method =='f1':
        
        f1=(2 * precision * recall) / (precision + recall+np.finfo(float).eps)
        ix = np.argmax(f1)
        return thresholds_pr[ix]
    elif method =='accuracy':
        accuracy_scores = []
        for thresh in thresholds_roc:
            accuracy_scores.append(accuracy_score(y, [m > thresh for m in proba]))
        accuracies = np.array(accuracy_scores)
        #max_accuracy = accuracies.max() 
        return thresholds_roc[accuracies.argmax()]


def univariate_statistics(x,y,test_names,labels=[1,0],test_parameters=None,sort_by=None):
    '''
        print(x.shape)
    print(y.shape)
    print(labels)
    print(labels[1])
    print(y==labels[1])
    print(y==labels[0])
    '''

    x=x.loc[y.index,:]
    sev=x.loc[y==labels[0],:]
    mm=x.loc[y==labels[1],:]
    
    #print(sev)
    #print(mm)
    res={}
    for test_name in test_names:
        res[test_name]=[]
    count=0
    for cpg in x.columns:
        #print(count)
        count+=1
        for k in range(len(test_names)):
            test_name=test_names[k]
            if test_name=='iqr':
                temp=iqr(sev[cpg])/iqr(mm[cpg])
            elif test_name=='levene':
                temp=levene(sev[cpg],mm[cpg],center=test_parameters[k]).pvalue
            elif test_name=='fligner':
                temp=fligner(sev[cpg],mm[cpg],center=test_parameters[k]).pvalue
            elif test_name=='anderson_ksamp':
                temp=anderson_ksamp([sev[cpg],mm[cpg]]).pvalue
            elif test_name=='ks_2samp':
                temp=ks_2samp(sev[cpg],mm[cpg]).pvalue
            elif test_name=='severe_var_smaller':
                temp=np.var(sev[cpg])<np.var(mm[cpg])
            elif test_name=='mannwhitneyu':
                temp=mannwhitneyu(sev[cpg],mm[cpg]).pvalue
            res[test_name].append(temp)

    df=pd.DataFrame.from_dict(res)
    df=df.set_index(x.columns)
    if not (sort_by is None):
        df=df.sort_values(by=sort_by)
    return df



class univariate_filter():
    def __init__(self,filter_name,filter_param,cut_off):
        self.filter_name = filter_name
        self.filter_param = filter_param
        self.cut_off = cut_off
        '''
        if purpose == 'Diagnosis':
            self.labels=['CD','Control']
        elif purpose == 'Prognosis':
            self.labels=['Severe CD','Mild/Moderate CD']
        else:
            raise Exception('purpose not supported or typo.')        
        
        '''
    def fit_transform(self,x,y,real):
        #joined=pd.concat([x,y.rename('target')],axis=1)
        
        res=univariate_statistics(x.loc[real],y.loc[real],[self.filter_name],test_parameters=[self.filter_param],sort_by=self.filter_name)
        #print(res.tail(20))
        #print(self.cut_off)
        #print(type(res))
        #print(res.iloc[:,0])
        temp=res[res.iloc[:,0]<self.cut_off]
        self.ind_keep=temp.index
        #print(self.ind_keep)
        output=x.loc[:,self.ind_keep]
        #print(output.shape)
        return output
    def transform(self,x):
        return x.loc[:,self.ind_keep]