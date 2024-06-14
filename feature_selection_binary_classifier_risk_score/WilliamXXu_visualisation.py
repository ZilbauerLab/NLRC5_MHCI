#probably no the most artistic


from WilliamXXu_utilities import pair_generator, hue_pair_generator
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind,bartlett, f_oneway
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
#from pingouin import welch_anova
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

'''if with_control:
        joined_processed=joined
    else:
        joined_processed=data_static.only_CD(joined)'''
def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)
def box_plot(x,entity,mode='function',figure_size=(12, 6),font_scale=1.2,title_font=20,dot_plot=True,order=None,fig_format='png',significance_level=True,name_plot=None,y_label='Beta Value',x_label=' ',pval_format='star',linear=False,colors=None,pval_pairs=None,multiple_pval_location=None,group_by='signature',significance_test='t',y_range=[0,1],show_plot=True,grid_line=True):
    '''
    make box plot, with optional dot plot and significance test. Plot saved under foler fig

    parameters:
    x - input of pandas DataFrame, raw numerical information (eg beta values). Contain a column named 'target' as the grouping information for each sample (eg. CD severity).
    entity - single mode: column id (eg. cpg id). function mode: a function, which turns the input matrix x into a column of numeric results, and also returns a name as the second return value
            multiple functions: a list of functions as described in function mode
    mode - 3 options: single, function, multiple fuctions. single plots value of a single column for the y-axis. function plots the value of a function applied to columns for the y-axis (eg. average or linear combination).
            multiple functions plots a list of functions as described in function mode, with hue used to distinguish functions/target groups
    dot_plot - whether to plot a dot plot or not
    fig_format: format of the figure, eg. 'png' or 'pdf' for matplotlib.pyplot savfig
    significance_level: include significance level bar or not
    name_plot: name of the plot. If not specified, will autogenerate the name
    y_label: label for the y-axis
    pval_format: format of the p-value, 'full' or 'simple' for statannotations.Annotator
    linear: y-axis linear scale or not
    pval_pairs: a list of pairs of p-values to be included in the significance test.  If not specified, will autogenerate all pairs
    multiple_pval_location: required for multiple comparison significance test (eg. anova). location of the p-value bar
    group_by: required for multiple functions mode. Grouping x-axis by 'signature' (calculated by different functions) or 'target' (like clinical phenotypes passed from input x)
    significance_test: name of significance test. Function mode:  t_welch (unequal variance), t_m_welch (welch t test on m values), t, mannwhitneyu, bartlett, anova. Multiple functions mode: t_welch, manwhitneyu
    y_range: y-axis range
    grid_line: whether to plot a grid line or not
    show_plot: if true, will show the plot
    '''
    
    
    #'single','function'
    #
    #print(x)
    
    
    pairs=pval_pairs
    beta_m=lambda x: np.log2(np.divide(x,(1-x)))


    '''
    ind_common=list(set(x.index).intersection(y.index))
    y=y.loc[ind_common,:]
    x=x.loc[ind_common,:]
    '''


    joined_processed=x
    print('n={}'.format(joined_processed.shape[0]))
    #print(joined_processed)
    #print(joined_processed.target)
    #print(joined_processed.target)
    def get_log_ax(figsize=(12, 6)):
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        fig.patch.set_alpha(1)
        #getattr(ax, set_scale)("log")
        return ax 
    
    if mode=='single':
        cpg=entity
        print(cpg)
        signature=joined_processed[cpg]
        plotting_parameters = {'data':joined_processed,'x':"target",'y':cpg}
        default_name=cpg#+' CpG Beta Value'
    elif mode=='function':
        print('start')
        signature,default_name=entity(joined_processed.drop(['target'],axis=1))
        joined_processed['signature']=signature
        
        plotting_parameters = {'data':joined_processed,'x':"target",'y':'signature'}
        print('ready')
    elif mode =='multiple functions':
        data_concat=[]
        all_names=[]
        for k in entity:
            temp=joined_processed.copy()
            signature,name=k(joined_processed.drop(['target'],axis=1))
            temp['signature']=signature
            temp['signature_type']=name
            all_names.append(name)
            data_concat.append(temp)
        default_name=', '.join(all_names)
        joined_processed=pd.concat(data_concat)
        if group_by=='target':
            plotting_parameters = {'data':joined_processed,'x':"target",'y':'signature','hue':'signature_type'}
            if pairs is not None:
            
            #print(pairs)
            #print(all_names)
                pairs=hue_pair_generator(pairs,all_names,mode='hue_constant')
            #print(pairs)
        else:
            plotting_parameters = {'data':joined_processed,'x':'signature_type','y':'signature','hue':"target"}   
            if pairs is not None:
                #print(pairs)
                #print(all_names)
                    pairs=hue_pair_generator(pairs,all_names,mode='x_constant')
                #print(pairs)
        
        #
        # 
        
    else:
        raise Exception('mode unsupported')
          
    sns.set(style="whitegrid")
    with sns.plotting_context("notebook", font_scale =font_scale):
    # Create new plot
        ax = get_log_ax(figsize=figure_size)

        '''
                if color is not None:
            sns.set_palette(sns.color_palette(color))
        '''

    # Plot with seaborn
        box=sns.boxplot(ax=ax, **plotting_parameters,showfliers = False,order=order,palette=colors)
        #box.add_legend(bbox_to_anchor=(1.05, 0), loc=2, borderaxespad=0.)
        #box.ax.legend(loc=2)
        #sns.move_legend(box, "upper left")
        if dot_plot ==True:

            '''
                     if color is not None:
                sns.set_palette(sns.color_palette(color))   
            '''

            swarm=sns.swarmplot(ax=ax, **plotting_parameters,dodge=True,order=order,palette=colors,edgecolor="black",linewidth=1.0)#,color=color)#,color=".25")
        #sns.move_legend(swarm, "upper left")
        #swarm.ax.legend(loc=2)
        ax.set(xlabel=x_label, ylabel=y_label)
        if not grid_line:
            ax.grid(False)
        
        if linear==True:
            box.set_yticks(y_range)
            if dot_plot ==True:
                swarm.set_yticks(y_range)
        if significance_level:
            print('start again')
            groups=list(set(joined_processed.target))
            #print(joined_processed.target)
            #print(groups)
            if (mode == 'multiple functions') and (pairs is None):
                pass
            else:
                pvals=[]
                if mode != 'multiple functions':
                    if multiple_pval_location is None:
                        if pairs is None:
                            pairs=pair_generator(groups,'no rep')
                        

                        for x in pairs:
                            #temp=list(x)
                            #test_res=mannwhitneyu(beta_values[joined_processed.ploting==temp[0]],beta_values[joined_processed.ploting==temp[1]])
                            one=signature[joined_processed.target==x[0]]
                            two=signature[joined_processed.target==x[1]]
                            if significance_test =='t_m_welch':
                                test_res=ttest_ind(beta_m(one),beta_m(two),equal_var=False)
                            elif significance_test =='t_welch':
                                test_res=ttest_ind(one,two,equal_var=False)
                            elif significance_test =='t':
                                test_res=ttest_ind(one,two,equal_var=True)
                            elif significance_test =='mannwhitneyu':
                                test_res=mannwhitneyu(one,two)
                            else:
                                raise Exception('significance test not supported')
                            pvals.append(test_res.pvalue)
                            print(test_res.pvalue)

                    else:
                        if pairs is None:
                            pairs=[groups]
                        
                        for x in pairs:
                            multi=[]
                            for k in x:
                                multi.append(list(signature[joined_processed.target==k]))
                            if significance_test=='bartlett':
                                test_res=bartlett(*multi)#,nan_policy='raise')
                            #print(multi)
                            #print(test_res)
                            elif significance_test=='anova':
                                test_res=f_oneway(*multi)
                            else:
                                raise Exception('test not supported!')
                        #print(test_res[0])
                            pvals.append(test_res.pvalue)  
                else:
                    for x in pairs:
                        if group_by=='target':
                        #test_res=mannwhitneyu(beta_values[joined_processed.ploting==temp[0]],beta_values[joined_processed.ploting==temp[1]])
                            double_select=lambda y: joined_processed.signature[(joined_processed.target==y[0]) & (joined_processed.signature_type==y[1])]

                        else:
                            double_select=lambda y: joined_processed.signature[(joined_processed.target==y[1]) & (joined_processed.signature_type==y[0])]
                        if significance_test =='t_welch':
                            test_res=ttest_ind(double_select(x[0]),double_select(x[1]),equal_var=False)
                        elif significance_test =='mannwhitneyu':
                            test_res=mannwhitneyu(double_select(x[0]),double_select(x[1]))
                        else:
                            raise Exception('test not supported!')                            
                        pvals.append(test_res.pvalue)                      
                # Add annotations
                print(pvals)
                print(pairs)
                if multiple_pval_location is not None:
                    pairs=multiple_pval_location
                annotator = Annotator(ax, pairs, order=order,**plotting_parameters)
                annotator.configure(text_format=pval_format)#,test_short_name=' ')
                #annotator.configure(text_format="full",test_short_name="2-sided Mann–Whitney U test")
                annotator.set_pvalues_and_annotate(pvals)

        '''
        ticks=0.05*np.array(range(10))
        ticklabels=[str(x) for x in ticks]
        ax.set_yticks(ticks)        
        '''
        if (linear==True) or (significance_level==False):
            plt.yscale('linear')

        if name_plot is None:
            name_plot=default_name+'. Significiance: '+str(significance_level)+'. Linear Scale: '+str(linear)+'. Significance Test: '+significance_test
        plt.title(label=name_plot,fontsize=title_font)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        '''
        def handle_legend(plot_method):
            plt.setp(plot_method.get_legend().get_texts(), fontsize='5')  
            #plt.setp(plot_method.get_legend().get_title(), fontsize='8') 
            #plot_method.legend(bbox_to_anchor= (100,1))
        handle_legend(box)
        handle_legend(swarm)
        '''
        if mode=='multiple functions':
            
            
            plt.legend(fontsize='6', title_fontsize='6')
            plt.legend(loc='upper left')
            
        plt.savefig('fig/'+name_plot+'.'+fig_format,format=fig_format)
                
                
        if show_plot==True:
            plt.show()
        plt.close()








def dimension_reduction(x,y,reduction_method='PCA',n_components=2,reducer_parameters=dict(),dimensions_to_plot=(0,1),plot=True,name_plot=None,plot_only_CD=True,legend_location='upper left',legend=None,show_plot=True,fig_format='png',class_labels=None):
    #produce dimensionality reduction, with the potential to plot
    
    #parameters:
    # x - input of pandas DataFrame, numeric values for dimensionality reduction (eg beta values). Contain the last column as the grouping information for each sample (eg. CD severity).
    # y - input of pandas DataFrame, clinical values as produced from data_processing.py. Contains a column of diagnosis_number with 1 representing Crohn disease
    # reduction_method - the method for dimensionality reduction. PCA, TSNE, UMAP
    # n_components - number of components after dimensionality reduction. only 2 can be plotted for now
    # plot - if true will plot and save under the foler fig
    # name_plot: name of the plot. If not specified, will autogenerate the name
    # plot_only_CD: to plot only CD samples or not
    # legend_location: location of the legend, eg. upper left
    # legend: legend description for the plot
    # show_plot: if true, will show the plot
    # fig_format: format of the figure, eg. 'png' or 'pdf' for matplotlib.pyplot savfig

    #returns:
    # embedding - same rows as x and y, but with n_components columns of embedding on the reduced dimensions
    # reducer - fitted dimesionality reducer, can be used for future projection. Check original packages

    if legend is None:
        if plot_only_CD:
            legend='Severity: 0 is severe'
        else:
            legend='Severity: 0 is TI, 1 is SC, 2 is DUO'

    if name_plot is None:
        name_plot=reduction_method+' Projection of Beta Values from '+str(x.shape[1]-1)+' CpGs. Only CD: '+str(plot_only_CD)
    if reduction_method=='PCA':
        reducer=PCA(n_components=n_components,**reducer_parameters)
    elif reduction_method=='UMAP':
        reducer = umap.UMAP(n_components=n_components,**reducer_parameters)
    elif reduction_method=='TSNE':
        reducer = TSNE(n_components=n_components,**reducer_parameters)
    else:
        raise ValueError('reducer_name not supported')
    embedding = reducer.fit_transform(x.iloc[:,0:-1])


    if plot==True:
        if plot_only_CD:
            switch=y['diagnosis_number']==1
        else:
            switch=y['diagnosis_number']>-10000
        fig, ax = plt.subplots()
        scatter=ax.scatter(
            embedding[switch, dimensions_to_plot[0]],
            embedding[switch, dimensions_to_plot[1]],
            #c=[sns.color_palette()[x] for x in y.prognosis_binary])
            c=x.iloc[:,-1].loc[switch])
        plt.gca().set_aspect('equal', 'datalim')

        if name_plot is None:
            name_plot=reduction_method+' Projection of Beta Values from '+str(x.shape[1]-1)+' CpGs. Only CD: '+str(plot_only_CD)
        plt.title(name_plot)#, fontsize=24)
        '''legend1 = ax.legend(*scatter.legend_elements(),
                            loc="upper left", title="Severity: 0 is severe")
                            '''


        if class_labels is None:
            legend1 = ax.legend(handles=scatter.legend_elements()[0],loc=legend_location,title=legend)#title="0 is TI, 1 is SC, 2 is DUO")
        else:
            legend1 = ax.legend(handles=scatter.legend_elements()[0],loc=legend_location,title=legend,labels=class_labels)
        ax.add_artist(legend1)
        plt.savefig('fig/'+name_plot+'.'+fig_format,format=fig_format)
        if show_plot==True:
            plt.show()
        plt.close()
        
        
    return embedding,reducer