import numpy as np
import itertools
import scipy
from scipy.stats import f

def get_key(val, x): 
    for key, value in x.items(): 
         if val == value: 
             return key   

    return "key doesn't exist"


def get_soe(mat, limit_for_soe_subset, soe_mat2, cy_second, denominator, new_variance, my_dictionary_soe_single, alpha_val=0.2, limit_for_step_two=None):  
    nf = mat.shape[1]
    total_soe_terms = nf*(nf+1)/2 
    
    # List of active SOE terms
    active_soe = []

    dfleft = limit_for_soe_subset # rank of all second order terms matrix

    actual_soe_df = np.linalg.matrix_rank(soe_mat2)

    if limit_for_step_two != None: # user specified
        if limit_for_step_two <= actual_soe_df:
            dummy5 = limit_for_step_two
        else:
            dummy5 = actual_soe_df
        concluding_statement = '\nLimit on the number of terms for subset selection - {} (user specified)'.format(limit_for_step_two)
    else:
        if total_soe_terms == limit_for_soe_subset:
            dummy5 = limit_for_soe_subset
            concluding_statement = '\nNo limit is set on the number of terms for subset selection since the design can estimate all second order terms'.format(dummy5)
        else:
            N_by_4 = int(mat.shape[0]/4)
            if soe_mat2.shape[1] <= limit_for_soe_subset or actual_soe_df <= N_by_4:
                dummy5 = actual_soe_df
                concluding_statement ='\nLimit on terms for subset selection - {} (Maximum number of second order terms jointly estimable of all the second order effects considered)'.format(dummy5)
            else:
                dummy5 = N_by_4
                concluding_statement = '\nLimit on terms for subset selection - {} (run size divided by 4)'.format(dummy5)
    
    # Skipped 'intercept' as response is centered
    
    ###################################################################
    ############    CODE FOR SUBSET SELECTION  #######################
    ###################################################################
    # m = subset size in consideration
    # id = last combination with lowest RSS
    # soe_mat2 does not have the intercept
    # dummy5 = # of df left for second order terms
    m = 0
    for n in range(dummy5):
        m += 1
        first = 0
        for combi in itertools.combinations(range(soe_mat2.shape[1]), m):
            dummy7 = np.ones((soe_mat2.shape[0],1))
            if m == 1:
                dummy8 = soe_mat2[:,combi[0]]
                dummy8.shape = (soe_mat2.shape[0],1)
                dummy7 = np.hstack((dummy7, dummy8))
                v_key = combi[0]
                resid = np.linalg.lstsq(dummy7, cy_second, rcond=None)[1]
                if combi[0] == 0:
                    l_rss = resid
                    l_v_key = v_key                
            else:
                dummy8 = soe_mat2[:,combi]
                dummy7 = np.hstack((dummy7, dummy8))
                resid = np.linalg.lstsq(dummy7, cy_second, rcond=None)[1]
                if resid.size == 0:
                    continue
                v_key = combi
                if first == 0:
                    l_rss = resid
                    l_v_key = v_key
                    first = 1
            if resid < l_rss:
                l_rss = resid
                l_v_key = v_key
        dfleft -= 1
        ftest = f.ppf(q= 1-alpha_val, dfn=(dfleft), dfd=denominator)
        F2 = (l_rss/(dfleft))/new_variance
        fcomp = F2/ftest
        fcomp = float(fcomp)
        if m == soe_mat2.shape[1] or fcomp < 1:
            break

    p_value_omars = 1-scipy.stats.f.cdf(fcomp, dfleft, denominator)
    
    if m == 1:
        active_soe.append(get_key(l_v_key, my_dictionary_soe_single))
    else:
        for s in range(len(l_v_key)):
            active_soe.append(get_key(l_v_key[s], my_dictionary_soe_single))

    return active_soe, p_value_omars, actual_soe_df, concluding_statement