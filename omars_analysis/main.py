import numpy as np
import itertools
import datetime
import sys
from scipy.stats import t
from scipy.stats import f

from omars_analysis.subset_selection import get_soe
from omars_analysis.generate_model_matrix_heredity import create_model_matrix_heredity

def hat_matrix(smat):
    hat = smat @ np.linalg.pinv(smat.T @ smat) @ smat.T
    return hat

def create_quadratic_interactions(smat):
    no_runs, nf = smat.shape[0], smat.shape[1]
    quad_ix = [i for i,x in enumerate(smat.T) if (abs(x)!=1).any()]
    
    # Quadratic
    if len(quad_ix) > 0:
        quadratic = np.absolute(smat[:,quad_ix])
        soe_mat = np.hstack((np.ones((no_runs,1)),quadratic))
    else:
        soe_mat = np.ones((no_runs,1))

    # Interactions
    all_interactions = []
    for combi in itertools.combinations(range(nf), 2):
        col = (smat[:,combi[0]]*smat[:,combi[1]]).reshape((no_runs,1))
        soe_mat = np.hstack((soe_mat, col))
        all_interactions.append(str(combi[0]+1)+'_'+str(combi[1]+1))

    #centering (important for later functions)
    soe_mat = soe_mat - soe_mat.mean(axis=0)

    return soe_mat[:,1:], all_interactions, quad_ix

def code(smat):
    a = np.ones((smat.shape[0],1))
    for i in smat.T:
        maxi = max(i)
        mini = min(i)
        avg = (maxi + mini)/2
        diff = (maxi - mini)/2
        temp = (i - avg)/diff
        temp = temp.reshape((smat.shape[0],1))
        a = np.hstack((a,temp))
       
    return a[:,1:]

def get_omars_analysis(smat, cy, qheredity: str, iheredity: str, effects_to_drop: list, full='y'):

    # mat (numpy array) - design matrix
    # cy (numpy array) - response (single response, the response will be centered in the function)
    # qheredity (string) - 
    # iheredity (string) - 
    # full (string) - ('y', 'n'); n - only main effects selection
    # effects_to_drop (list) - list of second order effects to be excluded, eg: ['1_1', '2_3']

    mat = code(smat)
    no_runs, nf = mat.shape[0], mat.shape[1]

    # Centering response in case not centered already
    cy = (cy - cy.mean(axis=0)).reshape(no_runs,1)

    # Calculate error variance
    intercept = np.ones((no_runs, 1))
    mat_qi, all_interactions, quad_ix = create_quadratic_interactions(mat)
    total_mat = np.hstack((np.hstack((intercept, mat)),mat_qi))
    
    new_error_variance_projection = np.identity(no_runs) - (total_mat @ np.linalg.pinv(total_mat))
    denominator = no_runs - np.linalg.matrix_rank(total_mat)
    print('\nInitial error degrees of freedom - {}'.format(denominator))

    new_variance = (cy.T @ new_error_variance_projection @ cy)/denominator
    print('\nInitial estimate of error variance - {}'.format(float(new_variance)))
    
    s = np.sqrt(new_variance)
    
    '''
    ANALYSIS OF ME (STEP 1)
    '''
    mebeta = (np.linalg.inv(mat.T @ mat) @ mat.T @ cy).reshape(1,nf)
    vdiags = np.matrix.diagonal(np.linalg.inv(mat.T @ mat))
    svdiags = np.sqrt(vdiags).reshape(1,nf)
    tvalue = np.absolute(mebeta/(s*svdiags))
    ttest = t.ppf(1 - (0.05/2), df=denominator)
    tcomp = tvalue/ttest
    
    indexrme = [i for i,x in enumerate(tcomp[0]) if (abs(x)>=1)] # index of active main effects    
    indexfme = [i for i,x in enumerate(tcomp[0]) if (abs(x)<1)] # index of inactive main effects
    active = mat[:, indexrme]
    inactive = mat[:, indexfme]
    new_me = [str(g+1) for g in indexrme]
    print('\nActive main effects - {}'.format(new_me))

    if full == 'n':
        return None
    
    '''
    END OF STEP 1
    '''
    # Update error variance
    new_variance_projection_updated = new_error_variance_projection + hat_matrix(inactive)
    numerator = cy.T @ new_variance_projection_updated @ cy
    new_denominator = denominator + len(indexfme)
    new_variance = numerator/new_denominator
    print('\nUpdated estimate of error variance - {}'.format(float(new_variance)))
    
    '''
    ANALYSIS OF SOE (STEP 2)
    '''
    # Second order proection matrix (make sure to only multiply with centered response since intercept is ignored)
    y_second_projection = np.identity(no_runs) - hat_matrix(mat) - new_error_variance_projection
    cy_second = (y_second_projection @ cy).reshape(no_runs,1)
    TSS = cy_second.T @ cy_second # Total sum of squares (for second order space)
    
    ####        F Test for second order effects
    o_dfleft = np.linalg.matrix_rank(mat_qi) # quadratics need to be centered in mat_qi
    ftest = f.ppf(q= 1-0.2, dfn=o_dfleft, dfd=new_denominator)
    F1 = (TSS/o_dfleft)/new_variance
    fcomp = F1/ftest
    '''
    CREATE SECOND ORDER MATRIX WITH HEREDITY FOR ANALYSIS (WITHOUT INTERCEPT)
    '''
    if fcomp >= 1:
        # standardize
        std = 'n'
        soe_mat2, my_dictionary_soe_single  = create_model_matrix_heredity(std, mat, mat_qi, quad_ix, active, inactive, indexrme, indexfme, qheredity, iheredity)
        if len(effects_to_drop) > 0:
            new_dict = {}
            reduced_matrix_indices = []
            ittr = 0
            for keyz in my_dictionary_soe_single.items():
                if keyz[0] in effects_to_drop:
                    continue
                else:
                    reduced_matrix_indices.append(keyz[1])
                    new_dict.update({keyz[0]:ittr})
                    ittr += 1
            soe_mat2 = soe_mat2[:,reduced_matrix_indices]
            my_dictionary_soe_single = new_dict
            
        
        # if quadratic_special_coding == 'y':
        #     relevant_q = []
        #     for dumzy in my_dictionary_soe_single:
        #         if dumzy.split('_')[0] == dumzy.split('_')[1]:
        #             relevant_q.append(int(dumzy.split('_')[0])-1)
                    
        #     if len(relevant_q)>0:
        #         temp_mat = mat[:,relevant_q]
        #         temp_mat = temp_mat - temp_mat.mean(axis=0) # center
        #         temp_mat = np.square(temp_mat) # square
        #         temp_mat = code(temp_mat) # normalize
        #         soe_mat2[:,0:len(relevant_q)] = temp_mat.copy()

        i_start = datetime.datetime.now()
        active_soe, p_value_omars = get_soe(mat, o_dfleft, soe_mat2, cy_second, new_denominator, new_variance, my_dictionary_soe_single)
        i_finish = datetime.datetime.now()
        timing = str((i_finish - i_start).total_seconds())
    else:
        active_soe = []
        p_value_omars = np.nan
    
    new_qe = []
    new_ie = []
    for d in active_soe:
        x = d.split('_')
        if x[0] == x[1]:
            new_qe.append(d)
        else:
            new_ie.append(d)

    print('\nActive interaction effects - {}'.format(new_ie))
    print('\nActive quadratic effects - {}'.format(new_qe))

    return p_value_omars