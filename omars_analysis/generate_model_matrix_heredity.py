import numpy as np
import itertools

'''
Inputs
qheredity = y,n
iheredity = s,w,n

'''

def create_model_matrix_heredity(std, mat, quad_ix, indexrme, indexfme, qheredity, iheredity):
    # my_dictionary_soe_ci - individual column indices for each soe column
    # my_dictionary_soe_single - location of a certain soe column in the soe matrix
    
    no_runs, nf = mat.shape[0], mat.shape[1]
    active = mat[:,indexrme]

    my_dictionary_soe_ci = {}
    my_dictionary_soe_single = {}
    if active.shape[1] == 0:# or active.shape[1] == 1:
        qheredity = 'n'
        iheredity = 'n'
    
    subtract = 2 # default
    if len(quad_ix) > 0: # CORRECT CONDITION and check if in indexrme
        if qheredity == 'y':
            # print('\nHeredity for quadratic effects assumed')
            ittr = 0
            for k,j in enumerate(active.T):
                if np.any(abs(j) != 1):
                    v_key = str(indexrme[k]+1)+"_"+str(indexrme[k]+1)
                    # current_dictionary = {v_key: [indexrme[k], indexrme[k]]}
                    # my_dictionary_soe_ci.update(current_dictionary)
                    current_dictionary2 = {v_key: ittr}
                    my_dictionary_soe_single.update(current_dictionary2)
                    ittr += 1
            if ittr == 0: # no three level active
                soe_mat2=np.ones((no_runs,1))
            else:
                subtract = 1                
                active_q = list(set(indexrme).intersection(quad_ix))
                soe_mat2 = np.absolute(mat[:,active_q])
        else:
            # print('\nNo Heredity assumed for quadratic effects')
            subtract = 1
            soe_mat2 = np.absolute(mat[:,quad_ix])
            for j in range(soe_mat2.shape[1]):
                v_key = str(quad_ix[j]+1)+"_"+str(quad_ix[j]+1)
                current_dictionary = {v_key: [quad_ix[j], quad_ix[j]]}
                my_dictionary_soe_ci.update(current_dictionary)
                current_dictionary2 = {v_key: j}
                my_dictionary_soe_single.update(current_dictionary2)
    else:
        soe_mat2=np.ones((no_runs,1))
        
    # Pattern - conventional
    interactions = []
    if iheredity == 's':
        # print('\nStrong Heredity assumed for interactions')    
        for combi in itertools.combinations(range(active.shape[1]), 2):
            iron = (active[:,combi[0]] * active[:,combi[1]]).reshape((no_runs,1))    
            soe_mat2 = np.hstack((soe_mat2, iron))
            v_key = str(indexrme[combi[0]]+1)+"_"+str(indexrme[combi[1]]+1)
            interactions.append(v_key)
            # current_dictionary = {v_key: [indexrme[combi[0]], indexrme[combi[1]]]}
            # my_dictionary_soe_ci.update(current_dictionary)
            current_dictionary2 = {v_key: soe_mat2.shape[1]-subtract} # subtract since in this loop a column is added
            my_dictionary_soe_single.update(current_dictionary2)
    # Pattern - Corrected
    elif iheredity == 'w':
        # index = indexrme + indexfme
        # print('\nWeak Heredity assumed for interactions')
        for combi in itertools.combinations(range(active.shape[1]), 2):
            iron = (active[:,combi[0]] * active[:,combi[1]]).reshape((no_runs,1))    
            soe_mat2 = np.hstack((soe_mat2, iron))
            v_key = str(indexrme[combi[0]]+1)+"_"+str(indexrme[combi[1]]+1)
            interactions.append(v_key)
            # current_dictionary = {v_key: [indexrme[combi[0]], indexrme[combi[1]]]}
            # my_dictionary_soe_ci.update(current_dictionary)
            current_dictionary2 = {v_key: soe_mat2.shape[1]-subtract}
            my_dictionary_soe_single.update(current_dictionary2)
        for combi2 in range(active.shape[1]):
            for combi3 in range(len(indexfme)):
                o = (active[:, [combi2]] * mat[:, [indexfme[combi3]]]).reshape((no_runs,1))
                soe_mat2 = np.hstack((soe_mat2, o))
                mini, maxi = min(indexrme[combi2], indexfme[combi3]), max(indexrme[combi2], indexfme[combi3])
                v_key = str(mini+1)+"_"+str(maxi+1)
                interactions.append(v_key)
                # current_dictionary = {v_key: [mini, maxi]}                   
                # my_dictionary_soe_ci.update(current_dictionary)
                current_dictionary2 = {v_key: soe_mat2.shape[1]-subtract}
                my_dictionary_soe_single.update(current_dictionary2)
    # Pattern - conventional
    else:
        # print('\nNo heredity assumed for interactions')
        for combi in itertools.combinations(range(nf), 2):
            iron = (mat[:,combi[0]] * mat[:,combi[1]]).reshape((no_runs,1))   
            soe_mat2 = np.hstack((soe_mat2, iron))
            v_key = str(combi[0]+1)+"_"+str(combi[1]+1)
            interactions.append(v_key)
            # current_dictionary = {v_key: [combi[0], combi[1]]}
            # my_dictionary_soe_ci.update(current_dictionary)
            current_dictionary2 = {v_key: soe_mat2.shape[1]-subtract}
            my_dictionary_soe_single.update(current_dictionary2)
    
    if std == 'y':
        soe_mat2 = (soe_mat2 - np.mean(soe_mat2, axis=0)) / np.std(soe_mat2, axis=0)
    else:
        soe_mat2 = soe_mat2 - np.mean(soe_mat2, axis=0)
    
    if subtract == 2:
        soe_mat2 = soe_mat2[:,1:]
        
    return soe_mat2, my_dictionary_soe_single
