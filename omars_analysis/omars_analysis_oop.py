import numpy as np
import pandas as pd
from omars_analysis.main import get_omars_analysis

class om_analysis:

    def __init__(self, smat_i: np.ndarray, cy_i: np.ndarray, alpha_i: list = [0.05, 0.2, 0.2], qheredity_i: str='n', iheredity_i: str='n', effects_to_drop_i: list=[], full_i: str='y', force_me_i: list=[], user_limit_for_step_two_i=None):
        
        list_of_outputs = get_omars_analysis(smat = smat_i, cy = cy_i, alpha = alpha_i, qheredity = qheredity_i, 
                                             iheredity = iheredity_i, effects_to_drop = effects_to_drop_i, full = full_i, force_me = force_me_i, user_limit_for_step_two=user_limit_for_step_two_i)

        self.collected_output = list_of_outputs
        self.design = smat_i

        if list_of_outputs['success'] == False:
            print('Analysis cannot continue due to lack of available degrees of freedom to estimate error variance')
            # think about app
        else:
            # assign to self
            self.intial_error_df = list_of_outputs['initial_df']
            self.intial_rmse = list_of_outputs['initial_rmse']
            self.active_me = list_of_outputs['active_me']

            if full_i == 'y':
                self.updated_rmse = list_of_outputs['updated_rmse']
                if list_of_outputs['active_soe']:
                    self.active_ie = list_of_outputs['active_ie']
                    self.active_qe = list_of_outputs['active_qe']
                    self.p_value_f_test_first = list_of_outputs['p_value_initial_f_test']
                    self.p_value_f_test_final = list_of_outputs['p_value_final_f_test']
                    self.analysis_time = list_of_outputs['analysis_time']
                else:
                    print('No active second order effects detected')

    def print_ME_p_values(self):
        df = pd.DataFrame (self.collected_output['me_p_values'], columns = ['p-value (threshold alpha - 0.05)'], index=[itr+1 for itr in range(self.design.shape[1])])
        return df

    def print_rank_statements(self):
        try:
            statements = self.collected_output['rank_statements']
            for itr in statements:
                print(itr)
        except:
            print('No rank statements available')
        return None
