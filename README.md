# Documentation for package omars-analysis

This package has one function 'get_omars_analysis'

First run the following:

```python
from omars_analysis.main import get_omars_analysis
```

Running the above command will make the function 'get_omars_analysis' available for use.

## Function usage

The function can be used with nine inputs (only the first two are required):

- smat
  
  the first input is the design matrix. This should be a numpy array. The design matrix need not be coded, however it must only consist of continuous variables. The function is meant to be used with designs that have all factors either with two levels or three levels per factor column. The design should **NOT** consist of headers. The design matrix should not consist of second order effects since this will be built internally in the function.

- cy
  
  This is the response. This should be a one dimensional column vector (1-D numpy).

- qheredity
  
  This is to specify heredity constraints for quadratic effects. The accepted inputs are 'y' or 'n' ('y'- strong heredity, 'n'- no heredity, 'n'- No heredity). The input must be a string in lowercase. The default is 'n'.

- iheredity
  
  This is to specify heredity constraints for two-factor interaction effects. The accepted inputs are 's', 'w' or 'n' ('s'- strong heredity, 'w'- weak heredity, 'n'- no heredity). The input must be a string in lowercase. The default is 'n'.

- alpha
  Here the three different alpha values can be specified (refer paper for mor information). The input must be a list of alpha values in float format. The default value for this parameter is [0.05, 0.2, 0.2].

- effects_to_drop
  
  This is to specify second order effects that must be excluded from the analysis. The input must be a list of strings. For example: ['1_1', '2_3']. This input specifies that the quadratic effect of the first factor and the interaction effect between factor two and three must be excluded from the second step of the analysis (subset selection). The default value for this parameter is an empty list.

  Note: Entering interaction between factor one and two should be represented as '1_2' and not as '2_1'. The smaller factor number should always come first.

- full
  
  'n' -  analysis is performed on the main effects only
  
  'y' - analysis is performed on the main effects and second order effects.

  The default is set to 'y'

- force_me
  
  Here certain main effects can be forced into the model. This can be used in cases where a main effect is statistically insignificant but by only a small margin. Forcing such main effects into the model can have an impact in reducing the updated estimate of the error.
  
  The input can be defined as a list of string values corresponding to the factor number that is to be forced. eg: ['3', '4']. The default value for this parameter is an empty list.

- user_limit_for_step_two
  
  The max limit on the number of second order effects to be fit can be specified using this parameter. The input should be an integer. The default value for this parameter is "None". If the user has specified a limit, then this limit will be considered, otherwise the limit on the terms is case dependent. More information is given below:
  - No limit is set if all second order terms for all factors are jointly estimable.
  - The limit is set to the number of second order terms specified using the heredity parameters (qheredity and iheredity), if this number is less than the maximum number of jointly estimable terms for all second order effects.
  - Otherwise, the limit will be always set to run size divided by four (refer paper).

## Output

The function will auttomatically print out the following:

- Initial error degrees of freedom available
- Initial estimate of the error variance
- Main effects that are declared active
- p-values for the different main effects tested during main effects selection
- Main effects that are forced into the model (if any)
- Updated estimate of the error variance taking into account inactive main effects
- p-value for the initial F-test (Step 4a from paper)
- Limit on the number of terms for subset selection (Step 4b from paper)
- Rank of matrix with only second order terms (this gives the possible maximum number of second order terms that can be fit during subset selection)
- p-value for the final F-test (Step 4b from paper)
- Active interaction effects
- Active quadratic effects

The function outputs one return value. This value is the p-value from the last failed F-test during the second order effects selection.

## Example code

```python
output = get_omars_analysis(smat=design_matrix, cy=response, alpha=[0.05, 0.2, 0.2], qheredity='n', iheredity='n', effects_to_drop=['1_3', '2_6'], full='y', force_me=['4'], user_limit_for_step_two=None)
```
