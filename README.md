# Documentation for package omars-analysis

This package has one function 'get_omars_analysis'

First run the following:

```python
from omars_analysis.main import get_omars_analysis
```

Running the above command will make the make the function 'get_omars_analysis' available for use.

## Function usage

The function requires six inputs:

- mat
  
  the first input is the design matrix. The design matrix need not be coded, however it must only consist of continuous variables. The function is meant to be used with designs that have all factors either with two levels or three levels. The design should **NOT** consist of headers. The design matrix should not consist of second order effects since this will be built internally in the function.
- cy
  
  This is the response.

- qheredity
  
  This is to specify heredity constraints for quadratic effects. The accepted inputs are 'y' or 'n' ('y'- strong heredity, 'n'- no heredity, 'n'- No heredity). The input must be a string in lowercase.
- iheredity
  
  This is to specify heredity constraints for two-factor interaction effects. The accepted inputs are 's', 'w' or 'n' ('s'- strong heredity, 'w'- weak heredity, 'n'- no heredity). The input must be a string in lowercase.

- effects_to_drop
  
  This is to specify second order effects that must be excluded from the analysis. The input must be a list of strings. For example: ['1_1', '2_3']. This input specifies that the quadratic effect of the first factor and the interaction effect between factor two and three must be excluded from the second step of the analysis (subset selection).

- full
  
  'n' -  analysis is performed on the main effects only
  
  'y' - analysis is performed on the main effects and second order effects.

  The default is set to 'y'

## Output

The function will auttomatically print out the following:

- Initial error degrees of freedom available
- The initial estimate of the error variance
- Main effects that are active
- Updated estimate of the error variance
- Active interaction effects
- Active quadratic effects

The function outputs one return value. This value is the p-value from the last failed F-test during the second order effects selection.

## Example code

```python
output = get_omars_analysis(mat=design_matrix, cy=response, qheredity='n', iheredity='n', effects_to_drop=[], full='y')
```
