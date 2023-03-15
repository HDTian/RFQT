# RFQT
This repository merging to the main branch illustrates the key messages to draw the **Random Forest of Q trees (RFQT)** and its relevant analysis for treatment heterogeneity studies with Mendelian randomization, it can achieve a lots of work, including
- Draw heterogeneity Q test 
- Build single Q tree  
- Build random forest of Q trees 
- Make treatment effect prediction
- Evaluate Variable Importance (VI) 
- draw permutation test


When drawing your treatment heterogeneity study, define the one-dimensional instrument *Z*, the exposure of interest *X*, the (high-dimensional) covariates *M*, and the outcome of interest *Y*. If the instrument is high-dimensional, first build the single instrument (e.g. by weighted gene score). The standard data structure regonized by RFQT is *>cbind(individual index, Z , X , Y, M)*.
## Toy example
Below is one example of data used for simulation considering 20 candidiate covariates where the first five covariates are the treatment effect modifiors  
$$X = 0.5Z + 0.5 \sum_{j=1}^{20} U_j  + \epsilon_X $$ 
$$M_j  =  b_j X + U_j \qquad j=1,2,\ldots,20 $$
$$Y = \left(   0.5+ \sum_{j=1}^{5} \gamma_j M_j  \right) X +  0.5 \sum_{j=1}^{20} U_j  + \epsilon_Y  $$
where $b_j$ is the effect of the exposure on the covariate $M_j$; when $b_j$ is not zero, it causes the collider bias conditioning on $M_j$. $\gamma_j$ is the strength of modification for $M_j$. One need to evaluate the heterogenous treatment effects based on different levels of effect modifiors. RFQT is a collider robust methods combing with data-adaptive methods that can let you estimate reliable individual treatment effects, find out the potential effect modifors from high-dimensional covariates, and many other relevant analysis.
## Build a single Q tree
Constructing a single Q tree caintains two iterative steps: 
1. Determine the splitting covariate
2. Split based on the candidate covariate
![Q-tree.pdf](https://github.com/HDTian/RFQT/files/10981826/Q-tree.pdf)
![Q-tree.pptx](https://github.com/HDTian/RFQT/files/10981826/Q-tree.pptx)
