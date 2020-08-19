# The transformed response variable is constructed to measure the spread in each group. Let

# zij = |yij - medianyj|


# where medianyj is the median of group j. The Brown–Forsythe test statistic is the model F statistic from a one way ANOVA on zij:


# F = (N-p)/(p-1) * sumjtop * nj(meanz.j - median z..) * power2 / sumjtop * sumitonj * (zij - median.j) * power2


# where p is the number of groups, nj is the number of observations in group j, and N is the total number of observations. Also meanz.j are the group means of the zij and median z.. is the overall mean of the zij. This F-statistic follows the F-distribution with degrees of freedom d1 = p - 1 and d2 = N - p under the null hypothesis.
using HypothesisTests

# X - data matrix (Size of matrix must be n-by-2; data=column 1, sample=column 2). alpha - significance level (default = 0.05).
X=[60.8 1;57.0 1;65.0 1;58.6 1;61.7 1;68.7 2;67.7 2;74.0 2;66.3 2;69.8 2; 102.6 3;102.1 3;100.2 3;96.5 3;87.9 4;84.2 4;83.1 4;85.7 4;90.3 4]

k = findmax(X[:,2])[1] # number of samples in one group

X[:,2]


diffij = groupI - groupJ
zij = abs.(diffij)
zijmedian = median(zij)
p = 2
N = 20

# numer = 
# denom = 

F = (N-p) / (p-1) * (numer) / (denom)