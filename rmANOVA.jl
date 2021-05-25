## code being ported from: ...\Dynamic_CSD\extra toolboxes\_rmANOVA\teg_repeated_measures_ANOVA.m 

# (example of use in ...\Dynamic_CSD\Deane2020\avrec_relres_stats.m)

home    = @__DIR__

# rawX is observation x variable-combination matrix.
# levels is vector of levels per factor.
# variables is a cell array of strings.

## test data:
levels    = [2,2] # 2 groups, 2 measurements
variables = ["Groups" "Measurement"]
rawX      = randn(10,4) # rows: 2 groups of 10 subjects each | columns: 2 groups of 2 measurements each - horizontally concatonated
#function rmAnova(rawX, levels=[2,2], variables=["Var1" "Var2"])

# perform a quick check
if size(datX,2) != prod(levels,dims=1)[1]
    error("Ensure the the number of columns equals the product of the levels. E.g.: If there are 2 groups and 2 measurements, there should be 4 columns. 2 goups and 4 measurements should have 8 columns. Horizontally concatonate the groups' vertically stacked measurement results.")
end

nSub, nVar = size(rawX,1), size(rawX,2)

# remove and store subject effects
subEffect = ones(nSub)
for iSub = 1:nSub
    fNonNaN = isfinite.(rawX[iSub,:]) # NaN will trigger "false"
    subEffect[iSub] = sum(rawX[iSub,fNonNaN] .* onesX[iSub,fNonNaN]) ./ sum(onesX[iSub, fNonNaN])
end
rawX_sub = subEffect * ones(1,size(rawX,2))
datX = rawX - rawX_sub

# preallocate ANOVA fields
# X1, factorStarts, nColsFactor, labels, cellsets = preANOVA(levels, variabls)


p_boncor = 0.05 / prod(levels,dims=1)  # bonferroni corrected p value


