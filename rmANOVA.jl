## code being ported from: ...\Dynamic_CSD\extra toolboxes\_rmANOVA\teg_repeated_measures_ANOVA.m 

# (example of use in ...\Dynamic_CSD\Deane2020\avrec_relres_stats.m)

home    = @__DIR__

# datX is observation x variable-combination matrix.
# levels is vector of levels per factor.
# variables is a cell array of strings.

## test data:
levels    = [2,2] # 2 groups, 2 measurements
variables = ["Groups" "Measurement"]
datX      = randn(10,4) # rows: 2 groups of 10 subjects each columns: 2 groups of 2 measurements each - horizontally concatonated

p_boncor = 0.05 / prod(levels,dims=1)  # bonferroni corrected p value

if size(datX,2) > prod(levels,dims=1)[1] # n of observations per cell
    b    = size(datX,2) / 2
    newX = datX[:,Int(b)+1:end]
    datX = datX[:,1:Int(b)] 
else
    newX = ones(size(datX))
end

nSub, nVar = size(datX,1), size(datX,2)
