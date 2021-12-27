## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DNBr)

data(data.example)
data(meta.example)
data(diffgenes.example)

# step1. Compute DNB's all modules of each group, return pre-result
# more parameter's details see ?DNBr::DNBcompute
a <- DNBcompute(
    data = data.example, 
    meta = meta.example, 
    diffgenes = diffgenes.example
)
a
a$A

# step2. Select top N modules of each group, and use those to recompute DNB module, return result
b <- DNBfilter(
  DNB_output = a, 
  ntop = 5
)
b
b$A

# step3. Check if the max number of exact group is equal to that you want before
# This is possible, when there are equal scores of modules for exact group, 
# or when there are not enough modules to fit your requirement.
maxRank_A <- getMaxRanking(b, group = "A")
maxRank_A 
# get 3 instead of 5

# step4. Select the best module for all groups
# group means which group the module is computed from
# ranking means the module ranking in your group when function DNBcompute processes
df_score <- ScoreExtract(
    object = b,
    ranking = maxRank_A,
    group = 'A'
)
df_score

# step5. Visualize the module's score
# input can either be the output of ScoreExtract (step4),
DNBplot(df_score, show = TRUE, save_pdf = FALSE)

# or directly object DNB_output (parameters same to ScoreExtract)
DNBplot(b, ranking = maxRank_A, group = "A", show = TRUE, save_pdf = FALSE)

# step6. Extract the whole result of DNB model for exact group
# for slot pre_result
df_pre_allscore <- resultAllExtract(
    object = b,
    group = "A",
    slot = "pre_result"
)
df_pre_allscore

# for slot result
df_allscore <- resultAllExtract(
    object = b,
    group = "A",
    slot = "result"
)
df_allscore

# slot rank refers to the ranking for each module inside its own group 
# slot rank_all refers to the ranking for each module in all modules of that group 

# *optional
# if you have an advanced list or more of module genes, to compute DNB for this module,
# here is a good resolution:
data(module_list.example)
b <- DNBcompute_custom(
    data = data.example, 
    meta = meta.example, 
    module_list = module_list.example
)
b
DNBplot(b, group = "D", ranking = 1, show = TRUE, save_pdf = FALSE)
# In this case, all groups (resources) will get the customized module lists, so group can be anyone (ignored) !


