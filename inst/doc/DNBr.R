## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DNBr)

data(data.example)
data(meta.example)
# data(diffgenes.example)

# step1. Compute DNB's all modules of each group, return pre-result
# data.example can be sparse matrix (dgCMatrix)
# more parameter's details see ?DNBr::DNBcompute
a <- DNBcompute(
    data = data.example, 
    meta = meta.example, 
    # diffgenes = diffgenes.example
)
a
a$A

# step2. Select top N modules of each group, and use those to recompute DNB module, return result
b <- DNBfilter(
  DNB_output = a, 
  ntop = 5
)
# in the best conditions, Module candidates should be number of group_number * ntop
# however, if some group has no enough Modules meeting requirement, the large ntop may cause some
# fake result, so check the @MODULEs@bestMODULE to sure that good Module.
b
b$A

# step3. Check if the max number of exact group is equal to that you want before
# This is possible, when there are equal scores of modules for exact group, 
# or when there are not enough modules to fit your requirement.
# Use function getMaxRanking to check and get the exact value.
maxRank_A <- getMaxRanking(b, group = "A")
maxRank_A 
# get 5
# sometime this value is actually less to the ntop you assigned in last step (here 5)
# in that case, you can change the minModule/maxModule or diffgenes in step1 to get pretty result

# step4. Select the best module for all groups
# group means which group the module is computed from
# ranking means the module ranking in your group when function DNBcompute processes
df_score <- ScoreExtract(
    object = b,
    ranking = maxRank_A, # 5
    group = 'A'
)
df_score

# step5. Visualize the module's score
# input can either be the output of ScoreExtract (step4),
DNBplot(df_score, show = TRUE, save_pdf = FALSE)

# or directly object DNB_output (parameters same to ScoreExtract)
DNBplot(b, ranking = maxRank_A, group = "A", show = TRUE, save_pdf = FALSE)

# step6. Extract the whole result of DNB model for exact group
# It can be written into files or for more analysis in R.
# for slot pre_result
df_pre_allscore <- resultAllExtract(
    object = b,
    group = "A",
    slot = "pre_result",
    # mess = FALSE # avoid message
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
# 1. Compute customized Modules by DNB
# if you have an advanced list or more of module genes, to compute DNB for this module,
# here is a good resolution:
data(module_list.example)
b <- DNBcompute_custom(
    data = data.example, 
    meta = meta.example, 
    module_list = module_list.example
)
b
DNBplot(b, group = "USER_CUSTOMIZED", ranking = 1, show = TRUE, save_pdf = FALSE)
# In this case, use group = 'USER_CUSTOMIZED' when DNBplot, but DNB_output's group didn't change

# 2. Cluster by customized function to find Modules
# after using pcc to compute the sample distance matrix already and time to find Modules,
# if you do not wanna use hierachical clustering, or that of default args, 
# you can design that function as you like
# here is an example:
cluster.fun <- function(
    d, # the first arg must be d (dist matrix), necessary
    # below not necessary, according to user's needs
    method = "average", # here different to the default hclust method ('complete')
    cutree_method = 'k', # same to cutree_method, but that of DNBcompute will fail to work
    cutree_cutoff = 10 # same to cutree_cutoff, that of DNBcompute will fail to work 
){
    clu_sle <- stats::hclust(d = d, method = method) # or you can design your own cluster-method
    if (cutree_method == 'h') {
        cut_sle <- stats::cutree(clu_sle, h = cutree_cutoff) # cutoff
    } else if (cutree_method == 'k') {
        tryCatch(
            cut_sle <- stats::cutree(clu_sle, k = cutree_cutoff),
            error = function(x) stop("Not proper value of k!")
        )
    }
    
    # return the result of Module group, a vector with name->sample-name and value->group-number-name
    return(cut_sle)
}

a <- DNBcompute(
    data = data.example, 
    meta = meta.example, 
    cluster_fun = cluster.fun, # the function object, instead of function name "'cluster.fun'"
    cluster_args = NULL, # default NULL, or you can add extra args (assign different value of args except d in each run) to better adapt
    # cluster_args = list(cutree_method = 'k', cutree_cutoff = 10),
    quiet = TRUE # not verbose this
)
b <- DNBfilter(
  DNB_output = a, 
  ntop = 5
)
DNBplot(b, group = "A", show = TRUE, save_pdf = FALSE)



