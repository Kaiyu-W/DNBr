# DNBr
This code is written by `R`, to compute the ___Dynamic Network Biomarkers___(__`DNB`__) model developed by ChenLab in CAS, China. 

If you're interested in this subject, please  
- go to http://sysbio.sibcb.ac.cn/sysbio/cb/chenlab/LuonanChen.htm or  
- email lnchen@sibs.ac.cn to let us know.  
  
To install this package, use the following command in R:  
```R
devtools::install_github("Kaiyu-W/DNBr")
```

Note that, add force_allgene = TRUE when DNBfilter() to re-produce the output by early version.   
This parameter is designed to control which genes are used as background (dnb_in/dnb_out when computing sd_in,pcc_in/pcc_out). Stay consistent with DNBcompute() for DNBfilter() as default. However, if using highly variable genes when DNBcompute(), there may be some missing genes when DNBfilter() because of different high variable genes among groups, and those missing ones will be removed from modules. Use force_allgene = TRUE to use all same genes with no module-in's gene difference but maybe different module-out's gene.   
Use high_cutoff = -1 when DNBcompute() to avoid selecting highly variable genes. How to set options depends on user's condition.    
More details see the codes by `vignette("DNBr")` or `?DNBr::function` in R, please.  
  
__Example__:
```R
library(DNBr)
data(data.example)
data(meta.example)

a <- DNBcompute(data.example, meta.example)
b <- DNBfilter(a, ntop = 5)
DNBplot(b, show = TRUE, save_pdf = FALSE)
```
