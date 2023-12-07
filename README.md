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
