# DNBr
The codes is written by `R`, to compute the ___Dynamic Network Biomarkers___(__`DNB`__) model developed by ChenLab in CAS, China. 

If you're interested in this subject, please  
- go to http://sysbio.sibcb.ac.cn/sysbio/cb/chenlab/LuonanChen.htm or  
- email lnchen@sibs.ac.cn to let us know.  
  
To install this package, use the following command in R:  
```R
devtools::install_github("Kaiyu-W/DNBr")
```

More details see the codes by `vignette("DNBr")` or `?DNBr::function` in R, please.  
  
__Example__:
```R
library(DNBr)
data(data.example)
data(meta.example)
data(diffgenes.example)

a <- DNBcompute(data.example, meta.example, diffgenes.example)
b <- DNBfilter(a, ntop = 5)
DNBplot(b, show = TRUE, save_pdf = FALSE)
```
