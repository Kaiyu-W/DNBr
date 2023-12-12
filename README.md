# DNBr
This code is written by `R`, to compute the ___Dynamic Network Biomarkers___(__`DNB`__) model developed by ChenLab in CAS, China. 
  
If you're interested in this subject, please  
- go to http://sysbio.sibcb.ac.cn/sysbio/cb/chenlab/LuonanChen.htm or  
- email lnchen@sibs.ac.cn to let us know.  
  
To install this package, use the following command in R:  
```R
devtools::install_github("Kaiyu-W/DNBr")
```
  
More details see the codes by `vignette("DNBr")` or `?DNBr::function` in R, please.  
  
__Log__:  
- Add chioce `high_cutoff = -1` for DNBcompute() to skip selecting highly variable genes.  
- Add parameter `force_allgene` for DNBfilter() to control which genes are used as background.  
  
- Note that,  
  - To stay consistent with DNBcompute(), DNBfilter(force_allgene = FALSE) as default.  
  - There may be some missing genes when DNBfilter(force_allgene = FALSE) because of different high variable genes among groups, and those missing ones will later be removed from modules during computation. Use `DNBcompute(..., high_cutoff = -1)` to avoid this.
  - Use `DNBfilter(..., force_allgene = TRUE)` to re-produce the output by early version.  
  
__Example__:
```R
library(DNBr)
data(data.example)
data(meta.example)

a <- DNBcompute(data.example, meta.example)
b <- DNBfilter(a, ntop = 5, force_allgene = TRUE)
DNBplot(b, show = TRUE, save_pdf = FALSE)
```
