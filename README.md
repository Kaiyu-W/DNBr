# DNBr
This code is written by `R`, to compute the ___Dynamic Network Biomarkers___(__`DNB`__) model developed by ChenLab in CAS, China. 

If you're interested in this subject, please  
- go to http://sysbio.sibcb.ac.cn/sysbio/cb/chenlab/LuonanChen.htm or  
- email lnchen@sjtu.edu.cn to let us know.  
  
To install this package, use the following command in R:  
```R
devtools::install_github("Kaiyu-W/DNBr")
```

Notes:
- Reproducibility: Set `force_allgene = TRUE` in `DNBfilter()` to reproduce outputs from earlier versions.
- Performance & Memory Strategy:
  - Recommended (Single-Process Multi-threading): Prioritize a multi-threaded BLAS backend (OpenBLAS or MKL) with `precomputed_cor = TRUE` and `fast = FALSE`. Since the bottleneck is the correlation matrix computation, a multi-threaded BLAS provides significant acceleration within a single process without extra memory overhead. 
  - Alternative (Multi-Process Parallelism): Using `future` with `fast = TRUE` can further accelerate correlation calculations via parallel forks (multicore on Unix-like systems). While this can further accelerate group-level computations, it drastically increases memory consumption due to data replication across workers. It may also lead to serialization errors (e.g., "long vectors not supported") when using `multisession` (common on Windows or when `avoid_multicore = TRUE`). _Note: Multi-process mode has not been extensively tested on Windows_.
  - Memory Warning: Please ensure your system has sufficient RAM, especially when analyzing large datasets with high gene counts, as pre-computing correlation matrices requires much memory.

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
