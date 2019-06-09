# DCDTr
Differential Network Inference with Compositional Data via Lasso Penalized D-Trace Loss

a)	The R code in folder "DCDTr\CDTrace" is an implementation of CDTr in our manuscript, which is available at https://github.com/coamo2/CDTr.

b)  The folder "DCDTr\dpm-master" an implementation of CDtrace, which is available at http://cran.r-project.org/web/packages/flare/index.html. Compile "dpm.c" using "R CMD SHLIB dpm.c" first!

c)	The folder "DCDTr\SpiecEasi-master" is an implementation of FGL and GGL (Zhao et al, 2014), which is available at https://github.com/zdk123/SpiecEasi by Kurtz et al (2015).

d) 	The folder  "DCDTr\gcoda-master” is an implementation of gCoda (Fang et al, 2017), which is available from https://github.com/huayingfang/gCoda by Fang et al (2017).

e)  The R package "Difdtl_1.0" is an implementation of DTL (Yuan et al, 2015), which is also an implementation of DCDTr in our manuscript.

e)	The R files named "DCDsim.R" and "simfunc" are for simulations under 6 different scenarios. The differential networks are estimated according to DCDTr, FGL, GGL and L1-M with the help of above-mentioned codes.
