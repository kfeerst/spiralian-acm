# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Sep 2017

# based on http://blog.phytools.org/2013/03/estimating-ancestral-states-when-tips.html
# and http://blog.phytools.org/2017/09/densitymap-like-plot-for-3-state.html
# script runs two stochastic character mappings, one with the ctenophore
#     sister tree and one with the sponge sister tree.

# set seed so script gives same answer each run
set.seed(420)

# if these packages are not installed 
# install.packages("maps")
# install.packages("ape")
# install.packages("phytools")

library(maps)
library(ape)
library(phytools)

# set variables; these can be adjusted if running on different data
mydir   <- "/Users/katie/Desktop/KATIE"
mynsim  <- 1000
mytree  <- "working_phy_pruned.tre"

model_cs <- "ER"
model_ed <- "ER"
model_lt_ard <- "ARD"
model_lt_er  <- "ER"

outpdf_cs <- "simmap_out_cs.pdf"
outpdf_ed   <- "simmap_out_ed.pdf"
outpdf_lt_ard   <- "simmap_out_lt_ard.pdf"
outpdf_lt_er    <- "simmap_out_lt_er.pdf"

dens_outpdf_cs <- "simmap_densityTree_cs.pdf"
dens_outpdf_ed <- "simmap_densityTree_ed.pdf"
dens_outpdf_lt_ard <- "simmap_densityTree_lt_ard.pdf"
dens_outpdf_lt_er  <- "simmap_densityTree_lt_er.pdf"

charmat_vals_cs <- "charmat_vals_cs.csv"
charmat_vals_ed  <- "egg_size_simmap_wdec.numbers.csv"
charmat_vals_lt  <- "larval_simmap_wdec.numbers.csv"

charmat_names_cs <- "charmat_names_cs.csv"
charmat_names_ed <- "egg_size_simmap_wdec.names.csv"
charmat_names_lt <- "larval_simmap_wdec.names.csv"

traitvals_cs <- c("equal","unequal")
traitvals_ed <- c("large","small")
traitvals_lt <- c("trochophore","non-trochophore")

setwd(mydir)

# function runs make.simmap and densityTree
mysimmap <- function(mytree,mypdfout,mypdfdensityout,cmv,cmn,traitvals,mod) {

    # set output pdf file
    pdf(file=mypdfout,width=8.5,height=11)

    # read in tree and matrix; matrix is broken up into data and rownames
    # outgroups (Monosiga_brevicollis and Proterospongia) have been removed
    anitree <- read.newick(mytree)
    anidata <- read.table(cmv,sep=",")
    rn <- read.table(cmn)
    animatrix <- as.matrix(anidata)
    rownames(animatrix) <- rn$V1
    colnames(animatrix) <- c(traitvals)

    # run make.simmap and print out information about the run
    SYM.simmap_trees <- make.simmap(anitree,animatrix,nsim=mynsim,model=mod)
    SYM.simmap_trees$loglike
    res_simmap <- describe.simmap(SYM.simmap_trees)
    print(res_simmap)
    print (res_simmap$ace)

    # make psuedo chrono tree for plotting
    anichrono <- chronopl(anitree, lambda = 0, age.min = 1)

    # plot tree with labels slightly offset so pies can go at tips
    plot(anichrono,label.offset=.01, cex=0.4)
    nodelabels(pie=res_simmap$ace,piecol=c("red","blue"),cex=0.2)
    tiplabels(pie=res_simmap$tips,piecol=c("red","blue"),cex=0.2)

    # reset pdf file for densityTree (getting blank 1st page so onfile=FALSE)
    dev.off()
    pdf(file=mypdfdensityout,onefile=FALSE,width=8.5,height=11)

    # set colors; adjust margins and fontsize; 
    colors <- setNames(c(c="red","blue"),traitvals)
    par(mar=c(1,0,0,1),cex=0.6)

    # run densityTree
    densityTree(SYM.simmap_trees,method="plotSimmap",lwd=3,colors=colors);

    cat("\n------------------------------------------------------\n")
}

# call our simmap function twice with different traitsets
mysimmap(mytree,outpdf_ed,dens_outpdf_ed,charmat_vals_ed,charmat_names_ed,traitvals_ed,model_ed)
mysimmap(mytree,outpdf_lt_ard,dens_outpdf_lt_ard,charmat_vals_lt,charmat_names_lt,traitvals_lt,model_lt_ard)
mysimmap(mytree,outpdf_lt_er,dens_outpdf_lt_er,charmat_vals_lt,charmat_names_lt,traitvals_lt,model_lt_er)
mysimmap(mytree,outpdf_cs,dens_outpdf_cs,charmat_vals_cs,charmat_names_cs,traitvals_cs,model_cs)

# print versions of all loaded modules so analysis can be reproduced
sessionInfo()
