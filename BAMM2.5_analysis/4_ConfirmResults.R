## libraries
library(RPANDA)
library(TreePar)
library(BAMMtools)
library(TESS)


### Data ---
### TreePar with Solanum ---
Solanum_incom_tree<-read.tree("./Data/Solanum_Sarkinen_non_negative.tre")
Solanum_sort_tree<-sort(getx(Solanum_incom_tree),decreasing=TRUE)## is the same as rev(sort(branching.times(Cetacea)))

# When estimating the the rate shift times t based on branching times x, 
# we allow the shift times to be 0.335 every 0.1 until max(branching times x) :
start <- 0.12
end <- max(Solanum_sort_tree)
grid <- 1
rho<-c(0.34,1)## How many rates are allowed. 0.5 == 50% of the species to simulate incomplete sampling

## To evaluate for three/four shifts of diversification the sampling =c(rho,1,1) where rho c(0.34,1) says that the 34% of the species at present are sampled
res_solanum_incom<-bd.shifts.optim(Solanum_sort_tree,c(rho,1,1),start=start, end=end ,grid=grid)
res_solanum_incom_tmp<-res_solanum_incom[[2]]

i<-1

test<-pchisq(2*(res_solanum_incom_tmp[[i]][1]-res_solanum_incom_tmp[[i+1]][1]),3)
test
## Since 0 shifts are not significant better than 0 then 0 continues to be the reference point 

#test if 2 shifts explain the tree significantly better than 0 shift:
i<-2
test<-pchisq(2*(res_solanum_incom_tmp[[i]][1]-res_solanum_incom_tmp[[i+1]][1]),3)
test

#test if 2 shifts explain the tree significantly better than 0 shift:
i<-3
test<-pchisq(2*(res_solanum_incom_tmp[[i]][1]-res_solanum_incom_tmp[[i+1]][1]),3)
test

## In this case 2 shifts are the ones that best explain the data

pdf("./Figures/Solanum_treePar.pdf")
bd.shifts.plot(list(res_solanum_incom[[2]]),2,ratemax=2,timemax = max(Solanum_sort_tree),plotturnover=TRUE)
legend("topright",c("lambda","mu"),lty=c(1,1),col=c("blue","red"))
dev.off()

save(res_solanum_incom,file="./Output/TreePar_Solanum_incomplete_sampl.RData")

### Fully resolved Solanum 
Solanum_fully_run1<-read.tree("./Data/Input_BAMM2.5/1/tree.tre")

Solanum_x<-sort(getx(Solanum_fully_run1),decreasing=TRUE)## is the same as rev(sort(branching.times(Cetacea)))

# When estimating the the rate shift times t based on branching times x, 
# we allow the shift times to be 0.335 every 0.1 until 35 :

start <- 0.12
end <- max(Solanum_x)
grid <- 1
rho<-c(1,1)

res_solanum_fully<-bd.shifts.optim(Solanum_x,c(rho,1,1),start=start, end=end,grid=grid)
res_solanum_fully_tmp<-res_solanum_fully[[2]]

bd.shifts.plot(list(res_solanum_fully[[2]]),2,ratemax=2,timemax = max(Solanum_x),plotturnover=TRUE)

i<-1
test<-pchisq(2*(res_solanum_fully_tmp[[i]][1]-res_solanum_fully_tmp[[i+1]][1]),3)
test
## Since 0 shifts are not significant better than 0 then 0 continues to be the reference point 

#test if 2 shifts explain the tree significantly better than 0 shift:
i<-2
test<-pchisq(2*(res_solanum_fully_tmp[[i]][1]-res_solanum_fully_tmp[[i+1]][1]),3)
test

save(res_solanum_fully,file="./Output/TreePar_Solanum_fully_run1.RData")

