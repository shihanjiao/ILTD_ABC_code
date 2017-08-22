clear
clc
tic

% This script is created for testing ILTD_ABC's performance 

NP=40;
MAXITER=1000;
dimension=10;
runno=20;
objfun='f01';
irange_l=-100;
irange_r=100;


for i=1
   
    [meanfitness(i,1), minium(i,1), maxmum(i,1), meanposition(i,1,:), SD(i,1),  dat2(i,1,:)] =iltdabc (NP,MAXITER,dimension,runno,irange_r(i),irange_l(i),objfun(i,:),i);
   
    meanfitness(i,:)
   
end


toc