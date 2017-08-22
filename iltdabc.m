function [meanfitness,minium,maxmum,meanposition,SD,data2]=iltdabc(NP,MAXITER,dimension,runno,irange_r,irange_l,objfun,tt,limit)

% intelligent learning with turbulent distribution in artificial bee
% colony and its application to multi-level image segmentation


FoodNumber=NP/2;
limit=100; 
evaluations=NP*MAXITER;
data1=zeros(runno,evaluations);

GlobalMins=zeros(1,runno);

for run=1:runno
  
rand('state',sum(run+tt*10))

for i=1:FoodNumber
    Foods(i,:) =irange_l+(irange_r-irange_l)*rand(1,dimension);
end

for i=1:FoodNumber
ObjVal(i)=feval(objfun,Foods(i,:));
Fitness=calculateFitness(ObjVal);
end


trial=zeros(1,FoodNumber);


BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);

iter=1;
fes=0;

memory=zeros(MAXITER,dimension);
memory_fit=zeros(1,MAXITER);

sign_step=(-1).^round(rand(FoodNumber,dimension));


while ((iter <= MAXITER)),
        
    [fit_sort,index]=sort(ObjVal);
    
    for i=1:FoodNumber
        if i~=index(FoodNumber) || iter<2
        sol=Foods(i,:);
        j=fix(rand*dimension)+1;
        mbest=sum(Foods)/FoodNumber;
        orientation=sign(mbest(j)-Foods(i,j))*sign(GlobalParams(j)-Foods(i,j));
        
        r1=rand;
        r2=rand;

        beta=r1+sin(2*pi*r1)/(2*pi);
        alpha=r2+sin(2*pi*r2)/(2*pi);

        if orientation==-1
            step=sign_step(i,j)*abs((GlobalParams(j)-Foods(i,j))*beta-(mbest(j)-Foods(i,j))*alpha);
        else
            step=(mbest(j)-Foods(i,j))*alpha+(GlobalParams(j)-Foods(i,j))*beta; 
        end
        
        sol(j)=Foods(i,j)+step;
       
        
        ind=find(sol<irange_l);
        sol(ind)=irange_l;
        ind=find(sol>irange_r);
        sol(ind)=irange_r;
       
        ObjValSol=feval(objfun,sol);
        FitnessSol=calculateFitness(ObjValSol);
   
       if (ObjValSol<ObjVal(i)) 
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
            sign_step(i,j)=sign(step);
       else
           sign_step(i,j)=-sign(step);
           trial(i)=trial(i)+1; 
       end
       
        else
            j=fix(rand*dimension)+1;
            
            temp(1:2,:)=memory(iter-1:iter,:);
            trial_fit(1:2)=memory_fit(iter-1:iter);
            ii=randperm(2);
            mean_pos=sum(temp)/2;
            for d=1:dimension
                std_pos(d)=std(temp(:,d));
                if std_pos(d)<0.1
                    std_pos(d)=mean_pos(d);
                end
            end
            x=rand(1,dimension);
            zz=(-1).^round(rand(1,dimension)).*(rand(1,dimension)-(2/(pi))*sin(pi*x))*2;

           
            
            sol=Foods(i,:);
            sol=mean_pos+zz.*std_pos;
            ind=find(sol<irange_l);
            sol(ind)=irange_l;
            ind=find(sol>irange_r);
            sol(ind)=irange_r;
       
        ObjValSol=feval(objfun,sol);
        FitnessSol=calculateFitness(ObjValSol);
   
       if (ObjValSol<ObjVal(i)) 
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
            
       else
           
           trial(i)=trial(i)+1; 
       end
        end
            
       
       if ObjVal(i)<GlobalMin
           GlobalMin=ObjVal(i);
           GlobalParams=Foods(i,:);
       end
       
       fes=fes+1;
       data1(run,fes)=GlobalMin;
         
    end

    [fit_sort,index]=sort(ObjVal);
    
  pro=Fitness./sum(Fitness);
  prob=zeros(1,NP/2);
  for ii=1:NP/2
      for iii=1:ii
          prob(ii)=prob(ii)+pro(iii);
      end
  end
  
      H=(FoodNumber*ObjVal)./(sum(ObjVal));
      pool=[];
      for i=1:FoodNumber
          if H(i)<=1
              pool=[i,pool];
          end
      end
  
  if isempty(pool)
      H=rand(1,FoodNumber);
      pool=index;
  end
  
 for i=1:FoodNumber

    j=fix(rand*dimension)+1;
    sol=Foods(i,:);
    
    r1=fix(rand*length(pool))+1;
    sol(j)=Foods(pool(r1),j)+(H(pool(r1)))*(Foods(pool(r1),j)-Foods(i,j));

    ind=find(sol<irange_l);
    
        sol(ind)=irange_l;
        ind=find(sol>irange_r);
        sol(ind)=irange_r;
       
        ObjValSol=feval(objfun,sol);
        FitnessSol=calculateFitness(ObjValSol);
   
       if (ObjValSol<ObjVal(i)) 
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
       else
           trial(i)=trial(i)+1; 
       end
       if ObjVal(i)<GlobalMin
           GlobalMin=ObjVal(i);
           GlobalParams=Foods(i,:);
       end
       
       fes=fes+1;
       data1(run,fes)=GlobalMin;
        
end; 

        
        

memory(iter,:)=GlobalParams;
memory_fit(iter)=GlobalMin;

                  
%%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       
ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    Bas(ind)=0;
    sol=irange_l+(irange_r-irange_l)*rand(1,dimension);
    ObjValSol=feval(objfun,sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
    if ObjVal(ind)<GlobalMin
           GlobalMin=ObjVal(ind);
           GlobalParams=Foods(ind,:);
    end  
    
    fes=fes+1;
    data1(run,fes)=GlobalMin;
end;

if fes>=evaluations
    break;
end

iter=iter+1;

end % End of ABC

GlobalMins(run)=GlobalMin;
GlobalParam(run,:)=GlobalParams;
end; 

sum1=0;
data=zeros(1,dimension);
for i=1:runno
    sum1=sum1+GlobalMins(i);
    data=data+GlobalParam(i,:);
end
meanfitness=sum1/runno;%输出平均值
minium=min(GlobalMins(:));
maxmum=max(GlobalMins(:));
SD=std(GlobalMins(:));
meanposition=data/runno;

data3=0;
  data2=zeros(1,evaluations);
  for i =1:evaluations
      for j=1:runno
          data3=data3+data1(j,i);
      end
      data2(i)=data3/runno;
      data3=0;
  end

data5=data1(:,1:evaluations); 

end

function fFitness=calculateFitness(fObjV)
fFitness=zeros(size(fObjV));
ind=find(fObjV>=0);
fFitness(ind)=1./(fObjV(ind)+1);
ind=find(fObjV<0);
fFitness(ind)=1+abs(fObjV(ind));
end

