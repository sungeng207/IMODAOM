%___________________________________________________________________%
%  Multi-Objective Dragonfly Algorithm (MODA) source codes demo     %
%                           version 1.0                             %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%   S. Mirjalili, Dragonfly algorithm: a new meta-heuristic         %
%   optimization technique for solving single-objective, discrete,  %
%   and multi-objective problems, Neural Computing and Applications %
%   DOI: http://dx.doi.org/10.1007/s00521-015-1920-1                %
%___________________________________________________________________%

clc;
clear;
close all;

bs_num=8; 
uav_num=16; 
cluster_num=4;

uav_origin=[0.427062347759695,0.955372102971021,0.724247033520084,0.580891712363774,0.540257907420780,0.705441191564081,0.00502888330112106,0.782515778836242,0.926859573385509,0.00829565739263971,0.824628342666394,0.767335868027880,0.997136895348686,0.227653072058219,0.919542206194465,0.641999305521780,7.38418767360661,12.0508166401978,98.1596208916383,49.6799421779149,2.24136544401531,5.38315489695107,14.0873797601841,89.3474311740186,46.5820073116048,56.0856738133542,49.4456347383399,6.77854884858151,89.7646542373691,28.8565308881261,26.9046814753517,59.4194175430363,47.5879024056308,36.8311038073403,65.5611090496500,93.8200405804263,62.0425185783994,28.2840091614325,20.5181251932987,43.9134081593786,2.72502196608166,87.6184335960346,61.0092246274455,20.3592380287419,51.9916824790401,5.38243021501204,86.2187431116878,44.2934666053869,89.9231313097571,83.3401165511889,76.3135368534335,60.5978347118252,84.2573618067275,73.3270757440223,77.5940116085254,67.1302176736209,80.2231478168383,78.9735710082504,69.7910978255359,74.9390363957855,81.0447845777941,64.4490864632728,89.9044221513849,79.3794655977227];

BSs=[2500,2500,0;  2500,1500,0;  5500,2550,0;  2500,130,0;  3000,3000,0;  3500,0,0;   3500,1500,0;  1750,1750,0
     -2500,2500,0;  -2500,1500,0;  -5500,2550,0;  -2500,30,0;  -3000,3000,0;  -3500,0,0;   -3500,1500,0;  -1750,1750,0
    2500,-2500,0;  2500,-1500,0;  5500,-2550,0;  2500,-30,0;  3000,-3000,0;  3500,0,0;   3500,-1500,0;  1750,-1750,0
     -2500,-2500,0;  -2500,-1500,0;  -5500,-2550,0;  -2500,-130,0;  -3000,-3000,0;  -3500,-400,0;   -3500,-1500,0;  -1750,-1750,0
     ];

Evsdp=[3400,1850,0;  -2100,2631,0;  2666,-223,0;  -3000,-1567,0];


obj_no=3;
max_iter=500;
N=30;
ArchiveMaxSize=30;
dim = 4*uav_num*cluster_num+cluster_num+cluster_num; 
lb=cell2mat(struct2cell(load('_data_lb')));
ub=cell2mat(struct2cell(load('_data_ub')));
ObjectiveFunction=@EO_evaluate_objective; 


%% main program

% Initial parameters of the MODA algorithm

Archive_X=zeros(30,dim);
Archive_F=ones(30,obj_no)*inf;

Archive_member_no=0;

r=(ub-lb)/2;
V_max=(ub(1)-lb(1))/10;

Food_fitness=inf*ones(1,obj_no);
Food_pos=zeros(dim,1);

Enemy_fitness=-inf*ones(1,obj_no);
Enemy_pos=zeros(dim,1);
X=Improved_initialize_variables(N,dim,ub,lb, bs_num,uav_num,cluster_num, BSs, Evsdp, uav_origin);
X=X';
fitness=zeros(N,2);

DeltaX=MODA_initialization(N,dim,ub,lb);

position_history=zeros(N,max_iter,dim);

for iter=1:max_iter
    
    for i=1:N
        X(:,i)=Improved_detect_coll(X(:,i)',dim, uav_num, cluster_num)';
    end
    
    r=(ub-lb)/4+((ub-lb)*(iter/max_iter)*2);
    
    w=0.9-iter*((0.9-0.2)/max_iter);
    
    my_c=0.1-iter*((0.1-0)/(max_iter/2));
    if my_c<0
        my_c=0;
    end
    
    if iter<(3*max_iter/4)
        s=my_c;             % Seperation weight
        a=my_c;             % Alignment weight
        c=my_c;             % Cohesion weight
        f=2*rand;           % Food attraction weight
        e=my_c;             % Enemy distraction weight
    else
        s=my_c/iter;        % Seperation weight
        a=my_c/iter;        % Alignment weight
        c=my_c/iter;        % Cohesion weight
        f=2*rand;           % Food attraction weight
        e=my_c/iter;        % Enemy distraction weight
    end
    
    for i=1:N %Calculate all the objective values first
        
        Particles_F(i,:)=ObjectiveFunction(X(:,i)', dim, uav_num, bs_num, cluster_num, BSs, Evsdp,uav_origin);
        
        if MODA_dominates(Particles_F(i,:),Food_fitness)
            Food_fitness=Particles_F(i,:);
            Food_pos=X(:,i);
        end
        
        if MODA_dominates(Enemy_fitness,Particles_F(i,:))
            if all(X(:,i)<ub') && all( X(:,i)>lb')
                Enemy_fitness=Particles_F(i,:);
                Enemy_pos=X(:,i);
            end
        end
    end
    
    [Archive_X, Archive_F, Archive_member_no]=MODA_UpdateArchive(Archive_X, Archive_F, X, Particles_F, Archive_member_no);
    
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=MODA_RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
        [Archive_X, Archive_F, Archive_mem_ranks, Archive_member_no]=MODA_HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=MODA_RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
    end
    
    Archive_mem_ranks=MODA_RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
    
    % Chose the archive member in the least population area as foods
    % to improve coverage
    index=MODA_RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
    Food_fitness=Archive_F(index,:);
    Food_pos=Archive_X(index,:)';
       
    % Chose the archive member in the most population area as enemies
    % to improve coverage
    index=MODA_RouletteWheelSelection(Archive_mem_ranks);
    if index==-1
        index=1;
    end
    Enemy_fitness=Archive_F(index,:);
    Enemy_pos=Archive_X(index,:)';

    for i=1:N
        index=0;
        neighbours_no=0;
        
        clear Neighbours_V
        clear Neighbours_X
        % Find the neighbouring solutions
        for j=1:N
            Dist=MODA_distance(X(:,i),X(:,j));
            if (all(Dist<=r) && all(Dist~=0))
                index=index+1;
                neighbours_no=neighbours_no+1;
                Neighbours_V(:,index)=DeltaX(:,j);
                Neighbours_X(:,index)=X(:,j);
            end
        end
        
        % Seperation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.1)
        S=zeros(dim,1);
        if neighbours_no>1
            for k=1:neighbours_no
                S=S+(Neighbours_X(:,k)-X(:,i));
            end
            S=-S;
        else
            S=zeros(dim,1);
        end
        
        % Alignment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.2)
        if neighbours_no>1
            A=(sum(Neighbours_V')')/neighbours_no;
        else
            A=DeltaX(:,i);
        end
        
        % Cohesion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.3)
        if neighbours_no>1
            C_temp=(sum(Neighbours_X')')/neighbours_no;
        else
            C_temp=X(:,i);
        end
        
        C=C_temp-X(:,i);
        
        % Attraction to food%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.4)
        Dist2Attraction=MODA_distance(X(:,i),Food_pos(:,1));
        if all(Dist2Attraction<=r)
            F=Food_pos-X(:,i);
            iter;
        else
            F=0;
        end
        
        % Distraction from enemy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.5)
        Dist=MODA_distance(X(:,i),Enemy_pos(:,1));
        if all(Dist<=r)
            E=Enemy_pos+X(:,i);
        else
            E=zeros(dim,1);
        end
        
        for tt=1:dim
            if X(tt,i)>ub(tt)
                X(tt,i)=lb(tt);
                DeltaX(tt,i)=rand;
            end
            if X(tt,i)<lb(tt)
                X(tt,i)=ub(tt);
                DeltaX(tt,i)=rand;
            end
        end
        
        
        Archive_BS_num=randi(Archive_member_no);
        Archive_BS=Archive_X(Archive_BS_num,dim-2*cluster_num+1: dim);
        
        
        X_BS=X(dim-2*cluster_num+1: dim,i);
        
        Food_BS=Food_pos(dim-2*cluster_num+1: dim);
        
        if any(Dist2Attraction>r)
            if neighbours_no>1
                for j=1:dim
                    DeltaX(j,i)=w*DeltaX(j,i)+rand*A(j,1)+rand*C(j,1)+rand*S(j,1);
                    if DeltaX(j,i)>V_max
                        DeltaX(j,i)=V_max;
                    end
                    if DeltaX(j,i)<-V_max
                        DeltaX(j,i)=-V_max;
                    end
                    X(j,i)=X(j,i)+DeltaX(j,i);
                end
                
            else
                X(:,i)=X(:,i)+MODA_Levy(dim)'.*X(:,i);
                DeltaX(:,i)=0;
            end
        else
            
        for j=1:dim-bs_num
                
                DeltaX(j,i)=s*S(j,1)+a*A(j,1)+c*C(j,1)+f*F(j,1)+e*E(j,1) + w*DeltaX(j,i);
                if DeltaX(j,i)>V_max
                    DeltaX(j,i)=V_max;
                end
                if DeltaX(j,i)<-V_max
                    DeltaX(j,i)=-V_max;
                end
                X(j,i)=X(j,i)+DeltaX(j,i);
            end
        end
        
        X(dim-2*cluster_num+1:dim,i)= MODA_updateOrder(Food_BS,X_BS,Archive_BS,bs_num, cluster_num);
        X = improved_MFO1( X , iter, max_iter, Food_pos, uav_num, BSs, dim, uav_origin,cluster_num, bs_num);
        
        Flag4ub=X(:,i)>ub';
        Flag4lb=X(:,i)<lb';  
        X(:,i)=(X(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
       
    end
    
    disp(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
end

figure(1)

str=datestr(now,30);

solution=Archive_X;
save(['__SOLUTION_',str],'solution')

moda=Archive_F;
save(['__PS_',str],'moda')

fileID = fopen('__RESULTS.txt','a');  % w:write，覆盖原有内容  a:add，在原有内容后面增加
fprintf(fileID,[num2str(max(moda(:,1))),'~',num2str(min(moda(:,1))),'\n']);
fprintf(fileID,[num2str(max(moda(:,2))),'~',num2str(min(moda(:,2))),'\n']);
fprintf(fileID,[num2str(max(moda(:,3))),'~',num2str(min(moda(:,3))),'\n\n']);
fclose(fileID);

plot3(Archive_F(:,1),Archive_F(:,2),Archive_F(:,3),'ko','MarkerSize',8,'markerfacecolor','k');
title('IMODAOM');
grid on;
xlabel('f1');
ylabel('f2');
zlabel('f3');

saveas(gcf,['__FIG_',str,'.fig']);

