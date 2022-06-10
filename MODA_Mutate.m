%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function X=MODA_Mutate(Food_X,X,uav_num,bs_num)

    solution_num=4*uav_num;
    
    for n=1:bs_num
        
        uav_location=Food_X(1, (n-1)*solution_num+1:n*solution_num );
       
        
        for m=1: solution_num
                
            if m<uav_num+1
                X( (n-1)*uav_num*4+m)= uav_location(m)+rand*(X( (n-1)*uav_num*4+m)-uav_location(m));
            else
                X( (n-1)*uav_num*4+m)= uav_location(m);
            end
             
        end
                    

end