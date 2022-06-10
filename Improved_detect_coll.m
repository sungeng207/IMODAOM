function new_solution = Improved_detect_coll( solution, dim, uav_num, cluster_num )

    uavs=solution(1:dim-2*cluster_num);
    for i=1:cluster_num
     
        uav=uavs((i-1)*uav_num*4+1:i*4*uav_num);
        uav=handle_NaN(uav,uav_num);
        flag=detect_collision(uav,uav_num);
        while sum(flag)>0
            uav=handle_collision(uav,uav_num,flag);
            flag=detect_collision(uav,uav_num);
        end
        uavs((i-1)*uav_num*4+1:i*4*uav_num)=uav;
    end
    solution(1:dim-2*cluster_num)=uavs;
    new_solution=solution;
    
end

function f = detect_collision(uav,uav_num)
    
    f=zeros(1,uav_num);
    for i=1:uav_num
        for j = i+1:uav_num
            x1=uav(uav_num+i);
            x2=uav(uav_num+j);
            y1=uav(2*uav_num+i);
            y2=uav(2*uav_num+j);
            z1=uav(3*uav_num+i);
            z2=uav(3*uav_num+j);  
            if (sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 )<0.5 )
%             if (sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 )< lambda/2 )
                f(i)=1;
            end
        end
    end
end


function f = handle_collision(uav,uav_num,flag)

    improved_rate=1/sqrt(100^2+100^2+30^2);
    
    uav_x=uav(uav_num+1:uav_num*2);
    uav_y=uav(uav_num*2+1:uav_num*3);
    uav_z=uav(uav_num*3+1:uav_num*4);
    
    for i=1:uav_num
        
        if flag(i)==1
            dis2=(uav_x- uav(uav_num+i)*ones(1,uav_num)).^2 +(uav_y- uav(2*uav_num+i)*ones(1,uav_num)).^2+(uav_z- uav(3*uav_num+i)*ones(1,uav_num)).^2;
            [~,id]=sort(dis2,'ascend'); 
            uav(uav_num+i)=uav(uav_num+i)- improved_rate*rand()*uav(uav_num+id(2))-improved_rate*rand()*uav(uav_num+id(3))-improved_rate*rand()*uav(uav_num+id(4));
        end
    end
    
    flag=detect_collision(uav,uav_num);
    
    if sum(flag)>0
        for i=1:2*uav_num
            xy(i)=100*rand();
        end

        for i=2*uav_num+1:3*uav_num
            xy(i)=60+30*rand();
        end
        uav(uav_num+1:4*uav_num)=xy;
    end
    
    f=uav;
    
end

function f=handle_NaN(uav,uav_num)
    if sum(uav(1:uav_num))==0
        uav(1:uav_num)=1;
    end
    f=uav;
end