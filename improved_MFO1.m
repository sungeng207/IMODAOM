function new_X = improved_MFO1( X , iter, max_iter, Elite_position, uav_num, BSs, V, UAV_origin,cluster_num, bs_num)


    BS_selection=X(V-2*cluster_num+1:V-cluster_num);
    BS_index=X(V-cluster_num+1:V);
    

    for k=1: cluster_num
        
        uav_access=BS_index(k);
        
        num1=4*uav_num*(uav_access-1)+1;
        num2=4*uav_num*uav_access;
        
        uav_location=X(num1 : num2);
        
        BS=BSs((uav_access-1)*bs_num+BS_selection(uav_access),:);
        
        
        uav_center_x=mean(uav_location(uav_num+1:uav_num*2) );
        uav_center_y=mean(uav_location(uav_num*2+1:uav_num*3));
        uav_center_z=mean(uav_location(uav_num*3+1:uav_num*4));

        if 1==1

            for i=1:uav_num

                uav_x=uav_location(i+uav_num);
                uav_y=uav_location(i+2*uav_num);
                uav_z=uav_location(i+3*uav_num);

                if uav_x > 100
                    uav_x=100;
                elseif uav_x<0
                    uav_x=0;
                end

                if uav_y > 100
                    uav_y=100;
                elseif uav_y<0
                    uav_y=0;
                end

                if uav_z > 90
                    uav_z=90;
                elseif uav_z<60
                    uav_z=60;
                end

                improve_rate1=1/sqrt(BS(1)^2+BS(2)^2+uav_z^2);
                improve_rate2=1/sqrt(100^2+100^2+30^2);
                improve_rate3=1/sqrt(100^2+100^2+30^2);

                uav_location(i+uav_num) = uav_x+ improve_rate1*(rand() * (BS(1)- uav_x)) + improve_rate2*(rand() * (UAV_origin(i+uav_num)- uav_x))+ improve_rate3*(rand()*(uav_center_x-uav_x));
                uav_location(i+2*uav_num) = uav_y+ improve_rate1*(rand() *(BS(2)- uav_y)) + improve_rate2*(rand() * (UAV_origin(i+2*uav_num)- uav_y))+improve_rate3*(rand()*(uav_center_y-uav_y));
                uav_location(i+3*uav_num) = uav_z- improve_rate1*(rand() * (BS(3)- uav_z)) + improve_rate2*(rand() * (UAV_origin(i+3*uav_num)- uav_z))+improve_rate3*(rand()*(uav_center_z-uav_z));


            end

        end
        
        X(num1 : num2)=uav_location;
        UAV_origin=uav_location;
    
    end 
    
    new_X=X;

end
    
