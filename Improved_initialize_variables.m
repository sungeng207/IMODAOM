
function f = Improved_initialize_variables(N, V, min_range, max_range, bs_num,uav_num,cluster_num, BSs, Evsdp, uav_origin)
 
    q=4;
    j=3;
    n=(q^j-1)/(q-1);
    A=improve_OA_permut(q,n,j);
    A2=q*ones(64,21);
    A2=A2-A;
    A3=[A,A2];
    A3=A3/4;
    A4=[A3;A3;A3;A3];

    V=V-2*cluster_num;
    min = min_range;
    max = max_range;
    
    for i = 1 : N
        for j = 1 : V
            f(i,j) = min(j) + (max(j) - min(j))* A4(j,i);
        end
        f(i,V+1:V+cluster_num) = round(rand(1,cluster_num)*(bs_num-1) + 1);
        f(i,V+cluster_num+1:V+cluster_num+cluster_num) = randperm(4);
    end

end
