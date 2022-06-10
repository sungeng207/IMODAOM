function [ BS ] = MODA_updateOrder(Elite_BS,X_BS,Archive_BS,bs_num, cluster_num)

    for i=1: cluster_num
        if rand()>0.7
            BS(i)= Elite_BS(i);
        elseif rand()>0.3
            BS(i)= X_BS(i);
        else
            BS(i)= randi(bs_num);
        end
    end

    if rand()>0.7
        BS(cluster_num+1: 2*cluster_num)= tsp_crossover(Elite_BS(cluster_num+1: 2*cluster_num)', X_BS(cluster_num+1: 2*cluster_num)',1)';
    elseif rand()>0.3
        BS(cluster_num+1: 2*cluster_num)= X_BS(cluster_num+1: 2*cluster_num);
    else
        BS(cluster_num+1: 2*cluster_num)= randperm(cluster_num);
    end
    
end

function [childPath] = tsp_crossover(parent1Path, parent2Path, prob)

    random = rand();
    if prob >= random
        [l, length] = size(parent1Path);
        childPath = zeros(l,length);
        setSize = floor(length/2) -1;
        offset = randi(setSize);
        for i=offset:setSize+offset-1
            childPath(1,i) = parent1Path(1,i);
        end
        iterator = i+1;
        j = iterator;
        while any(childPath == 0)
            if j > length
                j = 1;
            end
            if iterator > length
               iterator = 1;
            end
            if ~any(childPath == parent2Path(1,j))
                childPath(1,iterator) = parent2Path(1,j);
                iterator = iterator + 1;
            end
            j = j + 1;
        end
        else
        childPath = parent1Path;
    end
end

