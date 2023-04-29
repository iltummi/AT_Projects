%% functions
function q_out =Q_method(acc,mag)
A1 = [0 0 -9.81]';
A2 =  [27.555000000000000  -2.416900000000000 -16.084900000000001]';

U = [A1 A2];
    B = [acc, mag];

    
    C = zeros(3,3);
    D = zeros(3,1);

    wbj = 1;
    for j =1:1:2
        bj = B(:,j);
        uj = U(:,j);
        C = C + wbj*bj*uj';
        D = D + wbj*cross(bj,uj);
    end

    K = [C+C'-trace(C)*eye(3), D;D', trace(C)];

    [V,D] = eig(K);
    max_D = D(1,1);
    max_V = V(:,1);
    for l=1:1:4
        if D(l,l) > max_D
            max_D = D(l,l);
            max_V = V(:,l);
        end
    end
    q_out = [max_V(end) max_V(1:3)'];
   

  

end