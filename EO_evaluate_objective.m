function f = EO_evaluate_objective(uavs, V,  uav_num, bs_num, cluster_num, BSs, Evsdp, uav_origin)
    
    BS_selection=uavs(V-2*cluster_num+1:V-cluster_num);
    BS_index=uavs(V-cluster_num+1:V);
    
    f=zeros(1,3);
    
    for i=1:cluster_num
        
        
        uav_access=BS_index(i);
        num1=4*uav_num*(uav_access-1)+1;
        num2=4*uav_num*uav_access;
        uav=uavs( num1 : num2 );
        BS=BSs((uav_access-1)*bs_num+BS_selection(uav_access),:);
        [O1, O2]= Objective12_secrecy_SLL(uav, uav_num, BS, Evsdp);
        
        f(1)=f(1)-O1;
        f(2)=f(2)+O2;
        f(3) = f(3)+Objective3_energy(uav, uav_origin, uav_num);
        
        uav_origin=uav;
        
    end
    
end


function [f1, f2] = Objective12_secrecy_SLL(uavs, uav_num, user, Evsdp)

    f=0.9E9;                      
    lambda = 3E8/f;               
    k = 2*pi/lambda;             
    I0=1;                          
    B=20E6;                        
    P_uav=0.1;                
    Ko=(lambda/4/pi)^2;             
    sigma2=(10^-18.5)*B;
    alpha=2.7;                      
    
    %%   obtain coordinate and angel
    
    uav_center=[mean(uavs(uav_num+1:2*uav_num)),mean(uavs(2*uav_num+1:3*uav_num)),mean(uavs(3*uav_num+1:4*uav_num)) ];
    for n=1:uav_num
        uavs(uav_num+n)=uavs(uav_num+n)-uav_center(1);
        uavs(2*uav_num+n)=uavs(2*uav_num+n)-uav_center(2);
        uavs(3*uav_num+n)=uavs(3*uav_num+n)-uav_center(3); 
    end   
    
    uav_user = [user(1)-uav_center(1),user(2)-uav_center(2),user(3)-uav_center(3)];    
    [user_an_1,user_an_2]=obtain_angel(uav_user);
    user_an_1=user_an_1/180*pi; 
    user_an_2=user_an_2*pi/180;
    
    uav_ed=zeros(4,3);
    Evsdp_angel=zeros(4,2);
    for i=1:4
       uav_ed(i,:)= [Evsdp(i,1)-uav_center(1), Evsdp(i,2)-uav_center(2), Evsdp(i,3)-uav_center(3)]; 
       [Evsdp_angel(i,1), Evsdp_angel(i,2)]=obtain_angel(uav_ed(i,:));
       Evsdp_angel(i,1)=Evsdp_angel(i,1)/180*pi; 
       Evsdp_angel(i,2)=Evsdp_angel(i,2)/180*pi; 
    end

    
    %% calculate AF, G, R
    Ntheta=181; Nphi=361;         
    theta = linspace(0,pi,Ntheta);      %theta
    phi = linspace(-pi,pi,Nphi);        %phi
    [Phi,Theta] = meshgrid(phi,theta); 
    CT=cos(Theta); ST=sin(Theta); CP=cos(Phi); SP=sin(Phi);
   
    Rx=ST.*CP; Ry=ST.*SP; Rz=CT;        %unit radial vector components
    phi_inc=user_an_2;              % incident signal of interest 
    theta_inc=user_an_1;
  
    In=uavs(1,1:uav_num);
    Pi=In.^2*P_uav;
	dis_x=uavs(1,uav_num+1:uav_num*2);
	dis_y=uavs(1,uav_num*2+1:uav_num*3);
	dis_z=uavs(1,uav_num*3+1:uav_num*4);
    
    AF_Uniform=zeros(Ntheta,Nphi);
    Ip= I0 * exp(-1i*k*(sin(theta_inc)*cos(phi_inc)*dis_x+sin(theta_inc)*sin(phi_inc)*dis_y+cos(theta_inc)*dis_z));
    Inp=Ip.*In;
    for n=1:uav_num
        AF_Uniform = AF_Uniform+Inp(n)*exp(1i*k*(dis_x(n)*Rx+dis_y(n)*Ry+dis_z(n)*Rz));
    end
    
    G_denominator= sum(sum( (abs(AF_Uniform).^2)*pi/180*2*pi/360.*ST));
    G_value=4*pi*(abs(AF_Uniform).^2)/G_denominator;

    h=abs(uav_center(3));
    d2d_user=sqrt(uav_user(1)^2+uav_user(2)^2);
    ri_user=sqrt(uav_user(1)^2+uav_user(2)^2+uav_user(3)^2);
    
    adpl_los=0.501;
    adpl_nlos=0.00501;
    
    d1_1=460*log10(h)-700;
    d1_2=18;
    d1=max(d1_1,d1_2);
    p1=4300*log10(h)-3800;
    
    P_los_user=d1/d2d_user+exp(-d2d_user/p1)*(1-d1/d2d_user);
    NP_los_user=1-P_los_user;
    pass_loss_user=P_los_user*adpl_los+NP_los_user*adpl_nlos; 
    
    angel_1=int32(user_an_1/pi*180)+1;
    angel_2=int32(user_an_2/pi*180)+181;
    R_user= B * log2(1+(ri_user^-alpha)*sum(Pi)*Ko*pass_loss_user*G_value(angel_1,angel_2)/sigma2);
   
    R_evsdp=zeros(1,4);
    
    for i=1:4
        
        angel_1=int32(Evsdp_angel(i,1)/pi*180)+1;
        angel_2=int32(Evsdp_angel(i,2)/pi*180)+181;
        
        d2d_evsdp=sqrt(uav_ed(i,1)^2+uav_ed(i,2)^2); 
        ri_evsdp=sqrt(uav_ed(i,1)^2+uav_ed(i,2)^2+uav_ed(i,3)^2);    
        
        P_los_evsdp=d1/d2d_evsdp+exp(-d2d_evsdp/p1)*(1-d1/d2d_evsdp);  
        NP_los_evsdp=1-P_los_evsdp;
        pass_loss_evsdp=P_los_evsdp*adpl_los+NP_los_evsdp*adpl_nlos;
        R_evsdp(i)=B * log2(1+((ri_evsdp^-alpha)*sum(Pi)*Ko*G_value(angel_1,angel_2)*pass_loss_evsdp/sigma2));
  
    end
 
    f1=R_user-max(R_evsdp);

    [ML, maxSLL] = cal_maxSLL(AF_Uniform);
    
    f2=20*log10( abs(maxSLL/ML ) );

end

function f = Objective3_energy(uav,uav_origin, uav_num)

    dis = zeros(1,uav_num);

    for i= 1: uav_num
        x1=uav(uav_num+i);
        x2=uav_origin(uav_num+i);
        y1=uav(2*uav_num+i);
        y2=uav_origin(2*uav_num+i);
        dis(i) = sqrt( (x1-x2)^2 + (y1-y2)^2 ); 
    end

    form_time=sum(dis)/uav_num/10.2;

    V_ad = dis/ form_time;

    Utip=120;
    v0=4.03;
    d0=0.6;
    s=0.05;
    rho=1.225;
    A=0.503;
    delta=0.012;
    omega=300;
    R=0.4;
    W=20;
    k=0.1;
    g=9.8;

    Po=delta/8*rho*s*A*(omega^3)*(R^3);
    Pi=(1+k)*sqrt(W^3)/sqrt(2*rho*A);
    Pv= Po* ( 1+3*(V_ad/Utip).^2 ) + Pi*( sqrt( sqrt( 1+ (V_ad.^4/4/v0.^4) ) - V_ad.^2/2/v0.^2 ) ) +0.5*d0*rho*s*A*V_ad.^3;

    E_xy=Pv*form_time;

    for i=1:uav_num
        if V_ad(i)<4
            E_xy(i) = (126 * dis(i) / 10.2) + (168 * (form_time- dis(i) / 10.2));
        end
    end

    uav1_z=uav(uav_num*3+1:uav_num*4);
    uav2_z=uav_origin(uav_num*3+1:uav_num*4);
    move_z=uav2_z-uav1_z;
    P_z(move_z>=0)=196.000;
    P_z(move_z<0)=139.750;
    E_z=sum(P_z.*abs(move_z));

    f = sum(E_xy) +E_z;

end


%obtain angel
function [user_an_1, user_an_2] = obtain_angel(user)

    uav_center=[0,0,0];
    
    molecule=sqrt( (uav_center(1)-user(1))^2+ (uav_center(2)-user(2))^2 ) ;
    denominator= abs(uav_center(3)-user(3));
    
    if user(3)<uav_center(3)
        user_an_1=pi-atan(molecule/denominator);
    else
        user_an_1=atan(molecule/denominator);
    end
    
    side1=user(2)-uav_center(2);
    side2=user(1)-uav_center(1) ;
    
    if side1>0 && side2>0
        user_an_2=atan(abs(side1)/abs(side2) );
    elseif side1<0&& side2>0
        user_an_2=pi-atan(abs(side1)/abs(side2) );
    elseif side1<0&& side2<0    
        user_an_2=-pi+atan(abs(side1)/abs(side2) );
    else
        user_an_2=-atan(abs(side1)/abs(side2) );
    end

    user_an_1=round(user_an_1/pi*180);
    user_an_2=round(user_an_2/pi*180);

end

% calculate SLL for 3D
function [ML, maxSLL] = cal_maxSLL( R )

    if size(R,1)>1
        R=max(R);
    end

    [~,R_length]=size(R);
    crest=[];
    for n=1:R_length

        if n-1<1
            former=R_length;
        else
            former=n-1;
        end

        if n+1>R_length
            latter=1;
        else
            latter=n+1;
        end

        if R(n)>R(former) && R(n)>R(latter)       
            crest(size(crest,2)+1)=R(n);
        end

    end

    crest_sort=sort(crest,'descend');
    
    if size(crest_sort,2)>0
    
        ML=crest_sort(1);
        if size(crest_sort,2)>1
            maxSLL=crest_sort(2);
        else 
            maxSLL=crest_sort(1);
        end
        
    else 
        ML=0;
        maxSLL=0;
    end

end


 



