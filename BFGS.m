function [tau1, tau2, x, psi, t, Q,u0] = BFGS(tau1_0, tau2_0, h0, u0, B, g, l, a_max, x0, k, xf, T)

%STEP 1 - initial conditions
ep1 = 1e-7;
ep2 = 1e-7;
ep3 = 1e-7;
tau1 = tau1_0;
tau2 = tau2_0;
tau1_s = tau1;
tau2_s = tau2;
Qx = q_cost_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0, k, T);
R = 1;
gradQ_s = [tau1, tau2];

iteration = 0;
spike_generated=0;
while iteration < 1000
    iteration = iteration + 1;
    %STEP 2 - first STOP condition
    
    [gradQ, x, psi, t, Q, cn1, cn2] = gradient(tau1, tau2, h0, u0, B, g, l, a_max, x0, xf, k, T);
    if gradQ*gradQ' <= ep1,
        % calculating switching function
        u=zeros(length(t),2);
        u(1:cn1(1)-1,1)=u0(1);
        u(1:cn2(1)-1,2)=u0(2);
        for iu=1:length(cn1)-2,
            u(cn1(iu):cn1(iu+1)-1,1)=-u(cn1(iu)-1);
        end
        for iu=1:length(cn2)-2,
            u(cn2(iu):cn2(iu+1)-1,2)=-u(cn2(iu)-1);
        end
        u(cn1(end-2):end,1)=-u(cn1(end-2)-1);
        u(cn2(end-2):end,2)=-u(cn2(end-2)-1);
        
        
        [psi4_max, psi4_max_it]=max(abs(psi(psi(:,4).*u(:,1)<0,4)));
        t_psi4_max=t(psi(:,4).*u(:,1)<0);
        [psi8_max, psi8_max_it]=max(abs(psi(psi(:,8).*u(:,2)<0,8)));
        t_psi8_max=t(psi(:,8).*u(:,2)<0);
        
        % generating spikes
        mutiplier=10;
        if ~isempty(psi4_max)
            if(psi4_max==psi(end,4) && tau1(end)<T-mutiplier*ep3)
                tau1(end+1)=T-2*ep3;
                spike_generated=1;
                
            elseif(psi4_max==psi(1,4) && tau1(1)>mutiplier*ep3)
                tau1=[ep3*2, tau1];
                spike_generated=1;
                u0(1)=-u0(1);
                elseif(t_psi4_max(psi4_max_it)>mutiplier*ep3 && t_psi4_max(psi4_max_it)<T-mutiplier*ep3)
                [min_t_tau1 min_t_tau1_it]=min(abs(tau1- t_psi4_max(psi4_max_it)));
                if min_t_tau1>mutiplier*ep3,
                    tau1(end+1)=t_psi4_max(psi4_max_it)-ep3;
                    tau1(end+1)=t_psi4_max(psi4_max_it)+ep3;
                    tau1=unique(tau1);
                    spike_generated=1;
                end
            end
        end
        if ~isempty(psi8_max)
            if(psi8_max==psi(end,8) && tau2(end) < T-mutiplier*ep3)
                tau2(end+1)=T-2*ep3;
                spike_generated=1;
            elseif(psi8_max==psi(1,8) && tau2(1)>mutiplier*ep3)
                tau2=[ep3*2, tau2];
                u0(2)=-u0(2);
                spike_generated=1;
            elseif(t_psi8_max(psi8_max_it)>mutiplier*ep3 && t_psi8_max(psi8_max_it)<T-mutiplier*ep3)
                [min_t_tau2 min_t_tau2_it]=min(abs(tau2- t_psi8_max(psi8_max_it)));
                if min_t_tau2>mutiplier*ep3,
                    tau2(end+1)=t_psi8_max(psi8_max_it)-ep3;
                    tau2(end+1)=t_psi8_max(psi8_max_it)+ep3;
                    tau2=unique(tau2);
                    spike_generated=1;
                end
            end
        end
        if(spike_generated)
            %             disp('spike')
            R=1;
            spike_generated=0;
            continue;
        end
        % stop condition
        %         disp('STEP2 stop');
        break;
    end
    %STEP 3 - setting direction as vector d
    if R == 1
        W = eye(length(tau1)+length(tau2));
    else
        s = [tau1, tau2] - [tau1_s, tau2_s];
        s=s';
        r = gradQ - gradQ_s;
        r=r';
        W = W+(r*r')/(s'*r)-(W*(s*s')*W)/(s'*W*s);
    end
    d = -W^(-1)*gradQ';
    d=d';
    %STEP 4 - check if direction is worth searching on
    if d'*gradQ >= 0
        R = 1;
        continue;
    end
    %STEP 5 - saving old previous values
    tau1_s = tau1;
    tau2_s = tau2;
    gradQ_s = gradQ;
    max_contraction = 10;
    d1=[0, d(1:length(tau1)), 0];
    d2=[0, d(length(tau1)+1:end), 0];
    dd1=diff(d1);
    dd2=diff(d2);
    dtau1=diff([0 tau1 T]);
    dtau2=diff([0 tau2 T]);
    lambda1 =min([-dtau1(dd1<0)./dd1(dd1<0), 1]);
    lambda2 =min([-dtau2(dd2<0)./dd2(dd2<0), 1]);
    lambda = min([lambda1, lambda2]);
    while max_contraction > 0,% && Qn > Qx
        max_contraction = max_contraction-1;
        tau1_ = tau1+lambda*d(1:length(tau1));
        tau2_ = tau2+lambda*d((length(tau1)+1):(length(tau1)+length(tau2)));
        Qn = q_cost_BB(h0, tau1_, tau2_, u0, B, g, l, a_max, x0, k, T);
        if(Qn < Qx)
            break;
        end
        lambda = lambda / 2;
    end
    
    tau1 = tau1_;
    tau2 = tau2_;
    
    % reducing tau length
    changed=0;
    dtau1=diff([0 tau1 T]);
    dtau2=diff([0 tau2 T]);
    [min_dtau1, min_dtau_it1]=min(abs(dtau1));
    [min_dtau2, min_dtau_it2]=min(abs(dtau2));
    if(min_dtau1<ep3)
        if(min_dtau_it1==length(dtau1))
            tau1(end)=[];
        elseif(min_dtau_it1==1)
            tau1(1)=[];
            u0(1)=-u0(1);
        else
            tau1((min_dtau_it1-1):min_dtau_it1)=[];
        end
        changed=1;
    end
    if(min_dtau2<ep3)
        if(min_dtau_it2==length(dtau2))
            tau2(end)=[];
        elseif(min_dtau_it2==1)
            tau2(1)=[];
            u0(2)=-u0(2);
        else
            tau2((min_dtau_it2-1):min_dtau_it2)=[];
        end
        changed=1;
    end
    if(changed==1),
        %         disp('reduction')
        R=1;
        continue;
    end
    Qx = Qn;
    % STEP 7 - checking STOP conditions
    if norm([tau1, tau2]-[tau1_s, tau2_s]) < ep2
        if R == 1
            %             disp('STEP7 stop');
            break;
        else
            R = 1;
            continue;
        end
    else
        R = 0;
    end
end
Q=Qx;
