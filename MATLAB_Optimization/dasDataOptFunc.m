function output = pyroOptFunc(V,Y1,optType) 


%% ga Optimization Initialization
if optType == 0 %Dual Rate
    
    options.PopulationSize = 100;
    options.Display = 'iter';
    options.FunctionTolerance = 1e-8;
    nvars = 4;
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    %    Test 1
    lbx = [1e2, 60e3,  1e5,   120e3]; % Row
    ubx = [1e6, 120e3, 1e12,  300e3];

    FUNKY = @TGAOPT;

    [X_OPT,fval,exitflag,output] = ga(FUNKY,nvars,A,b,Aeq,beq,lbx,ubx,[],options);
    
else %Single Rate
    
    options.PopulationSize = 50;
    options.Display = 'iter';
    options.FunctionTolerance = 1e-8;
    nvars = 2;
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Test 1:
%     lbx = [1e2, 60e3];  %40 w/ 1e-8
%     ubx = [1e6, 120e3];
% %     %Test 2:
%     lbx = [1,   60e3];  %50 w/ 1e-8
%     ubx = [1e6, 120e3];
% %     %Test 3:
%     lbx = [1,    60e3];  %50 w/ 1e-8
%     ubx = [1e12, 230e3];
%     %Test 4:
%     lbx = [10,    60e3];  %50 w/ 1e-8
%     ubx = [1e12, 230e3];
%     %Test 5:
    lbx = [1,    60e2];  %50 w/ 1e-8
    ubx = [1e12, 230e3];

    FUNKY = @TGAOPT;

    [X_OPT,fval,exitflag,output] = ga(FUNKY,nvars,A,b,Aeq,beq,lbx,ubx,[],options);
end

%% Alpha Function 
    function out = alphaFunc(Mass_,mLength)

        for I = [1:1:mLength]

              f_(I,1) = 1 - (Mass_(1) - Mass_(I))/Mass_(1);
        end

        out = f_;
    end

%% F Model Function 
    function output = findingF(X)
        
        % Data Import
        master = dasData();
        T    = master{1,1}{1,V}; 
        dt   = master{1,3}{1,V};
        Mg   = master{1,4}{1,1};
        R    = 8.314;
        
        % Empty alpha matrix
        mLength = length(T);  
     
        % Complete devolatilisation check
        if optType == 0 %Singular Devol.
            M_ = Mg*Y1;
            Mass_ = zeros(mLength,1);
            Mass_(1,1) = M_;
            check = 0;
            
            for I = [2:1:mLength]

                dv_ = (Y1.*X(1).*exp(-X(2)./(R.*T(I))) + X(3).*exp(-X(4)./(R.*T(I)))).*M_.*dt(I);

                if dv_ < M_;

                    M_ = M_ - dv_;

                else dv_ > M_;

                    if check == 1;

                        M_ = 0;

                    else 

                        M_ = M_;  

                        check = 1;

                    end
                end
                
                Mass_(I,1) = M_;
            end            
       
        else
            M_ = Mg*Y1;
            Mass_ = zeros(mLength,1);
            Mass_(1,1) = M_;
            check = 0;
            
            for I = [2:1:mLength]

                dv_ = X(1).*exp(-X(2)./(R.*T(I))).*M_.*dt(I); 

                if dv_ < M_;

                    M_ = M_ - dv_;

                else dv_ > M_;

                    if check == 1;

                        M_ = 0;

                    else 

                        M_ = M_;  

                        check = 1;

                    end
                end
                Mass_(I,1) = M_;
            end           
            
        end

        % Convert mass to Alpha
        f_ = alphaFunc(Mass_,mLength);

        output = f_;
    end

    %% Main Function 
    function opt_output = TGAOPT(X)

        f_ = findingF(X);
        
        %% Data Import
        master = dasData();
        T      = master{1,1}{1,V}; 
        dadt   = master{1,2}{1,V};
        R      = 8.314;

        %% Optimization 
        if optType == 0 
            %Dual Rate
            dadt_opt =  Y1.*X(1).*exp(-X(2)./(R.*T)).*f_  + X(3).*exp(-X(4)./(R.*T)).*f_;
        else
            %Single Rate
            dadt_opt =  X(1).*exp(-X(2)./(R.*T)).*f_;
        end

        RSS = sum((dadt - dadt_opt).^2);  

        opt_output = RSS;

    end
        
     
      finalF = findingF(X_OPT);
      
      output = {X_OPT, finalF};

end
