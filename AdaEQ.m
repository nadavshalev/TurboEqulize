classdef AdaEQ < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        P;
        Q;
        
        Ncheck;
        Pcheck;
        Qcheck;
        
        Lr; % length of channel ISI
        Ld;
        Mu; % step size
        MaxErr;
        overSamp;
        
        Statistics;
        tmpData;
    end
    
    methods
        function obj = AdaEQ(Lr, Ld, dMu, overSamp, checkOut)
            
            % FF params
            obj.P = zeros(2*Lr+1,1);
            obj.P(Lr+1) = 1; % t0 symb in the middle
            
            % FB params
            obj.Q = zeros(2*Ld+1,1);
            obj.Q(Ld+1) = 0;  % t0 symb in the middle
            
            % save cahnel est length
            obj.Lr = Lr;
            obj.Ld = Ld;
            
            obj.Mu = dMu;
            obj.overSamp = overSamp;
            
            % checkout filter params iter number
            if nargin < 5
                obj.Ncheck = -1;
            end
            
        end
        
        function eqOut = turboEqualize_train(obj, chanOut,estOut,packLen)
            
            N = packLen;
           
            % padd buffers
            r = [zeros(obj.Lr,1);chanOut;zeros(obj.Lr,1);]; % FF buff
            d = [zeros(obj.Ld,1);estOut;zeros(obj.Ld,1);]; % FB buff
            
            eqOut = zeros(N,1);
            
            for i=1:N
                
                % set buff index
                ind = i+obj.Ld;
                indR = obj.overSamp * (i-1) + obj.Lr+1;
                
                % cut buff part - for conv with params
                R = flip(r(indR-obj.Lr:indR+obj.Lr));
                D = flip(d(ind-obj.Ld:ind+obj.Ld));
                
                % conv with params
                sn_ = obj.P.' * R - obj.Q.' * D;
                
                % calc error
                e = sn_-getHard(estOut(i));
                
                % update params
                obj.P = obj.P - obj.Mu*conj(R)*e;
                obj.Q = obj.Q + obj.Mu*conj(D)*e;
                obj.Q(obj.Ld+1) = 0; % make sure t0 is 0
                
                % get (soft) result
                eqOut(i) = sn_;
                
            end
            
            % checkout P and Q
            obj.Pcheck = obj.P;
            obj.Qcheck = obj.Q;
        end
        
        function eqOut = turboEqualize(obj, chanOut,estOut,packLen)     
            
            N = packLen;
            
            % get P and Q from checkout
            obj.P = obj.Pcheck;
            obj.Q = obj.Qcheck;
           
            % padd buffers
            r = [chanOut;zeros(obj.Lr,1);];
            d = [estOut;zeros(obj.Ld,1);];
            
            eqOut = zeros(N,1);
            
            for i=1:N
                
                % set buff index
                ind = i+obj.Ld;
                indR = obj.overSamp * (i-1) + obj.Lr+1;
                
                % cut buff part - for conv with params
                R = flip(r(indR-obj.Lr:indR+obj.Lr));
                D = flip(d(ind-obj.Ld:ind+obj.Ld));
                
                % conv with params
                sn_ = obj.P.' * R - obj.Q.' * D;
                
                % calc error (using hard of estimated prior)
                e = sn_-getHard(d(ind));
                
                % update params
                obj.P = obj.P - obj.Mu*conj(R)*e;
                obj.Q = obj.Q + obj.Mu*conj(D)*e;
                obj.Q(obj.Ld+1) = 0;
                
                % get (soft) result
                eqOut(i) = sn_;
                
                % chenge prior to est sn_
                d(ind) = getHard(sn_);
                
                % checkout P and Q
                if i == obj.Ncheck
                    obj.Pcheck = obj.P;
                    obj.Qcheck = obj.Q;
                end
                
            end
        end
        

    end
end

