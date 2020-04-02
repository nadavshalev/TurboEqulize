classdef AdaEQ < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        P;
        Q;
        
        Lr; % length of channel ISI
        Ld;
        Mu; % step size
        MaxErr;
        overSamp;
        
        Statistics;
        tmpData;
    end
    
    methods
        function obj = AdaEQ(Lr, Ld, dMu, overSamp)
            obj.P = zeros(2*Lr+1,1);
            obj.P(Lr+1) = 1;
            obj.Q = zeros(2*Ld+1,1);
            obj.Q(Ld+1) = 0;
            obj.Lr = Lr;
            obj.Ld = Ld;
            obj.Mu = dMu;
            obj.overSamp = overSamp;
            obj.Statistics = {};
            obj.Statistics.E = [];
            obj.Statistics.sn = [];
            obj.Statistics.symbHard = [];
            obj.Statistics.P = [];
            obj.Statistics.Q = [];
            
        end
        
        function eqOut = turboEqualize(obj, chanOut,estOut,packLen)      
            N = packLen;
           
            r = [chanOut;zeros(obj.Lr,1);];
            d = [estOut;zeros(obj.Ld,1);];
            obj.Q(obj.Ld+1) = 0;
            eqOut = zeros(N,1);
            
            err = zeros(N,1);
            sn = zeros(N,1);
            symbHard = zeros(N,1);
            
            for i=1:N
                ind = i+obj.Ld;
                indR = obj.overSamp * (i-1) + obj.Lr+1;
                
                R = flip(r(indR-obj.Lr:indR+obj.Lr));
                D = flip(d(ind-obj.Ld:ind+obj.Ld));
                
                sn_ = obj.P.' * R - obj.Q.' * D;
                e = sn_-getHard(d(ind));
                err(i) = e;
                sn(i) = sn_;
                symbHard(i) = getHard(estOut(i));
                
                obj.P = obj.P - obj.Mu*conj(R)*e;
                obj.Q = obj.Q + obj.Mu*conj(D)*e;
                obj.Q(obj.Ld+1) = 0;

                eqOut(i) = sn_;
                d(ind) = getHard(sn_);
                
            end
            obj.Statistics.E = [obj.Statistics.E;err];
            obj.Statistics.sn = [obj.Statistics.sn;sn];
            obj.Statistics.symbHard = [obj.Statistics.symbHard;symbHard];
            obj.Statistics.P = [obj.Statistics.P obj.P];
            obj.Statistics.Q = [obj.Statistics.Q obj.Q];
        end
        
        function eqOut = turboEqualize_train(obj, chanOut,estOut,packLen)
            
            N = packLen;
           
            r = [zeros(obj.Lr,1);chanOut;zeros(obj.Lr,1);];
            d = [zeros(obj.Ld,1);estOut;zeros(obj.Ld,1);];
            obj.Q(obj.Ld+1) = 0;
            eqOut = zeros(N,1);
            
            err = zeros(N,1);
            sn = zeros(N,1);
            symbHard = zeros(N,1);
            
            for i=1:N
                ind = i+obj.Ld;
                indR = obj.overSamp * (i-1) + obj.Lr+1;
                
                R = flip(r(indR-obj.Lr:indR+obj.Lr));
                D = flip(d(ind-obj.Ld:ind+obj.Ld));
                
                sn_ = obj.P.' * R - obj.Q.' * D;
                e = sn_-getHard(estOut(i));
                err(i) = e;
                sn(i) = sn_;
                symbHard(i) = getHard(estOut(i));
                
                obj.P = obj.P - obj.Mu*conj(R)*e;
                obj.Q = obj.Q + obj.Mu*conj(D)*e;
                obj.Q(obj.Ld+1) = 0;

                eqOut(i) = sn_;
                
            end
            obj.Statistics.E = [obj.Statistics.E;err];
            obj.Statistics.sn = [obj.Statistics.sn;sn];
            obj.Statistics.symbHard = [obj.Statistics.symbHard;symbHard];
            obj.Statistics.P = [obj.Statistics.P obj.P];
            obj.Statistics.Q = [obj.Statistics.Q obj.Q];
        end
      
        function eqOut = normalEqualize_run(obj, chanOut,packLen)
             
            N = packLen;
           
            r = [chanOut;zeros(obj.Lr,1);]; 
            d = [getHard(chanOut(1:obj.overSamp:end));zeros(obj.Ld,1);];
            obj.Q(obj.Ld+1) = 0;
            eqOut = zeros(N,1);
            
            err = zeros(N,1);
            sn = zeros(N,1);
            symbHard = zeros(N,1);
            
            for i=1:N
                ind = i+obj.Ld;
                indR = obj.overSamp * (i-1) + obj.Lr+1;
                
                R = flip(r(indR-obj.Lr:indR+obj.Lr));
                D = flip(d(ind-obj.Ld:ind+obj.Ld));
                
                sn_ = obj.P.' * R - obj.Q.' * D;
                
                e = sn_-getHard(sn_);
                
                err(i) = e;
                sn(i) = sn_;
                symbHard(i) = getHard(chanOut(i));
                
                obj.P = obj.P - obj.Mu*conj(R)*e;
                obj.Q = obj.Q + obj.Mu*conj(D)*e;
                obj.Q(obj.Ld+1) = 0;
                eqOut(i) = sn_;
                d(ind) = getHard(sn_);
            end
            obj.Statistics.E = [obj.Statistics.E;err];
            obj.Statistics.sn = [obj.Statistics.sn;sn];
            obj.Statistics.symbHard = [obj.Statistics.symbHard;symbHard];
            obj.Statistics.P = [obj.Statistics.P obj.P];
            obj.Statistics.Q = [obj.Statistics.Q obj.Q];
        end

    end
end

