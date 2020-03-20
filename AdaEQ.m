classdef AdaEQ < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        P;
        Q;
        
        L; % length of channel ISI
        Mu; % step size
        MaxErr;
        
        Statistics;
        tmpData;
    end
    
    methods
        function obj = AdaEQ(L, dMu)
            obj.P = zeros(2*L+1,1);
            obj.P(L+1) = 1;
            obj.Q = zeros(2*L+1,1);
            obj.Q(L+1) = 0;
            obj.L = L;
            obj.Mu = dMu;
            obj.Statistics = {};
            obj.Statistics.E = [];
            obj.Statistics.sn = [];
            obj.Statistics.symbHard = [];
            obj.Statistics.P = [];
            obj.Statistics.Q = [];
            
        end
        
%         e is hard
        function eqOut = turboEqualize(obj, chanOut,estOut,mu)
             if nargin < 4
                mu = obj.Mu;
             end
            
            N = length(estOut)-obj.L;
           
            r = [chanOut;zeros(obj.L,1);];
            d = [estOut;zeros(obj.L,1);];
            obj.Q(obj.L+1) = 0;
            eqOut = zeros(N,1);
            
            err = zeros(N,1);
            sn = zeros(N,1);
            symbHard = zeros(N,1);
            
            for i=1:N
                ind = i+obj.L;
                R = flip(r(ind-obj.L:ind+obj.L));
                D = flip(d(ind-obj.L:ind+obj.L));
                
                sn_ = obj.P.' * R - obj.Q.' * D;
                e = sn_-getHard(d(ind));
                err(i) = e;
                sn(i) = sn_;
                symbHard(i) = getHard(estOut(i));
                
                obj.P = obj.P - mu*conj(R)*e;
                obj.Q = obj.Q + mu*conj(D)*e;
                obj.Q(obj.L+1) = 0;

                eqOut(i) = sn_;
                d(ind) = getHard(sn_);
                
            end
            obj.Statistics.E = [obj.Statistics.E;err];
            obj.Statistics.sn = [obj.Statistics.sn;sn];
            obj.Statistics.symbHard = [obj.Statistics.symbHard;symbHard];
            obj.Statistics.P = [obj.Statistics.P obj.P];
            obj.Statistics.Q = [obj.Statistics.Q obj.Q];
        end
        
        function eqOut = turboEqualize_train(obj, chanOut,estOut,mu)
             if nargin < 4
                mu = obj.Mu;
             end
            
            N = length(estOut);
           
            r = [zeros(obj.L,1);chanOut;zeros(obj.L,1);];
            d = [zeros(obj.L,1);estOut;zeros(obj.L,1);];
            obj.Q(obj.L+1) = 0;
            eqOut = zeros(N,1);
            
            err = zeros(N,1);
            sn = zeros(N,1);
            symbHard = zeros(N,1);
            
            for i=1:N
                ind = i+obj.L;
                R = flip(r(ind-obj.L:ind+obj.L));
                D = flip(d(ind-obj.L:ind+obj.L));
                
                sn_ = obj.P.' * R - obj.Q.' * D;
                e = sn_-getHard(estOut(i));
                err(i) = e;
                sn(i) = sn_;
                symbHard(i) = getHard(estOut(i));
                
                obj.P = obj.P - mu*conj(R)*e;
                obj.Q = obj.Q + mu*conj(D)*e;
                obj.Q(obj.L+1) = 0;

                eqOut(i) = sn_;
%                 d(ind) = getHard(sn_);
                
            end
            obj.Statistics.E = [obj.Statistics.E;err];
            obj.Statistics.sn = [obj.Statistics.sn;sn];
            obj.Statistics.symbHard = [obj.Statistics.symbHard;symbHard];
            obj.Statistics.P = [obj.Statistics.P obj.P];
            obj.Statistics.Q = [obj.Statistics.Q obj.Q];
        end
        
        
        
        function eqOut = dfeEqualize(obj, chanOut, mu)
             if nargin < 4
                mu = obj.Mu;
             end
            
            N = length(chanOut)-obj.L - (obj.L-1);
           
            r = [chanOut;zeros(obj.L,1);];
            d = [getHard(chanOut);zeros(obj.L,1);];
            obj.Q(obj.L+1) = 0;
            eqOut = zeros(N,1);
            
            err = zeros(N,1);
            sn = zeros(N,1);
            symbHard = zeros(N,1);
            
            for i=1:N
                ind = i+obj.L;
                R = flip(r(ind-obj.L:ind+obj.L));
                D = flip(d(ind-obj.L:ind+obj.L));
                
                sn_ = obj.P.' * R - obj.Q.' * D;
                e = sn_- getHard(sn_);
                err(i) = e;
                sn(i) = sn_;
                symbHard(i) = getHard(d(ind));
                
                obj.P = obj.P - mu*conj(R)*e;
                obj.Q = obj.Q + mu*conj(D)*e;
                obj.Q(obj.L+1) = 0;

                eqOut(i) = sn_;
                d(ind) = getHard(sn_);
                
            end
            obj.Statistics.E = [obj.Statistics.E;err];
            obj.Statistics.sn = [obj.Statistics.sn;sn];
            obj.Statistics.symbHard = [obj.Statistics.symbHard;symbHard];
            obj.Statistics.P = [obj.Statistics.P obj.P];
            obj.Statistics.Q = [obj.Statistics.Q obj.Q];
        end
        
        
        
        
        function eqOut = normalEqualize_train(obj, chanOut,estOut,mu)
            if nargin < 4
                mu = obj.Mu;
            end
             
            N = length(estOut);
           
            r = [zeros(obj.L,1);chanOut;zeros(obj.L,1);]; 
            d = [zeros(obj.L,1);estOut;zeros(obj.L,1);];
            obj.Q(obj.L+1) = 0;
            eqOut = zeros(N,1);
            for i=1:N
                ind = i+obj.L;
                R = flip(r(ind-obj.L:ind+obj.L));
                D = flip(d(ind-obj.L:ind+obj.L));
                
                sn = obj.P.' * R - obj.Q.' * D;
                
                e = sn-estOut(i);
                
                obj.P = obj.P - mu*conj(R)*e;
                obj.Q = obj.Q + mu*conj(D)*e;
                obj.Q(obj.L+1) = 0;
                eqOut(i) = sn;
            end
        end
        function eqOut = normalEqualize_run(obj, chanOut,mu)
            if nargin < 4
                mu = obj.Mu;
            end
             
            N = length(chanOut)-(obj.L-1)-obj.L;
           
            r = [chanOut;zeros(obj.L,1);]; 
            d = [getHard(chanOut);zeros(obj.L,1);];
            obj.Q(obj.L+1) = 0;
            eqOut = zeros(N,1);
            for i=1:N
                ind = i+obj.L;
                R = flip(r(ind-obj.L:ind+obj.L));
                D = flip(d(ind-obj.L:ind+obj.L));
                
                sn = obj.P.' * R - obj.Q.' * D;
                
                e = sn-getHard(sn);
                
                obj.P = obj.P - mu*conj(R)*e;
                obj.Q = obj.Q + mu*conj(D)*e;
                obj.Q(obj.L+1) = 0;
                eqOut(i) = sn;
                d(ind) = getHard(sn);
            end
            
        end
        
    end
end

