%% single channel FxLMS for ANC

classdef ANC_algorithm
    properties
        Wc           % control filter dim: wlen x 1
        wlen         % control filter length
        SecP         % secondary path dim: slen x 1
        slen         % secondary path length
        Dis          % disturbance signal
        Ref          % reference signal
        N            % duration
        yc           % control signal
    end

    methods
        % Initialization
        function obj = ANC_algorithm(wLen,sLen,SecondaryPath,Dis,Ref)
            obj.Wc = zeros(wLen,1);    % Initialize control filter weights
            obj.wlen = wLen;            % Set control filter length
            obj.SecP = SecondaryPath;   % Assign secondary path
            obj.slen = sLen;            % Set secondary path length
            obj.Dis = Dis;              % Assign disturbance signal
            obj.Ref = Ref;              % Assign reference signal
            obj.N = length(Dis);            % duration
            obj.yc = zeros(obj.N,1);        % Initialize control signal
        end

        % FxLMS algorithm
        function [e,obj] = ANC_FxLMS(obj,muw)
            e = zeros(obj.N,1);

            xc = zeros(obj.wlen,1);   % reference signal buffer; delay line for FIR control filter
            xs = zeros(obj.slen,1);   % reference signal buffer; delay line for FIR secondary path
            xf = zeros(obj.wlen,1);   % filtered reference signal vector for control filter updating
            ys = zeros(obj.slen,1);   % control signal buffer; delay line for FIR secondary path

            % iteration
            for i = 1 : obj.N
                xc = [obj.Ref(i);xc(1:end-1)];  % update reference vector
                obj.yc(i) = obj.Wc'*xc;         % generate control signal

                ys = [obj.yc(i);ys(1:end-1)];   % update control signal vector
                y = ys'*obj.SecP;               % anti-noise

                e(i) = obj.Dis(i) - y;          % residual error signal

                xs = [obj.Ref(i);xs(1:end-1)];  % update reference vector
                fx = xs'*obj.SecP;              % filtered reference signal
                xf = [fx;xf(1:end-1)];          % update filtered reference signal vector
                obj.Wc = obj.Wc + muw*xf*e(i);  % control filter update
            end
        end

        
    end

end