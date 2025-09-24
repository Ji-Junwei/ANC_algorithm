%% McFxLMS with arbitrary channel
% (J x K x M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Structure of McFxLMS
%
%                 M x J
%              +-----------+                       +   
% x(n) ---+--->|   P(z)    |--d(n)----------------> sum --+---> e(n)
%      J  |    +-----------+                          ^-   |  M
%         |                                           |    |
%         |        \ K x J                M x K     y(n)   |     
%         |    +-----------+    K     +-----------+   |    |
%         +--->|   W(z)    |--yc(n)-->|   S(z)    |---+    |
%         |    +-----------+          +-----------+        |
%         |            \                                   |
%         |             \----------------\                 |
%         |       M x K                   \                |
%         |    +-----------+          +-----------+        |
%         +--->|  Est_S(z) |--xf(n)-->|   LMS     |<-------+
%              +-----------+          +-----------+        
% 
%       x(n):  reference siganl/noise source
%       d(n):  disturbance signal
%       xf(n): filtered refence signal
%       yc(n): control signal
%       y(n):  anti noise
%       e(n):  error signal
%       P(Z):  Primary path
%       W(Z):  control filters
%       Est_S(Z): estimated Secondary path
%       S(Z):  Secondary path (real)
%      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef MultiChannelFxLMS

    properties
        Wc       % control filters K x J x wlen
        wlen     % length of control filter
        SecP     % secondary path M x K x slen
        slen     % length of secondary path
        Ref      % reference signal J x T
        Dis      % disturbance signal M x T
        yc       % control signal K x T
        y        % antinoise M x T
        Err      % error signal M x T
        J        % number of reference
        K        % number of secondary source
        M        % number of error
    end

    methods
        % initialization
        function obj =  MultiChannelFxLMS(wLen,Secondarypath,sLen,ref,dis,numref,numsec,numerr)
            T = size(ref,2);     % duration
            obj.J = numref;
            obj.K = numsec;
            obj.M = numerr;
            obj.Wc = zeros(numsec,numref,wLen);
            obj.wlen = wLen;
            obj.SecP = Secondarypath;
            obj.slen = sLen;
            obj.Ref = ref;
            obj.Dis = dis;
            obj.yc = zeros(numsec,T);
            obj.y = zeros(numerr,T);
            obj.Err = zeros(numerr,T);

        end

        % McFxLMS
        function [obj] = McFxLMS_controller(obj,muw)
            T = size(obj.Ref,2);

            xc = zeros(obj.J,obj.wlen);       % reference buffer for control filter
            ys = zeros(obj.K,obj.slen);       % control signal buffer for secondary path
            xf = zeros(obj.J,obj.slen);       % reference buffer for secondary path
            Xsf = zeros(obj.J,obj.K,obj.M,obj.wlen);  % filtered reference buffer
            
            for n = 1:T
                % reference buffer update
                xc = [obj.Ref(:,n) xc(:,1:end-1)];

                % control signal
                for kk = 1:obj.K
                    bb = reshape(obj.Wc(kk,:,:),[obj.J,obj.wlen]);
                    obj.yc(kk,n) = sum(sum(bb .* xc,2),1);
                end

                % anti noise
                ys = [obj.yc(:,n) ys(:,1:end-1)];
                for mm = 1:obj.M
                    cc = reshape(obj.SecP(mm,:,:),[obj.K,obj.slen]);
                    obj.y(mm,n) = sum(sum(cc .* ys,2),1);
                end

                % error
                obj.Err(:,n) = obj.Dis(:,n) - obj.y(:,n);

                % filtered reference
                xf = [obj.Ref(:,n) xf(:,1:end-1)];
                Xt = reshape(repmat(xf,(obj.K*obj.M),1),[obj.J obj.K obj.M obj.slen]);

                res = zeros(obj.J,obj.K,obj.M);
                for jj = 1:obj.J
                    dd = reshape(Xt(jj,:,:,:),[obj.K obj.M obj.slen]);
                    res(jj,:,:) = sum(dd .* permute(obj.SecP,[2 1 3]),3);
                end
                Xsf(:,:,:,2:end) = Xsf(:,:,:,1:end-1);
                Xsf(:,:,:,1) = res;

                % gradient and control filter update
                delta = zeros(obj.K,obj.J,obj.wlen);
                Xsf1 = permute(Xsf,[2 1 4 3]); % K x J x wlen x M
                for mm = 1:obj.M
                    delta = delta + obj.Err(mm,n) * Xsf1(:,:,:,mm);
                end
                obj.Wc = obj.Wc + muw * delta;
            end
    
        end


     end
end
