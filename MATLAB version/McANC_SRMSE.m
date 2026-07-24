%% Conventional Multichannel ANC (1 x K x M)  single reference, multiple secondary sources and error sensors (SRMSE)
% 1 reference, K secondary sources, M error sensors


classdef McANC_SRMSE
    properties
        Wc    % control filter (K * wlen)
        wlen  % length of control filter
        SecP  % secondary path estimates (M * k * slen)
        slen  % length of secondary path
        Dis   % Disturbance
        Ref   % Reference
        yc    % control filter output
        Nums  % number of secondary sources, K
        Nume  % number of errors, M

    end

    methods
        % initialize (N:duration; PrimaryPath: M x plen)
        function obj = McANC_SRMSE(wLen,SecondaryPath,sLen,Nums,Nume,Dis,Ref)
            N   = length(Ref);
            obj.wlen = wLen;
            obj.Wc   = zeros(Nums,wLen);
            obj.SecP = SecondaryPath;
            obj.slen = sLen;
            obj.yc   = zeros(Nums,N);
            obj.Dis  = Dis;
            obj.Ref  = Ref;
            obj.Nume = Nume;
            obj.Nums = Nums;
        end
        % 1 x 2 x 2
        function [e,obj] = McFxLMS_SRMSE_122(obj,muw)
            N  = length(obj.Ref);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            for i = 1:N
                xc = [obj.Ref(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:));      % y1' received by error 1
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:));      % y2' received by error 2
                
                e(1,i) = obj.Dis(1,i)-y1;                    % error 1
                e(2,i) = obj.Dis(2,i)-y2;                    % error 2
                
                xs = [obj.Ref(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i));  % update control filter 1
                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i));  % update control filter 2
            end
        end

        % 1 x 4 x 4
        function [e,obj] = McFxLMS_SRMSE_144(obj,muw)
            N  = length(obj.Ref);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            for i = 1:N
                xc = [obj.Ref(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update

                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:));
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:));
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:));
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:));
                
                e(1,i) = obj.Dis(1,i)-y1;                    % error 1
                e(2,i) = obj.Dis(2,i)-y2;                    % error 2
                e(3,i) = obj.Dis(3,i)-y3;                    % error 3
                e(4,i) = obj.Dis(4,i)-y4;                    % error 4
                
                xs = [obj.Ref(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(1,3,:) = [sum(xs.*reshape(obj.SecP(3,1,:),[1,obj.slen])) reshape(xf(1,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s31
                xf(1,4,:) = [sum(xs.*reshape(obj.SecP(4,1,:),[1,obj.slen])) reshape(xf(1,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s41

                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22
                xf(2,3,:) = [sum(xs.*reshape(obj.SecP(3,2,:),[1,obj.slen])) reshape(xf(2,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s32
                xf(2,4,:) = [sum(xs.*reshape(obj.SecP(4,2,:),[1,obj.slen])) reshape(xf(2,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s42

                xf(3,1,:) = [sum(xs.*reshape(obj.SecP(1,3,:),[1,obj.slen])) reshape(xf(3,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s13
                xf(3,2,:) = [sum(xs.*reshape(obj.SecP(2,3,:),[1,obj.slen])) reshape(xf(3,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s23
                xf(3,3,:) = [sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen])) reshape(xf(3,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s33
                xf(3,4,:) = [sum(xs.*reshape(obj.SecP(4,3,:),[1,obj.slen])) reshape(xf(3,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s43

                xf(4,1,:) = [sum(xs.*reshape(obj.SecP(1,4,:),[1,obj.slen])) reshape(xf(4,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s14
                xf(4,2,:) = [sum(xs.*reshape(obj.SecP(2,4,:),[1,obj.slen])) reshape(xf(4,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s24
                xf(4,3,:) = [sum(xs.*reshape(obj.SecP(3,4,:),[1,obj.slen])) reshape(xf(4,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s34
                xf(4,4,:) = [sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen])) reshape(xf(4,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s44

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(1,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(1,4,:),[1,obj.wlen])*e(4,i));  % update control filter 1
                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(2,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(2,4,:),[1,obj.wlen])*e(4,i));  % update control filter 2
                obj.Wc(3,:) = obj.Wc(3,:) + muw*(reshape(xf(3,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(3,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(3,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(3,4,:),[1,obj.wlen])*e(4,i));  % update control filter 1
                obj.Wc(4,:) = obj.Wc(4,:) + muw*(reshape(xf(4,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(4,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(4,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(4,4,:),[1,obj.wlen])*e(4,i));  % update control filter 2
            end
        end
            
        % 1 x 6 x 6
        function [e,obj] = McFxLMS_SRMSE_166(obj,muw)
            N  = length(obj.Ref);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            for i = 1:N
                xc = [obj.Ref(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);         % control output y5
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);         % control output y6

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));

                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));

                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));

                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));

                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));

                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                
                e(1,i) = obj.Dis(1,i)-y1;                    % error 1
                e(2,i) = obj.Dis(2,i)-y2;                    % error 2
                e(3,i) = obj.Dis(3,i)-y3;                    % error 3
                e(4,i) = obj.Dis(4,i)-y4;                    % error 4
                e(5,i) = obj.Dis(5,i)-y5;                    % error 5
                e(6,i) = obj.Dis(6,i)-y6;                    % error 6
                
                xs = [obj.Ref(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(1,3,:) = [sum(xs.*reshape(obj.SecP(3,1,:),[1,obj.slen])) reshape(xf(1,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s31
                xf(1,4,:) = [sum(xs.*reshape(obj.SecP(4,1,:),[1,obj.slen])) reshape(xf(1,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s41
                xf(1,5,:) = [sum(xs.*reshape(obj.SecP(5,1,:),[1,obj.slen])) reshape(xf(1,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s51
                xf(1,6,:) = [sum(xs.*reshape(obj.SecP(6,1,:),[1,obj.slen])) reshape(xf(1,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s61

                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22
                xf(2,3,:) = [sum(xs.*reshape(obj.SecP(3,2,:),[1,obj.slen])) reshape(xf(2,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s32
                xf(2,4,:) = [sum(xs.*reshape(obj.SecP(4,2,:),[1,obj.slen])) reshape(xf(2,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s42
                xf(2,5,:) = [sum(xs.*reshape(obj.SecP(5,2,:),[1,obj.slen])) reshape(xf(2,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s52
                xf(2,6,:) = [sum(xs.*reshape(obj.SecP(6,2,:),[1,obj.slen])) reshape(xf(2,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s62

                xf(3,1,:) = [sum(xs.*reshape(obj.SecP(1,3,:),[1,obj.slen])) reshape(xf(3,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s13
                xf(3,2,:) = [sum(xs.*reshape(obj.SecP(2,3,:),[1,obj.slen])) reshape(xf(3,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s23
                xf(3,3,:) = [sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen])) reshape(xf(3,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s33
                xf(3,4,:) = [sum(xs.*reshape(obj.SecP(4,3,:),[1,obj.slen])) reshape(xf(3,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s43
                xf(3,5,:) = [sum(xs.*reshape(obj.SecP(5,3,:),[1,obj.slen])) reshape(xf(3,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s53
                xf(3,6,:) = [sum(xs.*reshape(obj.SecP(6,3,:),[1,obj.slen])) reshape(xf(3,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s63

                xf(4,1,:) = [sum(xs.*reshape(obj.SecP(1,4,:),[1,obj.slen])) reshape(xf(4,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s14
                xf(4,2,:) = [sum(xs.*reshape(obj.SecP(2,4,:),[1,obj.slen])) reshape(xf(4,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s24
                xf(4,3,:) = [sum(xs.*reshape(obj.SecP(3,4,:),[1,obj.slen])) reshape(xf(4,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s34
                xf(4,4,:) = [sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen])) reshape(xf(4,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s44
                xf(4,5,:) = [sum(xs.*reshape(obj.SecP(5,4,:),[1,obj.slen])) reshape(xf(4,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s54
                xf(4,6,:) = [sum(xs.*reshape(obj.SecP(6,4,:),[1,obj.slen])) reshape(xf(4,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s64

                xf(5,1,:) = [sum(xs.*reshape(obj.SecP(1,5,:),[1,obj.slen])) reshape(xf(5,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s15
                xf(5,2,:) = [sum(xs.*reshape(obj.SecP(2,5,:),[1,obj.slen])) reshape(xf(5,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s25
                xf(5,3,:) = [sum(xs.*reshape(obj.SecP(3,5,:),[1,obj.slen])) reshape(xf(5,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s35
                xf(5,4,:) = [sum(xs.*reshape(obj.SecP(4,5,:),[1,obj.slen])) reshape(xf(5,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s45
                xf(5,5,:) = [sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen])) reshape(xf(5,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s55
                xf(5,6,:) = [sum(xs.*reshape(obj.SecP(6,5,:),[1,obj.slen])) reshape(xf(5,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s65

                xf(6,1,:) = [sum(xs.*reshape(obj.SecP(1,6,:),[1,obj.slen])) reshape(xf(6,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s16
                xf(6,2,:) = [sum(xs.*reshape(obj.SecP(2,6,:),[1,obj.slen])) reshape(xf(6,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s26
                xf(6,3,:) = [sum(xs.*reshape(obj.SecP(3,6,:),[1,obj.slen])) reshape(xf(6,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s36
                xf(6,4,:) = [sum(xs.*reshape(obj.SecP(4,6,:),[1,obj.slen])) reshape(xf(6,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s46
                xf(6,5,:) = [sum(xs.*reshape(obj.SecP(5,6,:),[1,obj.slen])) reshape(xf(6,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s56
                xf(6,6,:) = [sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen])) reshape(xf(6,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s66

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(1,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(1,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(1,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(1,6,:),[1,obj.wlen])*e(6,i));  % update control filter 1

                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(2,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(2,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(2,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(2,6,:),[1,obj.wlen])*e(6,i));  % update control filter 2

                obj.Wc(3,:) = obj.Wc(3,:) + muw*(reshape(xf(3,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(3,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(3,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(3,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(3,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(3,6,:),[1,obj.wlen])*e(6,i));  % update control filter 3

                obj.Wc(4,:) = obj.Wc(4,:) + muw*(reshape(xf(4,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(4,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(4,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(4,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(4,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(4,6,:),[1,obj.wlen])*e(6,i));  % update control filter 4

                obj.Wc(5,:) = obj.Wc(5,:) + muw*(reshape(xf(5,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(5,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(5,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(5,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(5,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(5,6,:),[1,obj.wlen])*e(6,i));  % update control filter 5

                obj.Wc(6,:) = obj.Wc(6,:) + muw*(reshape(xf(6,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(6,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(6,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(6,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(6,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(6,6,:),[1,obj.wlen])*e(6,i));  % update control filter 4
            end
        end

        % Arbitrary channel
        function [e,obj] = McFxLMS_SRMSE_ANC(obj,muw)
            N = length(obj.Ref);
            K  = obj.Nums;            % # loudspeakers (secondary sources)
            M  = obj.Nume;            % # error mics
            Lw = obj.wlen;            % control-filter length
            Ls = obj.slen;            % secondary path length

            e        = zeros(M, N);        % error signals
            xc       = zeros(1, Lw);       % x buffer for control filter
            xs       = zeros(1, Ls);       % x buffer for SecP prefiltering
            ys       = zeros(K, Ls);       % per-loudspeaker output ring buffer (for SecP)
            xf       = zeros(K, M, Lw);    % filtered-x buffers (k,m)

            for i = 1:N
                xc = [obj.Ref(i) xc(1:end-1)];
                obj.yc(:,i) = (obj.Wc*xc.').';            % control signal

                ys(:, 2:Ls) = ys(:, 1:Ls-1);
                ys(:, 1)    = obj.yc(:, i);

                % error signal
                for mIdx = 1:M
                    y_m = 0;
                    for kIdx = 1:K
                        s_mk = reshape(obj.SecP(mIdx, kIdx, :), 1, Ls);   % s_{m,k}: 1 x Ls
                        y_m  = y_m + sum(s_mk .* ys(kIdx, :));     % anti noise
                    end
                    e(mIdx, i) = obj.Dis(mIdx, i) - y_m;
                end

                % update filtered reference
                xs = [obj.Ref(i), xs(1:Ls-1)];
                for kIdx = 1:K
                    for mIdx = 1:M
                        s_mk = reshape(obj.SecP(mIdx, kIdx, :), 1, Ls);
                        xfm  = sum(xs .* s_mk);  % current filtered-x sample

                        xf(kIdx, mIdx, 2:Lw) = xf(kIdx, mIdx, 1:Lw-1);
                        xf(kIdx, mIdx, 1)    = xfm;
                    end
                end

                % control filter update
                for kIdx = 1:K
                    grad_k = zeros(1, Lw);        % gradient
                    for mIdx = 1:M
                        grad_k = grad_k + e(mIdx, i) * reshape(xf(kIdx, mIdx, :), 1, Lw);
                    end
                    obj.Wc(kIdx, :) = obj.Wc(kIdx, :) + muw * grad_k;
                end

            end
        end


    end
end


