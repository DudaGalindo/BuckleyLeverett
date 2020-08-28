function [index, Swf, dfwf, bf] = calSatFront(Sw, fw, dfw)
    %% calculate advance front saturation using Welge graphic method
    % Sw  saturation vector of displacing wetting phase fluid
    % fw  fractional flow vector of displacing wetting phase fluid
    % dfw  fractional flow derivate vector of displacing wetting phase fluid
    % index  index of Swf in Sw
    % Swf  advance front saturation
    % dfwf  fractional flow derivative at Swf
    % bf  intercept value for tangent line
    index = 0;
    ns = size(Sw,1);
    df = 1/(Sw(ns)-Sw(1));
    dlt_Sw = Sw(ns)-Sw(ns-1);
    
    if nargin == 2
    % calculate derivate
        for i = 3:(ns-2)
            dfw(i) = (fw(i+1)-fw(i-1))/(2*dlt_Sw);
        end
        dfw(2) = (-11*fw(2)+18*fw(3)-9*fw(4)+2*fw(5))/(6*dlt_Sw);
        dfw(1) = abs(2*dfw(2)-dfw(3));
        dfw(ns-1) = -(-11*fw(nt-1)+18*fw(nt-2)-9*fw(nt-3)+2*fw(nt-4))/(6*dlt_Sw);
        dfw(ns) = abs(2*dfw(nt-1)-dfw(nt-2));
    % calculate Swf and its index
        for i = 2:ns-1
            dfds = fw(i)/(Sw(i)-Sw(1));
            if (dfds < dfw(i1)) && (dfds > dfw(i+1)) && (dfds >= df)
                index = i;
                Swf = Sw(i);
                dfwf = dfds;
                bf = fw(i)-dfwf*Sw(i);
                break;
            end
        end
        if index == 0
            error('Can’t find a tangent line through point (Swc,0). Decrease dlt_Sw!');
        end
    elseif nargin == 3
    % calculate Swf and its index
        for i = 2:ns-1
            dfds = fw(i)/(Sw(i)-Sw(1));
            if (dfds < dfw(i-1)) && (dfds > dfw(i+1)) && (dfds >= df)
                index = i;
                Swf = Sw(i);
                dfwf = dfds;
                bf = fw(i)-dfwf*Sw(i);
                break;
            end
        end
        if index == 0
            error('Can’t find a tangent line through point (Swc,0). Decrease dlt_Sw!');
        end
        else
            error('the number of input arguments in function calSatFront is incorrect!');
    end
end
%% end function calSatFront()
