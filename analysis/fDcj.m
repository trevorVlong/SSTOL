function [F,dFda] = fDcj(a,Dcj)
        sz = size(a);
        if sz(2)>1
            error('input a must be a column vecotr');
        end
%linear function for effect of Dcj on lift
%   Function values and derivatives in terms of a() coefficients 
        if length(a) == 1
            F  = a*Dcj;
            dFda = Dcj;
            %not a constant, but a special case of linear with no offset
        elseif length(a) == 2
            F  = a(1) + a(2)*Dcj;
            dFda = [1;Dcj];
            
        elseif length(a) == 3
            F  = a(1) + a(2)*Dcj + a(3)*Dcj^2;
            dFda = [1;Dcj;Dcj^2];
            
        elseif length(a) == 4
            F  = a(1) + a2*Dcj + a3*Dcj^2 + a4*Dcj^3;
            dFda = [1;Dcj;Dcj^2;Dcj^3];
        else
            error('dimension of a not supported')
            
        end

end

