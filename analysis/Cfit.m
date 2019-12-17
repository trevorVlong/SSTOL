function [a,e2c,NRMSE] = Cfit(Fdata,coeff,odr)
%   produces a curve fit for cl,cx,cm of blown flapped wing with c_ =
%   c_(alpha,Dcj)
%   Detailed explanation goes here
    if coeff >3 || coeff <1
        error('must fit lift, drag, or moment')
    end
    a1 = linspace(-10,10,10);
    a2 = linspace(-10,10,10);
    a3 = linspace(-10,10,10);
%     da1 = 1e-6;
%     da2 = 1e-4;
%     da3 = 1e-2;
%     
    
    
    iter = 0;
%     while iter < 100 
%         iter = iter+1;
%         if iter > 1
%              J = -diag([dfdahat1;dfdahat2*alfai(nalfa);dfdahat3*alfai(nalfa)^3]);
%              delta = diag(inv(J).*Fhat);
%              a1 = a1+da1*delta(1:length(a1));
%              a2 = a2+da2*delta(length(a1)+1:2*length(a1));
%              a3 = a3+da3*delta(2*length(a1)+1:3*length(a1));
%         end
      for ii = 1:length(a1)
        for jj = 1:length(a2)
            for kk = 1:length(a3)
                iter = iter+1;
                    
                for nalfa = 1:length(Fdata(1,coeff,:))
                        Fi = Fdata(:,1,nalfa);
                        alfai = Fdata(:,end,nalfa);
                    for nDcj = 1:length(Fdata(:,1,nalfa))
                        fhat1 = fDcj(a1(ii),Fdata(nDcj,6,nalfa));
                        fhat2 = fDcj(a2(jj),Fdata(nDcj,6,nalfa));
                        fhat3 = fDcj(a3(kk),Fdata(nDcj,6,nalfa));
                        if coeff == 1
                            Fhat             = (fhat1 + fhat2*alfai(nalfa) + fhat3*alfai(nalfa)^3);
                            e2(nalfa,nDcj)   = (Fi(nDcj) - Fhat)^2;
                        else
                            Fhat             = (fhat1 + fhat2*alfai(nalfa) + fhat3*alfai(nalfa)^2);
                            e2(nalfa,nDcj)   = (Fi(nDcj) - Fhat)^2 ;
                        end

                    end
                end
                e2c{iter} = e2;
                NRMSE(iter) = sum(sum(e2));
            end
        end
      end
   
   a = [a1;a2;a3]; 
end