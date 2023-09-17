function [lamda,ex,P0]=ekf(F,G,Q,R,P0,u,z,ex,lamda)
Xn=F*ex+G*u;

Zn=[atan( Xn(2)/sqrt(Xn(1)^2+Xn(3)^2) ),atan(-1*Xn(3)/Xn(1))]';
P=lamda*F*P0*F'+Q;
dh1_dx=-1*Xn(1)*Xn(2)/(Xn(1)^2+Xn(2)^2+Xn(3)^2)/sqrt(Xn(1)^2+Xn(3)^2);
dh1_dy=sqrt(Xn(1)^2+Xn(3)^2)/(Xn(1)^2+Xn(2)^2+Xn(3)^2);
dh1_dz=-1*Xn(2)*Xn(3)/(Xn(1)^2+Xn(2)^2+Xn(3)^2)/sqrt(Xn(1)^2+Xn(3)^2);
dh2_dx=Xn(3)/(Xn(1)^2+Xn(3)^2);
dh2_dy=0;
dh2_dz=-1*Xn(1)/(Xn(1)^2+Xn(3)^2);
H=[dh1_dx,dh1_dy,dh1_dz,0,0,0,0,0,0;dh2_dx,dh2_dy,dh2_dz,0,0,0,0,0,0];
K=P*H'/(H*P*H'+R);
% sigma=z-Zn;
% Nwsx=H*P*H'-H*F*P0*F'*H'-R;
% Mwsx=H*Q*H';
% Nwsx=H*P*H'-H*Q*H';
% beta=(sigma'*sigma-trace(H*Q*H'+R))/trace(H*F*P0*F'*H');
% Mwsx=beta*H*(P-Q)*H';
% C=trace(Nwsx)/trace(Mwsx);
%     if C>1
%     lamda=C;
%     else
%     lamda=1;
%     end
ex=Xn+K*(z-Zn);
%ex=Xn;
P0=(eye(9)-K*H)*P;
% a=ex'*ex/(trace(P0)+ex'*ex);
% % a=1-trace(P0)/(2*ex'*ex+1);
% ex=a*ex;  
end