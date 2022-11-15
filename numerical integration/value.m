clc;
clear


krho = 0:0.01:30;
n = 50;

main(krho,n);


%% main function
function main(krho,n)
% Coefficient setting
MU = 4*pi*1e-7;
EPSILON = 8.854e-12;
SIGMA = 0.1;
M_z = 1;
F = 20e3;
% F = 10e9;
OMEGA = 2*pi*F;
SIGMA_complex = SIGMA-1i*OMEGA*EPSILON;
POSITION.X = 0.5;
POSITION.Y = 0;
POSITION.Z = 1;
RHO = sqrt(POSITION.X^2+POSITION.Y^2);
R = sqrt(POSITION.X^2+POSITION.Y^2+POSITION.Z^2);
k2 = 1i*OMEGA*MU*SIGMA_complex;
k = sqrt(k2);
kz = @(x)((real(k2)-x.^2).^2+imag(k2).^2).^(1/4).*exp(1i./2.*atan2(imag(k2),real(k2)-x.^2));

% Objective function
E_phi2int = @(x)-(OMEGA*MU*M_z)/(4*pi)*(x.^2./kz(x)).*besselj(1,x.*RHO).*exp(1i.*abs(POSITION.Z).*kz(x));
H_rho2int = @(x)M_z./(4.*pi).*besselj(1,x.*RHO).*(abs(POSITION.Z)./POSITION.Z).*exp(1i.*abs(POSITION.Z).*kz(x)).*x.^2;
H_z2int = @(x)M_z./(4.*pi).*besselj(0,x.*RHO).*((1i.*x.^3)./kz(x)).*exp(1i.*abs(POSITION.Z).*kz(x));
exp_func = @(x) exp(1i.*abs(POSITION.Z).*kz(x));


% krho = Cheshev_interpolation(linspace(0,ImproperIntegralPoint(E_phi2int,xlist),3000));
% E_phi
E_phi_spline = integralcubicSpline(krho,E_phi2int(krho));
E_phi_GL = integralGL(krho,E_phi2int,n);
E_phi_splineGL = integralGLOfPoly(@(x)exp_func(x)./kz(x),E_phi2int,krho);
E_phi_Glag = integralGlag(@(x)E_phi2int(x).*exp(x),0,1,n);
E_phi_Adapt = integralAdapt(E_phi2int, krho(1), krho(end), krho(end)-krho(1));
% E_phi_Romberg = integralRomberg(E_phi2int, krho(1), krho(end));
E_phi_Value = -M_z/(4*pi)*1i*OMEGA*MU*POSITION.X*exp(1i*k*R)/R^3*(1i*k*R-1);


% krho = Cheshev_interpolation(linspace(0,ImproperIntegralPoint(H_rho2int,xlist),3000));
% H_rho
H_rho_spline = integralcubicSpline(krho,H_rho2int(krho));
H_rho_GL = integralGL(krho,H_rho2int,n);
H_rho_splineGL = integralGLOfPoly(exp_func,H_rho2int,krho);
H_rho_Glag = integralGlag(@(x)H_rho2int(x).*exp(x),0,1,n);
H_rho_Adapt = integralAdapt(H_rho2int, krho(1), krho(end), krho(end)-krho(1));
% H_rho_Romberg = integralRomberg(H_rho2int, krho(1), krho(end));
H_rho_Value = -POSITION.X*POSITION.Z*exp(1i*k*R)/(4*pi)*(k2/R^3+3*1i*k/R^4-3/R^5);

% krho = Cheshev_interpolation(linspace(0,ImproperIntegralPoint(H_z2int,xlist),3000));
% H_z
H_z_spline = integralcubicSpline(krho,H_z2int(krho));
H_z_GL = integralGL(krho,H_z2int,n);
H_z_splineGL = integralGLOfPoly(@(x)1i.*exp_func(x)./kz(x),H_z2int,krho);
H_z_Glag = integralGlag(@(x)H_z2int(x).*exp(x),0,1,n);
H_z_Adapt = integralAdapt(H_z2int, krho(1), krho(end), krho(end)-krho(1));
% H_z_Romberg = integralRomberg(H_z2int, krho(1), krho(end));
H_z_Value = exp(1i*k*R)/(4*pi*R)*(k2+1i*k/R-(k2*POSITION.Z^2+1)/R^2-3*1i*k*POSITION.Z^2/R^3+3*POSITION.Z^2/R^4);

% Result
i = 1;
filename =['./','result','%',datestr(now,29),'%',num2str(i),'.csv'];
while exist(filename,'file')
    i = i+1;
    filename =['./','result','%',datestr(now,29),'%',num2str(i),'.csv'];   
end
fid = fopen(filename,'w');
name  = char('Value','spline','GL','splineGL','Glag','Adapt');
E_phi = [E_phi_Value,E_phi_spline,E_phi_GL,E_phi_splineGL,E_phi_Glag,E_phi_Adapt];
H_rho = [H_rho_Value,H_rho_spline,H_rho_GL,H_rho_splineGL,H_rho_Glag,H_rho_Adapt];
H_z = [H_z_Value,H_z_spline,H_z_GL,H_z_splineGL,H_z_Glag,H_z_Adapt];
for col = 1:length(E_phi)+1
    if col == 1
        fprintf(fid,'%s,%s,%s,%s\n','method','Hx','Hz','Ey');
    else
        fprintf(fid,'%s,%g+%gi,%g+%gi,%g+%gi\n',name(col-1,:),H_rho(col-1),H_rho(col-1)/1i,H_z(col-1),H_z(col-1)/1i,E_phi(col-1),E_phi(col-1)/1i);
    end
end
end

