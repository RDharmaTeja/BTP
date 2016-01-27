% R D Teja rdarmateja@gmail.com
% maximum number of nodes
imax = 97;
jmax = 49;

% import data of node coordinates
% written node coordinates in txt seperatley
X = importdata('xCoord.txt');
Y = importdata('yCoord.txt');

% declaring area variables
xarea = zeros(imax,jmax-1,2);
yarea = zeros(imax-1,jmax,2);

%calculating area normal
%xarea
xarea(1:imax,1:jmax-1,1) = Y(1:imax,2:jmax)-Y(1:imax,1:jmax-1);
xarea(1:imax,1:jmax-1,2) = X(1:imax,1:jmax-1)-X(1:imax,2:jmax);
%yarea
yarea(1:imax-1,1:jmax,1) = Y(1:imax-1,1:jmax) - Y(2:imax,1:jmax);
yarea(1:imax-1,1:jmax,2) = X(2:imax,1:jmax)-X(1:imax-1,1:jmax);

%calculate area magnitude
xareamag(:,:) = sqrt(xarea(:,:,1).*xarea(:,:,1)+xarea(:,:,2).*xarea(:,:,2));
yareamag(:,:) = sqrt(yarea(:,:,1).*yarea(:,:,1)+yarea(:,:,2).*yarea(:,:,2));

% area normals
%xarea
xarean(1:imax,1:jmax-1,1) = xarea(1:imax,1:jmax-1,1)./xareamag(1:imax,1:jmax-1);
xarean(1:imax,1:jmax-1,2) = xarea(1:imax,1:jmax-1,2)./xareamag(1:imax,1:jmax-1);
%yarea
yarean(1:imax-1,1:jmax,1) = yarea(1:imax-1,1:jmax,1)./yareamag(1:imax-1,1:jmax);
yarean(1:imax-1,1:jmax,2) = yarea(1:imax-1,1:jmax,2)./yareamag(1:imax-1,1:jmax);

% declaring vol variables
vol = zeros(imax-1,jmax-1);
% calculating vol
for j = 1:jmax-1
    for i = 1:imax-1
        vol(i,j)= (0.5*abs((X(i+1,j+1)-X(i,j))*(Y(i,j+1)-Y(i,j))-(Y(i+1,j+1)-Y(i,j))*(X(i,j+1)-X(i,j)))+0.5*abs((X(i+1,j+1)-X(i,j))*(Y(i+1,j)-Y(i,j))-(Y(i+1,j+1)-Y(i,j))*(X(i+1,j)-X(i,j))));
    end
end


% state variables
rair = 1.4;%gamma
Rair = 287;%gas const
%Inflow conditions
pInfi = 1013.25; %#pressure in pascals
tInfi = 300; % temp
denInfi = pInfi/(Rair*tInfi);% density
mInfi = 0.4; % mach
uInfi = mInfi*sqrt(rair*pInfi/denInfi);%# u vel in m/s
vInfi = 0 ;%# v vel in m/s
toler = 1e-4;


% Variables q for primitive p,u,v,density decleration
% Everywhere we maintain q and U to make life easy
q = zeros(imax+1,jmax+1,4);
% U: charectaristic variables
U = zeros(imax+1,jmax+1,4);
%switch variables
[mp1,mm1,ap1,am1,bp1,bm1,dp1,dm1,cp1,cm1,vnp1,vnm1,rm1,rp1] = deal(zeros(imax,jmax-1)); %for fFlux
[mp2,mm2,ap2,am2,bp2,bm2,dp2,dm2,cp2,cm2,vnp2,vnm2,rm2,rp2] = deal(zeros(imax-1,jmax)); % for gFlux
%
a = zeros(imax+1,jmax+1); %sound speeds
%residual decleration
res1 = zeros(imax+1,jmax+1,4);
%# Initializing to inflow flow conditions
% Primitive variables
q(:,:,1) = pInfi;
q(:,:,2) = uInfi;
q(:,:,3) = vInfi;
q(:,:,4) = denInfi;
% Primitive to Characteristic variables
U(:,:,1) = q(:,:,4);
U(:,:,2) = q(:,:,4).*q(:,:,2);
U(:,:,3) = q(:,:,4).*q(:,:,3);
U(:,:,4) = (q(:,:,1)/(rair-1))+0.5*(q(:,:,4).*(q(:,:,2).^2+q(:,:,3).^2));


% Main iteration
% useful variables
%
tGlobal = 0 ;%# Global time step
resNorm = 0 ;%# First iteration residual norm
resNormI1 = 0;% # ith iteration res norm
cfl = 0.3 ;%# cfl number for calculating time steps
denUinfi = denInfi*uInfi;
conv(iterCount)=0;
%flux variables decleration
fFlux1 = zeros(imax,jmax-1,4);
gFlux1 = zeros(imax-1,jmax,4);

% ***************************************Boundary Conditions Start********************************%
% Inflow Boundary conditions
% Primitive variables
%q(1,:,1) = pInfi; %#supersonic
q(1,:,1) = q(2,:,1);%subsonic
q(1,:,2) = uInfi;
q(1,:,3) = vInfi;
q(1,:,4) = denInfi;
%  Conservative variables
U(1,:,1) = q(1,:,4);
U(1,:,2) = q(1,:,4).*q(1,:,2);
U(1,:,3) = q(1,:,4).*q(1,:,3);
U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
% Outflow bounday conditions
% Primitive variables
%q(imax+1,:,1) = q(imax,:,1); %# supersonic
q(imax+1,:,1) = pInfi;%subsonic
q(imax+1,:,2) = q(imax,:,2);
q(imax+1,:,3) = q(imax,:,3);
q(imax+1,:,4) = q(imax,:,4);
%  Conservative variables
U(imax+1,:,1) = q(imax+1,:,4);
U(imax+1,:,2) = q(imax+1,:,4).*q(imax+1,:,2);
U(imax+1,:,3) = q(imax+1,:,4).*q(imax+1,:,3);
U(imax+1,:,4) = (q(imax+1,:,1)/(rair-1))+0.5*(q(imax+1,:,4).*(q(imax+1,:,2).^2+q(imax+1,:,3).^2));
% Wall bounday conditions
% 'p',densi at G is equal to the 'p',densi at I
% Bottom wall
q(:,1,1) = q(:,2,1);
q(:,1,4) = q(:,2,4);
q(2:imax,1,2) = q(2:imax,2,2) - (2*(q(2:imax,2,2).*yarean(1:imax-1,1,1)+q(2:imax,2,3).*yarean(1:imax-1,1,2))).*yarean(1:imax-1,1,1);
q(2:imax,1,3) = q(2:imax,2,3) - (2*(q(2:imax,2,2).*yarean(1:imax-1,1,1)+q(2:imax,2,3).*yarean(1:imax-1,1,2))).*yarean(1:imax-1,1,2);
U(:,1,1) = q(:,1,4);
U(:,1,2) = q(:,1,4).*q(:,1,2);
U(:,1,3) = q(:,1,4).*q(:,1,3);
U(:,1,4) = (q(:,1,1)/(rair-1))+0.5*(q(:,1,4).*(q(:,1,2).^2+q(:,1,3).^2));
% Top wall
q(:,jmax+1,1) = q(:, jmax,1);
q(:,jmax+1,4) = q(:, jmax,4);
q(2:imax,jmax+1,2) = q(2:imax,jmax,2) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,1);
q(2:imax,jmax+1,3) = q(2:imax,jmax,3) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,2);
U(:,jmax+1,1) = q(:,jmax+1,4);
U(:,jmax+1,2) = q(:,jmax+1,4).*q(:,jmax+1,2);
U(:,jmax+1,3) = q(:,jmax+1,4).*q(:,jmax+1,3);
U(:,jmax+1,4) = (q(:,jmax+1,1)/(rair-1))+0.5*(q(:,jmax+1,4).*(q(:,jmax+1,2).^2+q(:,jmax+1,3).^2));
% ***********************************Boundary Conditions End ********************************%

for iters = 1:iterCount
    %*************************************** Sound speed and Mach number start ***************************** %
    % Sound speed
    a(1:imax+1,1:jmax+1) = (sqrt(rair*q(1:imax+1,1:jmax+1,1)./q(1:imax+1,1:jmax+1,4))); % sound speeds
    % Mach number
    mp1(1:imax,1:jmax-1) = (q(1:imax,2:jmax,2).*xarean(1:imax,1:jmax-1,1)+q(1:imax,2:jmax,3).*xarean(1:imax,1:jmax-1,2))./(0.5*(a(1:imax,2:jmax)+a(2:imax+1,2:jmax)));
    mm1(1:imax,1:jmax-1) = (q(2:imax+1,2:jmax,2).*xarean(1:imax,1:jmax-1,1)+q(2:imax+1,2:jmax,3).*xarean(1:imax,1:jmax-1,2))./(0.5*(a(1:imax,2:jmax)+a(2:imax+1,2:jmax)));
    mp2(1:imax-1,1:jmax) = (q(2:imax,1:jmax,2).*yarean(1:imax-1,1:jmax,1)+q(2:imax,1:jmax,3).*yarean(1:imax-1,1:jmax,2))./(0.5*(a(2:imax,1:jmax)+a(2:imax,2:jmax+1)));
    mm2(1:imax-1,1:jmax) = (q(2:imax,2:jmax+1,2).*yarean(1:imax-1,1:jmax,1)+q(2:imax,2:jmax+1,3).*yarean(1:imax-1,1:jmax,2))./(0.5*(a(2:imax,1:jmax)+a(2:imax,2:jmax+1)));
    %*************************************** Sound speed and Mach number end ***************************** %
    % ************************* Switches  Start ***********************************************%
    % betas
    bp1 = min(0,sign(abs(mp1)-1));
    bm1 = min(0,sign(abs(mm1)-1));
    bp2 = min(0,sign(abs(mp2)-1));
    bm2 = min(0,sign(abs(mm2)-1));
    % Alphas
    ap1 =0.5*(1+sign(mp1));
    am1 =0.5*(1-sign(mm1));
    ap2 =0.5*(1+sign(mp2));
    am2 =0.5*(1-sign(mm2));
    % cvl s
    cp1 =ap1.*(1+bp1).*mp1 - 0.25*bp1.*((1+mp1).^2);
    cm1 =am1.*(1+bm1).*mm1 + 0.25*bm1.*((1-mm1).^2);
    cp2 =ap2.*(1+bp2).*mp2 - 0.25*bp2.*((1+mp2).^2);
    cm2 =am2.*(1+bm2).*mm2 + 0.25*bm2.*((1-mm2).^2);
    % Ds
    dp1 =ap1.*(1+bp1)-0.25*bp1.*((1+mp1).^2).*(2-mp1);
    dm1 =am1.*(1+bm1)-0.25*bm1.*((1-mm1).^2).*(2+mm1);
    dp2 =ap2.*(1+bp2)-0.25*bp2.*((1+mp2).^2).*(2-mp2);
    dm2 =am2.*(1+bm2)-0.25*bm2.*((1-mm2).^2).*(2+mm2);
    % ************************* Switches  Ends ***********************************************%
    % ************************* Flux Construction Starts ***********************************************%
    %fFlux
    fFlux1(1:imax,1:jmax-1,1) = (q(1:imax,2:jmax,4).*a(1:imax,2:jmax).*cp1 + q(2:imax+1,2:jmax,4).*a(2:imax+1,2:jmax).*cm1).*xareamag;
    fFlux1(1:imax,1:jmax-1,2) = (U(1:imax,2:jmax,2).*a(1:imax,2:jmax).*cp1 + U(2:imax+1,2:jmax,2).*a(2:imax+1,2:jmax).*cm1 + (dp1.*q(1:imax,2:jmax,1)+dm1.*q(2:imax+1,2:jmax,1)).*xarean(:,:,1)).*xareamag;
    fFlux1(1:imax,1:jmax-1,3) = (U(1:imax,2:jmax,3).*a(1:imax,2:jmax).*cp1 + U(2:imax+1,2:jmax,3).*a(2:imax+1,2:jmax).*cm1 + (dp1.*q(1:imax,2:jmax,1)+dm1.*q(2:imax+1,2:jmax,1)).*xarean(:,:,2)).*xareamag;
    fFlux1(1:imax,1:jmax-1,4) = (q(1:imax,2:jmax,4).*a(1:imax,2:jmax).*cp1.*(3.5*q(1:imax,2:jmax,1)./q(1:imax,2:jmax,4)+0.5*(q(1:imax,2:jmax,2).^2+q(1:imax,2:jmax,3).^2))...
        + q(2:imax+1,2:jmax,4).*a(2:imax+1,2:jmax).*cm1.*(3.5*q(2:imax+1,2:jmax,1)./q(2:imax+1,2:jmax,4)+0.5*(q(2:imax+1,2:jmax,2).^2+q(2:imax+1,2:jmax,3).^2))).*xareamag;

    % gFlux
    gFlux1(1:imax-1,1:jmax,1) = (q(2:imax,1:jmax,4).*a(2:imax,1:jmax).*cp2 + q(2:imax,2:jmax+1,4).*a(2:imax,2:jmax+1).*cm2).*yareamag;
    gFlux1(1:imax-1,1:jmax,2) = (U(2:imax,1:jmax,2).*a(2:imax,1:jmax).*cp2 + U(2:imax,2:jmax+1,2).*a(2:imax,2:jmax+1).*cm2 + (dp2.*q(2:imax,1:jmax,1)+dm2.*q(2:imax,2:jmax+1,1)).*yarean(:,:,1)).*yareamag;
    gFlux1(1:imax-1,1:jmax,3) = (U(2:imax,1:jmax,3).*a(2:imax,1:jmax).*cp2 + U(2:imax,2:jmax+1,3).*a(2:imax,2:jmax+1).*cm2 + (dp2.*q(2:imax,1:jmax,1)+dm2.*q(2:imax,2:jmax+1,1)).*yarean(:,:,2)).*yareamag;
    gFlux1(1:imax-1,1:jmax,4) = (q(2:imax,1:jmax,4).*a(2:imax,1:jmax).*cp2.*(3.5*q(2:imax,1:jmax,1)./q(2:imax,1:jmax,4)+0.5*(q(2:imax,1:jmax,2).^2+q(2:imax,1:jmax,3).^2))...
        + q(2:imax,2:jmax+1,4).*a(2:imax,2:jmax+1).*cm2.*(3.5*q(2:imax,2:jmax+1,1)./q(2:imax,2:jmax+1,4)+0.5*(q(2:imax,2:jmax+1,2).^2+q(2:imax,2:jmax+1,3).^2))).*yareamag;
    % ************************* Flux Construction Ends ***********************************************%
    % ************************** No penetration condition Starts *************************************%
    % # Bottom
    gFlux1(1:imax-1,1,1)=0;
    gFlux1(1:imax-1,1,4)=0;
    gFlux1(1:imax-1,1,2)= (dp2(1:imax-1,1).*q(2:imax,2,1)+dm2(1:imax-1,1).*q(2:imax,2,1)).*yarea(1:imax-1,1,1);
    gFlux1(1:imax-1,1,3)= (dp2(1:imax-1,1).*q(2:imax,2,1)+dm2(1:imax-1,1).*q(2:imax,2,1)).*yarea(1:imax-1,1,2);
    gFlux1(1:imax-1,jmax,1)=0;
    gFlux1(1:imax-1,jmax,4)=0;
    gFlux1(1:imax-1,jmax,2)= (dp2(1:imax-1,jmax).*q(2:imax,jmax+1,1)+dm2(1:imax-1,jmax).*q(2:imax,jmax+1,1)).*yarea(1:imax-1,jmax,1);
    gFlux1(1:imax-1,jmax,3)= (dp2(1:imax-1,jmax).*q(2:imax,jmax+1,1)+dm2(1:imax-1,jmax).*q(2:imax,jmax+1,1)).*yarea(1:imax-1,jmax,2);
    % **************************No penetration condition Ends *************************************%
    %***************************Residual and resnorm starts*****************************************%
    %calculating residual
    res1(2:imax,2:jmax,:) = fFlux1(2:imax,1:jmax-1,:) - fFlux1(1:imax-1,1:jmax-1,:)+ gFlux1(1:imax-1,2:jmax,:) - gFlux1(1:imax-1,1:jmax-1,:);
    %calculating residual norm
    res11 = (res1(2:imax,2:jmax,1)/denInfi).^2;
    res12 = (res1(2:imax,2:jmax,2)/denUinfi).^2;
    res13 = (res1(2:imax,2:jmax,3)/denUinfi).^2;
    res14 = (res1(2:imax,2:jmax,4)/(denUinfi*uInfi)).^2;
    resNormI1 = sqrt(sum(res11(:))+sum(res12(:))+sum(res13(:))+sum(res14(:)));
    %***************************Residual and resnorm end*****************************************%
    % **************************Time step starts*************************************************%
    %essentials for time step calculations
    vnp1(1:imax,1:jmax-1) = q(1:imax,2:jmax,2).*xarean(1:imax,1:jmax-1,1)+ q(1:imax,2:jmax,3).*xarean(1:imax,1:jmax-1,2);
    vnm1(1:imax,1:jmax-1) = q(2:imax+1,2:jmax,2).*xarean(1:imax,1:jmax-1,1)+ q(2:imax+1,2:jmax,3).*xarean(1:imax,1:jmax-1,2);
    rp1 = abs(vnp1)+sqrt(rair*q(1:imax,2:jmax,1)./q(1:imax,2:jmax,4));
    rm1 = abs(vnm1)+sqrt(rair*q(2:imax+1,2:jmax,1)./q(2:imax+1,2:jmax,4));
    vnp2(1:imax-1,1:jmax) = q(2:imax,1:jmax,2).*yarean(1:imax-1,1:jmax,1)+ q(2:imax,1:jmax,3).*yarean(1:imax-1,1:jmax,2);
    rp2= abs(vnp2)+sqrt(rair*q(2:imax,1:jmax,1)./q(2:imax,1:jmax,4));
    vnm2(1:imax-1,1:jmax) = q(2:imax,2:jmax+1,2).*yarean(1:imax-1,1:jmax,1)+ q(2:imax,2:jmax+1,3).*yarean(1:imax-1,1:jmax,2);
    rm2 = abs(vnm2)+sqrt(rair*q(2:imax,2:jmax+1,1)./q(2:imax,2:jmax+1,4));
    % local time stepping
    % directly calculate time/volume
    tlv(2:imax,2:jmax) = (2*cfl)./(xareamag(1:imax-1,1:jmax-1).*rm1(1:imax-1,1:jmax-1) + xareamag(2:imax,1:jmax-1).*rp1(2:imax,1:jmax-1)+yareamag(1:imax-1,1:jmax-1).*rm2(1:imax-1,1:jmax-1)+yareamag(1:imax-1,2:jmax).*rp2(1:imax-1,2:jmax));
    %*************************** Time step ends *****************************************%
    %*************************** update sol start *****************************************%
    % Conservative update solution
    U(2:imax,2:jmax,1) = U(2:imax,2:jmax,1) - tlv(2:imax,2:jmax).*res1(2:imax,2:jmax,1);
    U(2:imax,2:jmax,2) = U(2:imax,2:jmax,2) - tlv(2:imax,2:jmax).*res1(2:imax,2:jmax,2);
    U(2:imax,2:jmax,3) = U(2:imax,2:jmax,3) - tlv(2:imax,2:jmax).*res1(2:imax,2:jmax,3);
    U(2:imax,2:jmax,4) = U(2:imax,2:jmax,4) - tlv(2:imax,2:jmax).*res1(2:imax,2:jmax,4);
    % Primitive update solution
    q(:,:,4) = U(:,:,1);
    q(:,:,2) = U(:,:,2)./U(:,:,1);
    q(:,:,3) = U(:,:,3)./U(:,:,1);
    q(:,:,1) = (rair-1)*(U(:,:,4)-0.5*(U(:,:,2).^2+U(:,:,3).^2)./U(:,:,1));
    %*************************** update sol ends *****************************************%
    % init res norm
    if iters==1
        resNorm = resNormI1;
    elseif (resNormI1/resNorm)< toler
        break
    end
    % ***************************************Boundary Conditions Start********************************%
    % Inflow Boundary conditions
    % Primitive variables
    %q(1,:,1) = pInfi; %#supersonic
    q(1,:,1) = q(2,:,1);%subsonic
    q(1,:,2) = uInfi;
    q(1,:,3) = vInfi;
    q(1,:,4) = denInfi;
    %  Conservative variables
    U(1,:,1) = q(1,:,4);
    U(1,:,2) = q(1,:,4).*q(1,:,2);
    U(1,:,3) = q(1,:,4).*q(1,:,3);
    U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
    % Outflow bounday conditions
    % Primitive variables
    %q(imax+1,:,1) = q(imax,:,1); %# supersonic
    q(imax+1,:,1) = pInfi;%subsonic
    q(imax+1,:,2) = q(imax,:,2);
    q(imax+1,:,3) = q(imax,:,3);
    q(imax+1,:,4) = q(imax,:,4);
    %  Conservative variables
    U(imax+1,:,1) = q(imax+1,:,4);
    U(imax+1,:,2) = q(imax+1,:,4).*q(imax+1,:,2);
    U(imax+1,:,3) = q(imax+1,:,4).*q(imax+1,:,3);
    U(imax+1,:,4) = (q(imax+1,:,1)/(rair-1))+0.5*(q(imax+1,:,4).*(q(imax+1,:,2).^2+q(imax+1,:,3).^2));
    % Wall bounday conditions
    % 'p',densi at G is equal to the 'p',densi at I
    % Bottom wall
    q(:,1,1) = q(:,2,1);
    q(:,1,4) = q(:,2,4);
    q(2:imax,1,2) = q(2:imax,2,2) - (2*(q(2:imax,2,2).*yarean(1:imax-1,1,1)+q(2:imax,2,3).*yarean(1:imax-1,1,2))).*yarean(1:imax-1,1,1);
    q(2:imax,1,3) = q(2:imax,2,3) - (2*(q(2:imax,2,2).*yarean(1:imax-1,1,1)+q(2:imax,2,3).*yarean(1:imax-1,1,2))).*yarean(1:imax-1,1,2);
    U(:,1,1) = q(:,1,4);
    U(:,1,2) = q(:,1,4).*q(:,1,2);
    U(:,1,3) = q(:,1,4).*q(:,1,3);
    U(:,1,4) = (q(:,1,1)/(rair-1))+0.5*(q(:,1,4).*(q(:,1,2).^2+q(:,1,3).^2));
    % Top wall
    q(:,jmax+1,1) = q(:, jmax,1);
    q(:,jmax+1,4) = q(:, jmax,4);
    q(2:imax,jmax+1,2) = q(2:imax,jmax,2) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,1);
    q(2:imax,jmax+1,3) = q(2:imax,jmax,3) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,2);
    U(:,jmax+1,1) = q(:,jmax+1,4);
    U(:,jmax+1,2) = q(:,jmax+1,4).*q(:,jmax+1,2);
    U(:,jmax+1,3) = q(:,jmax+1,4).*q(:,jmax+1,3);
    U(:,jmax+1,4) = (q(:,jmax+1,1)/(rair-1))+0.5*(q(:,jmax+1,4).*(q(:,jmax+1,2).^2+q(:,jmax+1,3).^2));
    % ***********************************Boundary Conditions End ********************************%
    conv(iters) = resNormI1/resNorm;
    %end of main iteration loop
end
% flow tangency condition
q(2:imax,1,2) = q(2:imax,2,2) - (2*(q(2:imax,2,2).*yarean(1:imax-1,1,1)+q(2:imax,2,3).*yarean(1:imax-1,1,2))).*yarean(1:imax-1,1,1);
q(2:imax,1,3) = q(2:imax,2,3) - (2*(q(2:imax,2,2).*yarean(1:imax-1,1,1)+q(2:imax,2,3).*yarean(1:imax-1,1,2))).*yarean(1:imax-1,1,2);
q(2:imax,jmax+1,2) = q(2:imax,jmax,2) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,1);
q(2:imax,jmax+1,3) = q(2:imax,jmax,3) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,2);
% calculating at nodes for post processing
qNodes(1:imax,1:jmax,:) = 0.25*(q(2:imax+1,2:jmax+1,:)+q(1:imax,2:jmax+1,:)+q(1:imax,1:jmax,:)+q(2:imax+1,1:jmax,:));
