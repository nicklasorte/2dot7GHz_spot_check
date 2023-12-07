function g=itu_ant_gain(app,deploy_type,phi,theta)

% Computes antenna gain pattern for 4G base stations following ITU-R F.1336
% (2/2014) with parameters from ITU-R M.2292. 
% Inputs:
%           deploy_type     type of base station per F.1336:
%                           'macro rural'
%                           'macro suburban'
%                           'macro urban'
%                           'micro urban'
%                           'indoor urban'
%           phi             azimuth angle (deg)
%           theta           elevation angle (deg)

omni = 0;

switch deploy_type          % Base station deployment type
    case 'macro rural'      % per M.2292, macro cells are sectorized
        g0 = 18;            % Peak gain (dBi)
        beta = 3;           % Downtilt (deg)
        phi3 = 65;          % 3-dB azimuth beamwidth
    case 'macro suburban'
        g0 = 16;
        beta = 6;
        phi3 = 65;          
    case 'macro urban'
        g0 = 16;
        beta = 10;
        phi3 = 65;          
    case 'micro urban'
        g0 = 5;             % Note that, per M.2292, last two types are omni
        omni = 1;
    case 'indoor urban'
        g0 = 0;
        omni = 1;
    case 'ligado'
        g0 = 17;
        beta = 3;
        phi3 = 72;
    otherwise
        error('Not a valid deployment type');
        return;
end

if omni            % Omni antenna - see Recommendation 2 in F.1336
    theta3 = 107.6*10.^(-g0/10);    
    k = 0.7;        % Use 0.7 for "typical" and 0 for "improved side-lobe"
    theta4 = theta3*sqrt(1-1/1.2*log10(k+1));
    
    g = g0 - (12*(theta/theta3).^2).*(abs(theta) < theta4) + (-12+10*log10(k+1)).*((abs(theta) < theta3) & (abs(theta)>= theta4)) + (-12+10*log10((abs(theta)/theta3).^(-1.5) + k)).*((abs(theta) <= 90) & (abs(theta) >= theta3));
    i = find(theta == 0);
    g(i) = g0;
else                % Sectorized antenna - see Recommendation 3.1 in F.1336
    theta3 = 31000./phi3.*10.^(-0.1*g0);
    theta_e = (90*(theta+beta)./(90+beta)).*(theta >= -beta) + (90*(theta+beta)./(90-beta)).*(theta < -beta);
    ka = 0.7;kp = 0.7;kh = 0.7;kv = 0.3; % Other parameters from M.2292

    xv = abs(theta_e)/theta3;
    g180 = -12 + 10*log10(1+8*kp)-15*log10(180/theta3);
    xh = 180/phi3;
    gh180 = -12*xh.*xh.*(xh<= 0.5) + (-12*xh.^(2-kh)-3*(1-.5.^(-kh))).*(xh> 0.5);
    gh180 = (gh180 >= g180).*gh180 + (gh180 < g180).*g180;
    xh = 0;
    gh0 = -12*xh.*xh.*(xh<= 0.5) + (-12*xh.^(2-kh)-3*(1-.5.^(-kh))).*(xh> 0.5);
    gh0 = (gh0 >= g180).*gh0 + (gh0 < g180).*g180;
    xh = abs(phi)/phi3;
    gh = -12*xh.*xh.*(xh<= 0.5) + (-12*xh.^(2-kh)-3*(1-.5.^(-kh))).*(xh> 0.5);
    gh = (gh >= g180).*gh + (gh < g180).*g180;

    R = (gh - gh180)./(gh0 - gh180);
    C = 10*log10(((180/theta3)^1.5*(4^(-1.5)+kv))./(1+8*kp))./log10(22.5/theta3);
    lambda_kv = 12-C*log10(4)-10*log10(4^(-1.5)+kv);
    xk = sqrt(1-0.36*kv);
    gv = (xv < xk).*(-12*xv.*xv) + (-12+10*log10(xv.^(-1.5) + kv)).*((xv < 4) & (xv>=xk)) + (-lambda_kv-C*log10(xv)).*((xv < 90/theta3) & (xv >= 4)) + g180.*(xv == 90/theta3);
    i = find(xv == 0);
    gv(i) = 0;

    g = g0 + gh + R.*gv;
end

