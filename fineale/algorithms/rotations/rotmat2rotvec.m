% Convert rotation matrix R to rotation vector V
%
% function V=rotmat2rotvec(R)
%
function V=rotmat2rotvec(R)
rm11    = R(1,1);
rm22    = R(2,2);
rm33    = R(3,3);
rmtrace = rm11 + rm22 + rm33;
maxt    = rmtrace;
imax    = 0;
if (rm11 > maxt)  maxt  = rm11; imax  = 1; end
if (rm22 > maxt)  maxt  = rm22; imax  = 2; end
if (rm33 > maxt)  maxt  = rm33; imax  = 3; end

if        (imax == 0) 
    q0 = 0.5 * sqrt(1 + rmtrace);
    q1 = (R(3,2) - R(2,3)) * 0.25 / q0;
    q2 = (R(1,3) - R(3,1)) * 0.25 / q0;
    q3 = (R(2,1) - R(1,2)) * 0.25 / q0;
elseif (imax == 1) 
    q1 = sqrt(0.5 * rm11 + 0.25 * (1 - rmtrace));
    q0 = (R(3,2) - R(2,3)) * 0.25 / q1;
    q2 = (R(1,2) + R(2,1)) * 0.25 / q1;
    q3 = (R(3,1) + R(1,3)) * 0.25 / q1;
    if (q0 < 0)  q0 = -q0; q1 = -q1; q2 = -q2; q3 = -q3; end
elseif (imax == 2) 
    q2 = sqrt(0.5 * rm22 + 0.25 *(1 - rmtrace));
    q0 = (R(1,3) - R(3,1)) * 0.25 / q2;
    q1 = (R(1,2) + R(2,1)) * 0.25 / q2;
    q3 = (R(3,2) + R(2,3)) * 0.25 / q2;
    if (q0 < 0)  q0 = -q0; q1 = -q1; q2 = -q2; q3 = -q3; end
elseif (imax == 3) 
    q3 = sqrt(0.5 * rm33 + 0.25 * (1 - rmtrace));
    q0 = (R(2,1) - R(1,2)) * 0.25 / q3;
    q2 = (R(3,2) + R(2,3)) * 0.25 / q3;
    q1 = (R(3,1) + R(1,3)) * 0.25 / q3;
    if (q0 < 0) q0 = -q0; q1 = -q1; q2 = -q2; q3 = -q3; end
end
qq = sqrt(q1*q1 + q2*q2 + q3*q3);
if (qq <= 0) 
    V=[0 0 0]';
else 
    if (qq < 0.5)
        rotvn = 2 * asin(qq); % this is more precise for smaller angles */
    else 
        rotvn = 2 * acos(q0); % this is more precise than sin for rotvn ~ PI */
    end
    V = rotvn * [ q1 q2 q3 ]' / qq;
end
