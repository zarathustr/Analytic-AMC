% A Simplified Analytic Attitude Determination Algorithm Using Accelerometer and Magnetometer
% Fast Accelerometer-Magnetometer Combination (FAMC) algorithm by Zhuohua Liu and Jin Wu
%
% author: Zhuohua Liu, Jin Wu
% e-mail: liuzhuohua@bupt.edu.cn; jin_wu_uestc@hotmail.com; 

function q = FAMC(Ab, Mb)

    ax = Ab(1);       ay = Ab(2);       az = Ab(3);
    mx = Mb(1);       my = Mb(2);       mz = Mb(3);

    mD = dot(Ab, Mb);
    mN = sqrt(1 - mD^2);

    B11 = (mN * mx) / 2; B13 = ax / 2 + (mD * mx) / 2;
    B21 = (mN * my) / 2; B23 = ay / 2 + (mD * my) / 2;
    B31 = (mN * mz) / 2; B33 = az / 2 + (mD * mz) / 2;

    tau = (B13 + B31);

    p1  = B33 - B11 + 1;
    A11 = -1 / p1;
    A12 = B21 / p1;
    A13 = tau / p1;

    p2  = (-B21^2 / p1 + B11 + B33 + 1);
    A21 =  -B21 / (p1 * p2);
    A22 =  -1 / p2;
    A23 = ( B23 + B21 * tau / p1) / p2;

    p3  = p1 - 2 + tau^2 / p1 + A23^2 * p2;

    A31 = (tau / p1 + B21 * A23 / p1) / p3;
    A32 = A23 / p3;
    A33 = 1 / p3;

    a = B23 * (A11 + A12 * (A21 + A23 * A31) + A13*A31) ...
        - (B13 - B31) * (A21 + A23 * A31) - A31 * B21;
    b = B23 * (A12 * (A22 + A23 * A32) + A13 * A32) ...
        - (B13 - B31) * (A22 + A23 * A32) - A32 * B21;
    
    c = B23 * (A13 * A33 + A12 * A23 * A33) - A33 * B21 - A23 * A33 * (B13 - B31);

    q = [-1; a; b; c];
    q = q ./ norm(q);
end
