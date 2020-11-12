function traj = OneImageTraj(traj, par)
    omega1 = par.omega1;
    omega2 = par.omega2;
    dwellTime = par.dwellTime; %us ->s
    spectralDwellTime = pi/omega1;
    time = 0:dwellTime:spectralDwellTime;
    for 