function frequency = getHzFromPPM(shifts, b0)
    gamma = getGamma("overTwoPi", true);
    frequency = shifts * gamma * b0 / 1e6;
end