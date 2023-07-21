function gamma=getgamma(nucleus)

%taken from io_loadspec_niimrs.m
switch nucleus
        case '1H'
            gamma = 42.577;
        case '2H'
            gamma = 6.536;
        case '3HE'
            gamma = -32.434;
        case '7LI'
            gamma = 16.546;
        case '13C'
            gamma = 10.708;
        case '19F'
            gamma = 40.052;
        case '23NA'
            gamma = 11.262;
        case '31P'
            gamma = 17.235;
        case '129XE'
            gamma = -11.777;
end
end