double gmst(double Mjd_UT1){
    double res = 0;

    double Secs = 86400;
    double pi=3.14159265358979;

    double MJD_J2000 = 51544.5;
    double Mjd_0 = floor(Mjd_UT1);
    

    double UT1   = Secs*(Mjd_UT1-Mjd_0); 
    
    double T_0   = (Mjd_0  -MJD_J2000)/36525;
    
    double T     = (Mjd_UT1-MJD_J2000)/36525;
    
    double gmstV=24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1+ (0.093104-6.2e-6*T)*T*T;
    
    res=(gmstV/Secs);
    printf("Valor frac %f \n", Frac(res));
    res=2*pi*Frac(res);

    

    return res;
}