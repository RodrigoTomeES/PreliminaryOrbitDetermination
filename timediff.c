//$Header$
//------------------------------------------------------------------------------
//                                   timediff
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of timediff function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "timediff.h"

//------------------------------------------------------------------------------
//  void timediff(double UT1_UTC,double TAI_UTC, double timediff_const[])
//------------------------------------------------------------------------------
/**
* Define constants
*
* @param - x
* @return - res
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void timediff(double UT1_UTC,double TAI_UTC,double * UT1_TAI, double * UTC_GPS, double * UT1_GPS, double * TT_UTC, double * GPS_UTC) {
    double TT_TAI  = 32.184;          // TT-TAI time difference [s]

    double GPS_TAI = -19.0;            // GPS-TAI time difference [s]

    double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    double UT1_TAIv = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]

    double UTC_GPSv = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

    double UT1_GPSv = UT1_TAIv-GPS_TAI;  // UT1-GPS time difference [s]

    double TT_UTCv  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

    double GPS_UTCv = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]



    *UT1_TAI = UT1_TAIv;
    *UTC_GPS = UTC_GPSv;
    *UT1_GPS = UT1_GPSv;
    *TT_UTC = TT_UTCv;
    *GPS_UTC = GPS_UTCv;
}