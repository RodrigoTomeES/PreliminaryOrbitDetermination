//$Header$
//------------------------------------------------------------------------------
//                                   Mjday
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of Mjday function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "Mjday.h"

//------------------------------------------------------------------------------
//  double Mjday(double year, double month, double day,
//               double hour, double min, double sec)
//------------------------------------------------------------------------------
/**
* Define constants
*
* @param - double year, month, day, hour, min, sec
* @return - double Mjd (Modified Julian Date)
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double Mjday(double year, double month, double day, double hour, double min, double sec) {
    // nargin - Number of arguments of the function
    // We can't get this function in C, so we always pass six params
    // https://stackoverflow.com/questions/4421681/how-to-count-the-number-of-arguments-passed-to-a-function-that-accepts-a-variabl

    /*
    if (nargin < 4) {
        hour = 0;
        min  = 0;
        sec  = 0;
    }
    */

    double y, m, b, c, a, jd, Mjd;

    y = year;
    m = month;
    b = 0;
    c = 0;

    if (m <= 2) {
       y = y - 1;
       m = m + 12;
    }

    if (y < 0) {
       c = -.75;
    }

    // check for valid calendar date
    if (year < 1582) {
       // null
    } else if (year > 1582) {
       a = fix(y / 100);
       b = 2 - a + floor(a / 4);
    } else if (month < 10) {
       // null
    } else if (month > 10) {
       a = fix(y / 100);
       b = 2 - a + floor(a / 4);
    } else if (day <= 4) {
       // null
    } else if (day > 14) {
       a = fix(y / 100);
       b = 2 - a + floor(a / 4);
    } else {
        printf("\n\n  This is an invalid calendar date!!\n");
        return 0;
    }

    jd = fix(365.25 * y + c) + fix(30.6001 * (m + 1));
    jd = jd + day + b + 1720994.5;
    jd = jd + (hour+min/60+sec/3600)/24;
    Mjd = jd - 2400000.5;

    return Mjd;
}