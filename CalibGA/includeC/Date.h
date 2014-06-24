#ifndef GREG_CALEND
#define GREG_CALEND

#include "Constants.h"

//////////////////////////////////////////////////////
////   Date: Gregorian calendar
//////////////////////////////////////////////////////

const int hours_in_day    = 24;
const int minutes_in_day  = hours_in_day * 60;
const int minutes_to_noon = (hours_in_day / 2) * 60;

int MOD(int x, int y) { 
    return x - (x/y)*y; 
    // return x & y;
}

int date_of_gregorian( const Date& date ) {
    int ym;
    if(date.month == 1 || date.month == 2) {
        ym = ( 1461 * ( date.year + 4800 - 1 ) ) / 4 +
                ( 367 * ( date.month + 10 ) ) / 12 -
                ( 3 * ( ( date.year + 4900 - 1 ) / 100 ) ) / 4;
    } else {
        ym = ( 1461 * ( date.year + 4800 ) ) / 4 +
                ( 367 * ( date.month - 2 ) ) / 12 -
                ( 3 * ( ( date.year + 4900 ) / 100 ) ) / 4;
    }

    int tmp = ym + date.day - 32075 - 2444238;
    
    return tmp * minutes_in_day + date.hour * 60 + date.min;
}


Date gregorian_of_date ( const int& minutes_since_epoch ) {
    int jul, l, n, i, j, d, m, y, daytime;

    jul = minutes_since_epoch / minutes_in_day;
    l   = jul + 68569 + 2444238;
    n   = ( 4 * l ) / 146097;
    l   = l - ( 146097 * n + 3 ) / 4;
    i   = ( 4000 * ( l + 1 ) ) / 1461001;
    
    l   = l - ( 1461 * i ) / 4 + 31;
    j   = ( 80 * l ) / 2447;
    d   = l - ( 2447 * j ) / 80;
    l   = j / 11;

    m   = j + 2 - ( 12 * l );
    y   = 100 * ( n - 49 ) + i + l;

    daytime = MOD( minutes_since_epoch, minutes_in_day );

    if ( daytime == minutes_to_noon )  return Date(y,m,d,12,0);
    else                               return Date(y,m,d, daytime / 60, MOD(daytime, 60));
}

bool check_date(const int& year, const int& month, const int& day) {
    bool tmp1, tmp2, tmp3;

    tmp1 = ( 1 <= day && 1 <= month && month <= 12 && 1980 <= year && year <= 2299 );
    tmp2 = ( day <= 28 );

    if        ( month == 2 ) {
        tmp3 = ( day == 29 && MOD(year, 4) == 0 && ( year == 2000 || (not (MOD(year, 100) == 0)) ) );
    } else if ( month == 4 || month == 6 || month == 9 || month == 11 ) {
        tmp3 = ( day <= 30 );
    } else {
        tmp3 = ( day <= 31 );
    }

    return tmp1 && (tmp2 || tmp3);
}

REAL days_between(const REAL& t1, const REAL& t2) {
    REAL diff = static_cast<REAL>(t1 - t2);
    return diff / minutes_in_day; 
}

REAL date_act_365(const REAL& t1, const REAL& t2) { 
    return days_between(t1, t2) / 365.0; 
}

bool leap(const int& y) { 
    return  ( MOD(y,4) == 0  && ( ( MOD(y,100)!=0 ) || (MOD(y,400)==0) ) ); 
}


int end_of_month( const int& year, const int& month ) {
    bool leap_year = leap(year);
    if      ( month == 2 && leap_year )                               return 29;
    else if ( month == 2 )                                            return 28;
    else if ( month == 4 || month == 6 || month == 9 || month == 11 ) return 30;
    else                                                              return 31;
}

REAL add_months ( const REAL& time, const REAL& rnbmonths ) {
    int nbmonths = static_cast<int>( rnbmonths );
    Date date    = gregorian_of_date( static_cast<int> ( time ) );

    date.month += nbmonths;
    date.year  += (date.month-1) / 12;
    date.month  = MOD( date.month-1, 12 ) + 1;
 
    if ( date.month <= 0 ) { date.year -= 1; date.month += 12; }

    { // Result!
        int month_end = end_of_month(date.year, date.month);
        date.day = ( date.day < month_end ) ? date.day : month_end;
        date.hour = 12;
        date.min  = 0;
    }

    return static_cast<REAL>( date_of_gregorian( date ) );
}

REAL add_years( const REAL& date, const REAL& nbyears ) {
    return add_months( date, nbyears * 12.0 );
}

void test_dates() {
    REAL tmp;

    fprintf(stderr, "add_months(min_date,1)==48240 : ");
    tmp = add_months(MIN_DATE, 1.0);
    if( equalEps(tmp, 48240.0) )   fprintf(stderr, "SUCCESS!\n\n");
    else                           fprintf(stderr, "res: %f, FAILS!\n\n", tmp);

    fprintf(stderr, "add_months(min_date,2)==90000 : ");
    tmp = add_months(MIN_DATE, 2.0);
    if( equalEps(tmp, 90000.0) )   fprintf(stderr, "SUCCESS!\n\n");
    else                           fprintf(stderr, "res: %f, FAILS!\n\n", tmp);


    fprintf(stderr, "add_years(min_date,1)==530640 : ");
    tmp = add_years(MIN_DATE, 1.0);
    if( equalEps(tmp, 530640.0) )  fprintf(stderr, "SUCCESS!\n\n");
    else                           fprintf(stderr, "res: %f, FAILS!\n\n", tmp);

    fprintf(stderr, "add_years(min_date,5)==2634480 : ");
    tmp = add_years(MIN_DATE, 5.0);
    if( equalEps(tmp, 2634480.0) ) fprintf(stderr, "SUCCESS!\n\n");
    else                           fprintf(stderr, "res: %f, FAILS!\n\n", tmp);
    
    fprintf(stderr, "days_between(max_date(),min_date()) == 116877.499306 : ");
    tmp = days_between( MAX_DATE, MIN_DATE ) ;
    if( equalEps(tmp, 116877.499305555) ) fprintf(stderr, "SUCCESS!\n\n");
    else                                  fprintf(stderr, "res: %f, FAILS!\n\n", tmp);

    fprintf(stderr, "act_365(max_date,min_date) == 320.212327 : ");
    tmp = date_act_365( MAX_DATE, MIN_DATE );
    if( equalEps(tmp, 320.2123268645) )   fprintf(stderr, "SUCCESS!\n\n");
    else                                  fprintf(stderr, "res: %f, FAILS!\n\n", tmp);
}
#endif // ifndef GREG_CALEND

