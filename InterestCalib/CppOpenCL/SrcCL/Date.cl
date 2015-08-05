//////////////////////////////////////////////////////
////   Date: Gregorian calendar
//////////////////////////////////////////////////////
#ifndef GREG_CALEND
#define GREG_CALEND

#define MINUTES_IN_DAY  1440   // hours_in_day * 60;
#define MINUTES_TO_NOON 720    // (hours_in_day / 2) * 60;

typedef struct DDD{
    int year;
    int month;
    int day;
    int hour;
    int min;
} Date;

Date mkDate(int y, int mo, int d, int h, int mn) {
    Date date;
    date.year  = y; 
    date.month = mo; 
    date.day   = d; 
    date.hour  = h; 
    date.min   = mn;
    return date;
}

inline
int MOD(int x, int y) { 
    return ( x - (x/y)*y );
}

inline
int date_of_gregorian( Date date ) {
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
    
    return tmp * MINUTES_IN_DAY + date.hour * 60 + date.min;
}

inline
Date gregorian_of_date ( int minutes_since_epoch ) {
    int jul, l, n, i, j, d, m, y, daytime;

    jul = minutes_since_epoch / MINUTES_IN_DAY;
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

    daytime = MOD( minutes_since_epoch, MINUTES_IN_DAY );

//    if ( daytime == MINUTES_TO_NOON )  return Date(y,m,d,12,0);
//    else                               return Date(y,m,d, daytime / 60, MOD(daytime, 60));

    return mkDate(y,m,d, daytime / 60, MOD(daytime, 60));
}

inline
bool check_date(int year, int month, int day) {
    bool tmp1, tmp2, tmp3;

    tmp1 = ( 1 <= day && 1 <= month && month <= 12 && 1980 <= year && year <= 2299 );
    tmp2 = ( day <= 28 );

    if        ( month == 2 ) {
        tmp3 = ( day == 29 && MOD(year, 4) == 0 && ( year == 2000 || (! (MOD(year, 100) == 0)) ) );
    } else if ( month == 4 || month == 6 || month == 9 || month == 11 ) {
        tmp3 = ( day <= 30 );
    } else {
        tmp3 = ( day <= 31 );
    }

    return tmp1 && (tmp2 || tmp3);
}

inline
REAL days_between(REAL t1, REAL t2) {
    return ( (t1 - t2) / MINUTES_IN_DAY ); 
}

inline
REAL date_act_365(REAL t1, REAL t2) { 
    return ( days_between(t1, t2) / 365.0 ); 
}

inline
bool leap(int y) { 
    return  ( MOD(y,4) == 0  && ( ( MOD(y,100)!=0 ) || (MOD(y,400)==0) ) ); 
}

inline
int end_of_month( int year, int month ) {
    bool leap_year = leap(year);
    if      ( month == 2 && leap_year )                               return 29;
    else if ( month == 2 )                                            return 28;
    else if ( month == 4 || month == 6 || month == 9 || month == 11 ) return 30;
    else                                                              return 31;
}

inline
REAL add_months ( REAL time, REAL rnbmonths ) {
    int nbmonths = convert_int(rnbmonths);   //(int) rnbmonths;
    Date date    = gregorian_of_date( convert_int( time ) );

    date.month += nbmonths;
    date.year  += (date.month-1) / 12;
    date.month  = MOD( date.month-1, 12 ) + 1;
 
    if ( date.month <= 0 ) { date.year -= 1; date.month += 12; }

    { // date.day = min( date.day, end_of_month(date.year, date.month) );
        int tmp = end_of_month(date.year, date.month);
        if(date.day > tmp) date.day = tmp;
    }

    date.hour = 12;
    date.min  = 0;

    return (REAL) ( date_of_gregorian( date ) );
}

inline
REAL add_years( REAL date, REAL nbyears ) {
    return add_months( date, nbyears * 12.0 );
}

#endif // ifndef GREG_CALEND

