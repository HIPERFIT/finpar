module Date (
              today
            , max_date
            , min_date
            , dateAct365
            , addMonths
            , addYears
            , test_dates
            ) where

import Constants
import Debug.Trace

type GregDate = (Int,Int,Int,Int,Int) -- (year,month,day,hour,min)

today :: Double
today = fromIntegral $ date_of_gregorian (2012, 1, 1, 12, 0)

--"2299-12-31T23:59:59"
max_date :: Double
max_date = 168307199.0

--"1980-01-01T12:00:00"
min_date :: Double
min_date = 3600.0

hours_in_day :: Int
hours_in_day    = 24
minutes_in_day  = hours_in_day * 60
minutes_to_noon = (hours_in_day `div` 2) * 60


date_of_gregorian :: GregDate -> Int
date_of_gregorian (year,month,day,hour,min) = 
    let ym = if month == 1 || month == 2
             then   ( 1461 * ( year  + 4800 - 1 ) ) `div` 4  +
                    ( 367  * ( month + 10       ) ) `div` 12 -
                    ( 3  * ( ( year  + 4900 - 1 ) `div` 100 ) ) `div` 4

             else   ( 1461 * ( year  + 4800 ) ) `div` 4  +
                    ( 367 *  ( month - 2    ) ) `div` 12 -
                    ( 3  * ( ( year  + 4900 ) `div` 100 ) ) `div` 4

        tmp = ym + day - 32075 - 2444238

    in  tmp * minutes_in_day + hour * 60 + min


gregorian_of_date :: Int -> GregDate
gregorian_of_date minutes_since_epoch = 
    let jul = minutes_since_epoch `div` minutes_in_day
        l1  = jul + 68569 + 2444238
        n   = ( 4 * l1 ) `div` 146097
        l2  = l1 - ( 146097 * n + 3 ) `div` 4
        i   = ( 4000 * ( l2 + 1 ) ) `div` 1461001

        l3  = l2 - ( 1461 * i ) `div` 4 + 31
        j   = ( 80 * l3 ) `div` 2447
        d   = l3 - ( 2447 * j ) `div` 80
        l   = j `div` 11

        m   = j + 2 - ( 12 * l )
        y   = 100 * ( n - 49 ) + i + l

        daytime = mod minutes_since_epoch minutes_in_day

    in  if ( daytime == minutes_to_noon )  
        then  (y, m, d, 12, 0)
        else  (y, m, d, daytime `div` 60, mod daytime 60)

check_date :: Int -> Int -> Int -> Bool
check_date year month day = 
    let tmp1 = (1 <= day) && (1 <= month) && (month <= 12) && (1980 <= year) && (year <= 2299)
        tmp2 = (day <= 28)
        tmp3 =  if (month == 2)
                then (day == 29) && ((mod year 4) == 0) && ( year == 2000 || (not ((mod year 100) == 0)) )
                else if ( month == 4 || month == 6 || month == 9 || month == 11 )
                then ( day <= 30 )
                else ( day <= 31 )
    in  tmp1 && (tmp2 || tmp3)


days_between :: Double -> Double -> Double
days_between t1 t2 = (t1 - t2) / (fromIntegral minutes_in_day)

dateAct365 :: Double -> Double -> Double
dateAct365 t1 t2 = (days_between t1 t2) / 365.0

leap :: Int -> Bool
leap y = ( (mod y 4) == 0 ) && ( not ( (mod y 100) == 0 ) || ( (mod y 400) == 0 ) )

end_of_month :: Int -> Int -> Int
end_of_month year month = 
    let leap_year = leap year
    in  if      ( month == 2 && leap_year ) then 29
        else if ( month == 2 )              then 28
        else if ( month == 4 || month == 6 || month == 9 || month == 11 ) then 30
        else 31

addMonths :: Double -> Double -> Double
addMonths time rnbmonths = 
    let nbmonths = truncate rnbmonths
        (year,month,day,hour,min) = gregorian_of_date (truncate time)
        
        month1  = month + nbmonths
        year''  = year  + (month1 - 1) `div` 12
        month'' = (mod (month1 - 1) 12) + 1
        
        (year''', month''') = if (month'' <= 0)
                              then (year''-1, month'' + 12)
                              else (year'',   month'')
        month_end = end_of_month year''' month'''
        day'''    = if ( day < month_end ) 
                    then day 
                    else month_end
    in  fromIntegral $ date_of_gregorian (year''', month''', day''', 12, 0)


addYears :: Double -> Double -> Double
addYears date nbyears = addMonths date (nbyears * 12.0)

test_dates :: Int -> Int
test_dates arg =
    let tmp1 = trace ("addMonths(min_date,1)==48240 : ") $ addMonths min_date 1.0
        arg1 = if equalEps tmp1 48240.0
               then trace ("... SUCCESS !") arg
               else trace ("... FAILS   !") arg

        tmp2 = trace ("addMonths(min_date,2)==90000 : ") $ addMonths min_date 2.0
        arg2 = if equalEps tmp2 90000.0
               then trace ("... SUCCESS !") arg1
               else trace ("... FAILS   !") arg1

        tmp3 = trace ("addYears(min_date,1)==530640 : ") $ addYears min_date 1.0
        arg3 = if equalEps tmp3 530640.0
               then trace ("... SUCCESS !") arg2
               else trace ("... FAILS   !") arg2

        tmp4 = trace ("addYears(min_date,5)==2634480 : ") $ addYears min_date 5.0
        arg4 = if equalEps tmp4 2634480.0
               then trace ("... SUCCESS !") arg3
               else trace ("... FAILS   !") arg3

        tmp5 = trace ("days_between(max_date(),min_date()) == 116877.499306 : ") $ days_between max_date min_date
        arg5 = if equalEps tmp5 116877.499305555
               then trace ("... SUCCESS !") arg4
               else trace ("... FAILS   !") arg4

        tmp6 = trace ("act_365(max_date,min_date) == 320.212327 : ") $ dateAct365 max_date min_date
        arg6 = if equalEps tmp6 320.2123268645
               then trace ("... SUCCESS !") arg5
               else trace ("... FAILS   !") arg5
        
    in  arg6

