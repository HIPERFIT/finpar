import datetime as D
import sys

"""
Date is an integer. Minimum resolution: minutes. Zero: 1979-12-30T00:00:00
A calendar day (YYYY-MM-DD) by default has HH:MM:SS = 12:00:00
This module uses datetime as a convenient tool to deal with Gregorian calendars (February 29th, months of different lengths, etc.)
  An ad-hoc, general implementation is present in the OCaml code
"""

date_zero=D.datetime(1979,12,30,0,0,0)

hours_in_day = 24.
minutes_in_day = hours_in_day * 60.
minutes_to_noon = minutes_in_day/2.
fminutes_in_day = float(minutes_in_day)

def days_between(t2,t1):
  return (t2 - t1) / fminutes_in_day

def act_365(t2,t1):
  return days_between(t2,t1)/365.

def _datetime2mins(d1,d0=date_zero):
  if sys.version_info < (2, 7): # legacy hack for python < 2.7
    dif = (d1-d0)
    return int( dif.seconds/60 + dif.days * 24 * 60)
  else:
    return int( (d1-d0).total_seconds()/60 )

def add_years(date,nbyears):
  return add_months(date, int(nbyears * 12))

def add_months(date,nbmonths):
  t0=D.timedelta(minutes=date)
  d0=date_zero+t0
  #
  d1_month=int(d0.month+nbmonths)%12
  d1_year=d0.year+int(nbmonths/12)
  d1=D.datetime(d1_year,d1_month,d0.day,d0.hour,d0.minute,d0.second)
  #
  return _datetime2mins(d1)

def of_string(date_str):
  """date_str: format: YYYY-MM-DDThh:mm:ss"""
  y,m,d=map(int,date_str.split("T")[0].split("-"))
  if "T" in date_str:
    hh,mm,ss=map(int,date_str.split("T")[1].split(":"))
  else:
    hh,mm,ss=12,0,0
  d1=D.datetime(y,m,d,hh,mm,ss)
  return _datetime2mins(d1)

max_date = of_string("2299-12-31T23:59:59")
min_date = of_string("1980-01-01T12:00:00")

# Sanity checks for compatibility with OCaml implementation
assert max_date==168307199
assert min_date==3600
assert add_months(min_date,1)==48240
assert add_months(min_date,2)==90000
assert add_years(min_date,1)==530640
assert add_years(min_date,5)==2634480
assert "%.6f" % days_between(max_date,min_date) == "116877.499306"
assert "%.6f" % act_365(max_date,min_date) == "320.212327"
