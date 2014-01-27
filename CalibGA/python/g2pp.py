from __future__ import division # division among integers is lifted to floating point division
import numpy as N

import Swaption_quotes as SW
import Math
import Date
#from Date import *

# Zero coupon baseline
r=0.03
today = Date.of_string("2012-01-01")

infinity=1e49
epsilon=1e-5

_DEBUG=10
_GRAPH3D=False


def zc(t):
  # zc: zero coupon. Uses global vars: should be zc(t, rate, t_reference)
  return N.exp (-r * Date.act_365(t,today))

# Normal distribution: cumulative distribution function by linear Hermite expansion
x_quadrature, w_quadrature = Math.gauss_hermite_coefficients
nquadrature = len(x_quadrature)

def swap_schedule2lvl(swap_schedule,min_date=Date.min_date,max_date=Date.max_date):
  lvl=0
  t0=max_date
  tn=min_date
  for j,(tf,tp) in enumerate(swap_schedule):
    # reduce over the dimensions of swap_schedule:
    #   minimize the date: t0, maximize the date: tn, accumulate the zero coupon weights: lvl
    lvl+=zc(tp)*Date.act_365(tp,tf)
    t0=min(t0,tf)
    tn=max(tn,tp)
  return (lvl,t0,tn)

def extended_swaption_of_swaption(today,zc,swaption):
  swaption_maturity_in_year=swaption['swaption_maturity_in_year']
  swap_frequency=swaption['swap_frequency']
  swap_term_in_year=swaption['swap_term_in_year']
  #
  frq=float(swap_frequency)
  maturity = Date.add_years(today,swaption_maturity_in_year)
  nschedule=int(12*swap_term_in_year/frq)
  #
  swap_schedule=[]
  for j in xrange(nschedule):
    swap_schedule.append( ( Date.add_months(maturity,j*frq),Date.add_months(maturity,(j+1)*frq) ) )
  #
  (lvl,t0,tn)=swap_schedule2lvl(swap_schedule)
  strike = (zc(t0)-zc(tn))/lvl
  return {'maturity':maturity, 'swap_schedule':swap_schedule, 'strike':strike}

def b_fun(z,tau):
  return (1-N.exp(-z*tau))/z

def t_fun(sigma,x,tau):
  expxtau=N.exp(-x*tau)
  exp2xtau=expxtau*expxtau
  return sigma*sigma/(x*x)*(tau+2/x*expxtau-1/(2*x)*exp2xtau-3/(2*x))

def bigv(a,b,rho,nu,sigma,tau):
  """\var(\int_t^T [x(u)+y(u)]du)"""
  if sigma==0: sigma=1e-10 # sanity check
  ba=b_fun(a,tau)
  bb=b_fun(b,tau)
  t1=t_fun(sigma,a,tau)
  t2=t_fun(nu,b,tau)
  t3=2*rho*nu*sigma/(a*b)*(tau-ba-bb+b_fun(a+b,tau))
  return t1+t2+t3,ba,bb

def bigmx(a,b,rho,nu,sigma,today,tmat,s,t):
  """x drift term in tmat-forward measure"""
  ts=Date.act_365(t,s)
  tmatt=Date.act_365(tmat,t)
  tmat0=Date.act_365(tmat,today)
  tmats=Date.act_365(tmat,s)
  t0=Date.act_365(t,today)
  s0=Date.act_365(s,today)
  return ( sigma*sigma/(a*a)+sigma*rho*nu/(a*b) ) *\
    (1-N.exp(-a*ts)) -\
    (sigma*sigma/(2*a*a) * ( N.exp(-a*tmatt)-N.exp(-a*(tmats+ts)) ) ) -\
    rho*sigma*nu/(b*(a+b)) * ( N.exp(-b*tmatt)-N.exp(-b*tmat0-a*t0+(a+b)*s0) )

def bigmy(a,b,rho,nu,sigma,today,tmat,s,t):
  """y drift term in tmat-forward measure"""
  ts=Date.act_365(t,s)
  tmatt=Date.act_365(tmat,t)
  tmat0=Date.act_365(tmat,today)
  tmats=Date.act_365(tmat,s)
  t0=Date.act_365(t,today)
  s0=Date.act_365(s,today)
  return ( nu*nu/(b*b)+sigma*rho*nu/(a*b) ) *\
    (1-N.exp(-b*ts)) -\
    (nu*nu/(2*b*b) * ( N.exp(-b*tmatt)-N.exp(-b*(tmats+ts)) ) ) -\
    rho*sigma*nu/(a*(a+b)) * ( N.exp(-a*tmatt)-N.exp(-a*tmat0-b*t0+(a+b)*s0) )

def black_price(today,zc,swaption,vol):
  """Black'76 swaption pricer"""
  swaption=extended_swaption_of_swaption(today=today,zc=zc,swaption=swaption)
  swap_schedule=swaption['swap_schedule']
  strike=swaption['strike']
  maturity=swaption['maturity']
  #
  sqrtt=Date.act_365(maturity,today)
  (lvl,t0,tn)=swap_schedule2lvl(swap_schedule)
  #
  s0=(zc(t0)-zc(tn))/lvl
  d1=N.log(s0/strike)/(vol*sqrtt)+0.5*vol*sqrtt
  d2=d1-vol*sqrtt
  #
  return lvl*(s0*Math.uGaussian_P(d1)-strike*Math.uGaussian_P(d2))

"""
def a_fun(end_date,a,b,rho,nu,sigma,today,maturity,zc_mat,v0_mat):
  # Brigo and Mercurio: defined top p. 148
  v0_end,dummyA,dummyB=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,today))
  vt_end,ba,bb=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,maturity))
  res=zc(end_date)/zc_mat*N.exp(0.5*(vt_end-v0_end+v0_mat))
  return res,ba,bb
"""

def pricer_of_swaption(today,zc,swaption,params):
  swaption=extended_swaption_of_swaption(today=today, zc=zc, swaption=swaption)
  swap_schedule=swaption['swap_schedule']
  strike=swaption['strike']
  maturity=swaption['maturity']
  tmat0=Date.act_365(maturity,today)
  schedulei=swap_schedule
  taui=map(lambda (start_date,end_date): Date.act_365(end_date,start_date), schedulei)
  lastindex=len(schedulei)-1 # last index has different formula: updated below
  ci=map(lambda tau: tau*strike, taui)
  ci[-1]+=1
  #
  n_schedi=len(schedulei)
  bai=N.zeros(n_schedi)
  bbi=N.zeros(n_schedi)
  aici=N.zeros(n_schedi)
  log_aici=N.zeros(n_schedi)
  scales=N.zeros(n_schedi)
  t1_cst=N.zeros(n_schedi)
  scale=N.zeros(n_schedi)
  #
  a,b,rho,nu,sigma = paramsDict2tuple(params)
  v0_mat,dummyA,dummyB = bigv(a,b,rho,nu,sigma,tau=Date.act_365(maturity,today))
  zc_mat=zc(maturity)
  #
  def a_fun(end_date):
    """Brigo and Mercurio: defined top p. 148"""
    v0_end,dummyA,dummyB=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,today))
    vt_end,ba,bb=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,maturity))
    res=zc(end_date)/zc_mat*N.exp(0.5*(vt_end-v0_end+v0_mat))
    return res,ba,bb
  #
  sigmax=sigma*N.sqrt(b_fun(z=2*a,tau=tmat0))
  sigmay=   nu*N.sqrt(b_fun(z=2*b,tau=tmat0))
  rhoxy=rho*sigma*nu/(sigmax*sigmay) * b_fun(a+b,tmat0)
  rhoxyc=1-rhoxy*rhoxy
  rhoxycs=N.sqrt(rhoxyc)
  t2=rhoxy/(sigmax*rhoxycs)
  sigmay_rhoxycs=sigmay*rhoxycs
  t4=rhoxy*sigmay/sigmax
  #
  mux=-bigmx(a,b,rho,nu,sigma,today=today,tmat=maturity,s=today,t=maturity)
  muy=-bigmy(a,b,rho,nu,sigma,today=today,tmat=maturity,s=today,t=maturity)
  for i in xrange(n_schedi):
    #a_fun(end_date=schedulei[i][1],a=a,b=b,rho=rho,nu=nu,sigma=sigma,today=today,maturity=maturity,zc_mat=zc_mat,v0_mat=v0_mat)
    _a,ba,bb=a_fun(end_date=schedulei[i][1])
    x=ci[i]*_a
    log_ac=N.log(x)
    aici[i]=x
    log_aici[i]=log_ac
    bai[i]=ba
    bbi[i]=bb
    #
    t3=muy-0.5*rhoxyc*sigmay*sigmay*bb
    cst=bb*(mux*t4-t3)
    t1_cst[i]=x*N.exp(cst)
    scale[i]=-(ba+bb*t4)
  #
  k=-3.71901648545568 # ugaussian_Pinv(k)=1e-4
  # Root finder
  def exact_yhat(x):
    lo=-infinity
    up=0
    for i in xrange(n_schedi):
      baix=bai[i]*x
      lo=max(lo,(log_aici[i]-baix)/bbi[i])
      up+=aici[i]*N.exp(-baix)
    s_up=up
    if n_schedi==1:
      return lo
    else:
      log_s=N.log(s_up)
      tmp=log_s/bbi[n_schedi-1]
      if tmp<=0:
        up=tmp
      else:
        tmp=log_s/bbi[0]
        if 0<=tmp:
          up=tmp
        else:
          #This case happens when all ai * ci are too close to 0. or x to high => to_solve x y is close to -1 for all y,
          #       thus the root is reached for y negative with high absolute value (tosolve x is a decreasing function of y)
          #
          up=-infinity
      yl=lo-epsilon
      yu=up+epsilon
      # y01 x = y0, y1 / phi(h_i(x, y0)) <= epsilon, 1 - phi(h_i(x, y1)) <= epsilon
      y0=sigmay * (rhoxy*(x-mux)/sigmax+k*rhoxycs) - rhoxyc/b + muy
      y1=sigmay * (rhoxy*(x-mux)/sigmax-k*rhoxycs) + muy
      if      y1 <= yl:
        return y1+1 # yhat is greater than y1 => 1 - phi(h_i(x, yhat)) < epsilon
      else:
        if yu <= y0:
          return y0-1 # yhat is lower than y0 => phi(h_i(x, yhat)) < epsilon
        else:
          for i in xrange(n_schedi):
            scales[i]=aici[i] * N.exp(-bai[i]*x)
          # equation at bottom of page 158
          def to_solve(yhat):
            sum=-1
            for i in xrange(n_schedi):
              sum+=scales[i]*N.exp(-bbi[i]*yhat)
            return sum
          root_lb=max(yl,y0)
          root_ub=min(yu,y1)
          (root, iteration, error)=Math.rootFinding_Brent(
            f=to_solve,
            lb=root_lb,ub=root_ub,
            tol=1e-4,iteration_max=1000
            )
          #
          if root==None: # not bracketed?
            if error<0: return y0-1
            else: return y1+1
          else:
            return root
  #
  def yhat(x):
    eps=0.5*sigmax
    f=exact_yhat(mux)
    df=0.5*(exact_yhat(mux+eps) - exact_yhat(mux-eps))/eps
    return f + df*(x-mux)
  #
  def integrand(x):
    t1=N.exp(-0.5 * ((x-mux)/sigmax)**2)
    h1=( (yhat(x)-muy)/sigmay_rhoxycs ) - t2*(x-mux)
    _t2=Math.uGaussian_P(-h1)
    #
    acc=0
    for i in xrange(n_schedi):
      h2 = h1 + bbi[i] * sigmay_rhoxycs
      acc+=t1_cst[i] * N.exp (scale[i]*x) * Math.uGaussian_P (-h2)
    #
    return t1*(_t2-acc)
  #
  sqrt2sigmax = N.sqrt(2)*sigmax
  sum=0
  for i in xrange(nquadrature): # parallel map, then fold
    sum+=w_quadrature[i]*integrand(sqrt2sigmax*x_quadrature[i]+mux)
  integral=sum/N.sqrt(N.pi)
  #
  return zc_mat*integral

def params2dict(a, b, sigma, nu, rho):
  return {'a':a,'b':b,'sigma':sigma,'nu':nu,'rho':rho}

def paramsDict2tuple(params):
  return params['a'], params['b'], params['rho'], params['nu'], params['sigma']

def prices_plot(swaption_quotes=SW.swaption_quotes,today=today,zc=zc):
  import pylab as P
  #
  xlabel='swaption maturity in year'
  ylabel='swap term in year'
  zlabel='swaption price'
  title="Prices"
  #
  _m=SW.maturities(swaption_quotes=swaption_quotes)
  _t=SW.swap_terms(swaption_quotes=swaption_quotes)
  ax_x=_m
  ax_y=_t
  #
  quotes=N.array( map( lambda swp: (swp[0][xlabel.replace(" ","_")],swp[0][ylabel.replace(" ","_")],swp[1]), swaption_quotes) )
  prices_Black=prices_Black76(today=today,zc=zc,swaption_quotes=swaption_quotes)
  # Black'76 prices
  ax=SW.plot3Dwireframe(x=quotes[:,0],y=quotes[:,1],z=prices_Black,ax_x=_m,ax_y=_t,xlabel=xlabel,ylabel=ylabel,zlabel=zlabel,title=title)
  # G2++ prices
  params=params2dict(a = 0.02453, b = 0.98376, sigma = 0.02398, nu = 0.11830, rho = -0.82400)
  prices_G2PP=[]
  for j,(swaption,vol) in enumerate(swaption_quotes):
    if _DEBUG >=10:
      print "# swaption %3d/%3d: %s" % (j,len(swaption_quotes)-1,swaption)
    prices_G2PP.append( pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params) )
  #
  _x,_y=N.meshgrid(ax_x,ax_y)
  ax.plot_wireframe(_x,_y,N.array(prices_G2PP).reshape((len(ax_x),len(ax_y))),colors=[(1,0,0,1)]*len(prices_G2PP))
  ax.plot(quotes[:,0],quotes[:,1],prices_G2PP,"rd",label="g2++")
  P.legend()

def prices_Black76(today=today,zc=zc,swaption_quotes=SW.swaption_quotes):
  black76_prices=[]
  for (swaption,vol) in SW.swaption_quotes:
    black76_prices.append(black_price(today,zc,swaption,vol))
  return black76_prices

def test_swaptions(swaption_quotes=SW.swaption_quotes,today=today,zc=zc,scale=1e4):
  # Black'76 prices
  prices_Black=prices_Black76(today=today,zc=zc,swaption_quotes=swaption_quotes)
  # G2++ prices
  params=params2dict(a = 0.02453, b = 0.98376, sigma = 0.02398, nu = 0.11830, rho = -0.82400)
  for j,(swaption,vol) in enumerate(swaption_quotes):
    g2pp_price=pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)
    diff_perc=100*(g2pp_price-prices_Black[j])/prices_Black[j]
    print "# swaption %3d/%3d: %80s: g2pp: %8.3f, Black76: %8.3f (diff: %6.2f%%)" % (j,len(swaption_quotes)-1,swaption,g2pp_price*scale,prices_Black[j]*scale, diff_perc)

#=========================================================================
# SANITY CHECKS
#=========================================================================
# sanity checks with OCaml code
# b_fun 3.24 1.362
assert "%.6f" % b_fun(3.24,1.362)=="0.304901"
# bigv {g_a=0.02; g_b=0.02; g_sigma=0.04; g_nu=0.01; g_rho=0.} 1.12
assert "%.6f %.6f %.6f" % bigv(a=0.02, b=0.02, sigma=0.04, nu=0.01, rho=0.0, tau=1.12) == "0.000783 1.107549 1.107549"
# bigmx {g_a=0.02; g_b=0.02; g_sigma=0.04; g_nu=0.01; g_rho=0.} 9000 18000 400000 9000000
assert "%.6f" % bigmx(a=0.02, b=0.02, sigma=0.04, nu=0.01, rho=0.0, today=9000,tmat=18000,s=400000,t=9000000) == "-0.235607"
# bigmy {g_a=0.02; g_b=0.02; g_sigma=0.04; g_nu=0.01; g_rho=0.} 9000 18000 400000 9000000
assert "%.6f" % bigmy(a=0.02, b=0.02, sigma=0.04, nu=0.01, rho=0.0, today=9000,tmat=18000,s=400000,t=9000000) == "-0.014725"
# extended_swaption_of_swaption today zc {swaption_maturity_in_year = 10; swap_term_in_year = 4; swap_frequency = Freq_6M}
assert extended_swaption_of_swaption(today=today,zc=zc,swaption={'swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6}) == {'maturity':22094640, 'swap_schedule':[(22094640, 22355280), (22355280, 22620240), (22620240, 22880880), (22880880, 23145840), (23145840, 23407920), (23407920, 23672880), (23672880, 23933520), (23933520, 24198480)], 'strike':0.030226283149239714}
# black_price today zc {swaption_maturity_in_year = 10; swap_term_in_year = 4; swap_frequency = Freq_6M} 0.2454
swaption={'swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6}
assert "%.6f" % (black_price(today=today,zc=zc,swaption=swaption,vol=0.2454)*10000) == "654.142965"
# pricer_of_swaption today zc {swaption_maturity_in_year = 10; swap_term_in_year = 4; swap_frequency = Freq_6M} {g_a = 0.02453; g_b = 0.98376; g_sigma = 0.02398; g_nu = 0.11830; g_rho = -0.82400}
params=params2dict(a = 0.02453, b = 0.98376, sigma = 0.02398, nu = 0.11830, rho = -0.82400)
swaption={'swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6}
assert "%.3f" % (1e4*pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)) == "657.822"
swaption={'swaption_maturity_in_year': 30, 'swap_term_in_year': 30, 'swap_frequency': 6}
assert "%.3f" % (1e4*pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)) == "1902.976"

#=========================================================================
# MAIN
#=========================================================================
if __name__=='__main__':
  if _GRAPH3D:
    SW.swaption_plot(swaption_quotes=SW.swaption_quotes)
    prices_plot(swaption_quotes=SW.swaption_quotes,today=today,zc=zc)
    #
    P.show()
  else:
    print "# Date: today: %d, min_date: %d, max_date: %d" % (today,Date.min_date,Date.max_date)
    test_swaptions()
