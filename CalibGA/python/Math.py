# https://en.wikipedia.org/wiki/Brent%27s_method
from __future__ import division # division among integers is treated as floating point
import numpy as N


_DEBUG=False

#-------------------------------------------------------------------------
# Cumulative Distribution Function for a standard normal distribution
def uGaussian_P(x):
  u = x/N.sqrt(2)
  if u<0: e = -erf(-u)
  else:   e = erf(u)
  return 0.5*(1+e)

# polynomial expansion of the erf() function, with error<=1.5e-7
  # formula 7.1.26 (page 300), Handbook of Mathematical Functions, Abramowitz and Stegun
    # http://people.math.sfu.ca/~cbm/aands/frameindex.htm
def erf(x):
  p = 0.3275911
  a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
  t=1/(1+p*x)
  t2=t*t
  t3=t*t2
  t4=t2*t2
  t5=t2*t3
  return 1 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * N.exp(-(x*x))
#-------------------------------------------------------------------------
def brent(a,b,eval,eps=1e-9):
  return rootFinding_Brent(f=eval,lb=a,ub=b,tol=eps,iteration_max=10000)

def rootFinding_Brent(f,lb,ub,tol=1e-9,iteration_max=1000):
  infinity=1e49
  # Brent's method:
  a,b=lb,ub
  fa,fb=f(a),f(b)
  if fa*fb>=0:
    # root is not bracketed
    if fa>=0:
      return (None, 0, +infinity) # root not bracketed above
    else:
      return (None, 0, -infinity) # root not bracketed below
  if abs(fa)<abs(fb):
    # swap arguments
    fa,fb=fb,fa
    a,b=b,a

  # initialization
  c,fc=a,fa
  mflag=True
  iteration=0
  while not( fb==0 or abs(b-a)<tol ):
    if fa!=fc and fb!=fc:
      # inverse quadratic interpolation
      s1=(a*fb*fc)/( (fa-fb)*(fa-fc) )
      s2=(b*fa*fc)/( (fb-fa)*(fb-fc) )
      s3=(c*fa*fb)/( (fc-fa)*(fc-fb) )
      s=s1+s2+s3
    else:
      # secant rule
      s=b-fb*(b-a)/(fb-fa)
    if (
        not (s>=(3*a+b)/4 and s<=b)
        or (mflag==True and abs(s-b)>=abs(b-c)/2)
        or (mflag==False and abs(s-b)>=abs(c-d)/2)
        or (mflag==True and abs(b-c)<=abs(tol))
        or (mflag==False and abs(c-d)<=abs(tol))
      ):
      s=(a+b)/2
      mflag=True
    else:
      mflag=False
    fs=f(s)
    # d is assigned for the first time here: it's not used above because mflag is set
    d=c
    c,fc=b,fb
    if fa*fs<0:
      b,fb=s,fs
    else:
      a,fa=s,fs
    #
    if abs(fa)<abs(fb):
      # swap arguments
      fa,fb=fb,fa
      a,b=b,a
    #
    iteration+=1
    if iteration>=iteration_max:
      print "# ERROR: Brent's method not converged at iteration %d: error is: %.3f" % (iteration,fb)
      break
    #
    if _DEBUG==True:
      print "_Brent_Loop: a: %.2f, b: %.2f, c: %.2f, mflag: %s, s: %.2f, iter: %d" % (a,b,c,mflag,s,iteration)
  #
  error=fb
  root=b
  return (root, iteration, error)

#-------------------------------------------------------------------------
# Gaussian Quadrature with Hermite linear expansion: cmulative distribution function of Normal distribution
  
gauss_hermite_coefficients = \
  [
    # coefficients
    [
      0., 0.6568095668820999044613, -0.6568095668820997934390, -1.3265570844949334805563,  1.3265570844949330364670,  2.0259480158257567872226,
      -2.0259480158257558990442, -2.7832900997816496513337,  2.7832900997816474308877,  3.6684708465595856630159, -3.6684708465595838866591
    ],
    # weights
    [
      0.6547592869145917315876, 0.6609604194409607336169, 0.6609604194409606225946, 0.6812118810666693002887, 0.6812118810666689672217, 0.7219536247283847574252,
      0.7219536247283852015144, 0.8025168688510405656800, 0.8025168688510396775015, 1.0065267861723647957461, 1.0065267861723774522886
    ]
  ]
#=========================================================================
def equal(x1,x2,tol=1e-8):
  return N.abs(x1-x2)<=tol

def test():
  # Rootfinder.brent (-4.) (4./.3.) (fun x -> (x+.3.)*.(x-.1.)**2.) 1e-4 == -3
  print "# Brent test: ", equal(brent( a=-4, b=4/3, eval=lambda x: (x+3)*(x-1)**2 )[0], -3)
  # erf 0. == 0. ;; 100. *. erf (1./.sqrt 2.)
  print "# Erf test:", equal(erf(0),0) and equal(int(100*erf(1/N.sqrt(2))),68)
  # ugaussian_P 0. ;; ugaussian_P 1. +. ugaussian_P (-1.)
  print "# Gaussian test: ", equal(uGaussian_P(0),0.5) and equal(uGaussian_P(-1)+uGaussian_P(1),1)
#=========================================================================
if __name__=='__main__':
  test()
