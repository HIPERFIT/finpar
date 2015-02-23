module MathMod 
        (
          uGaussianP
        , uGaussianPwithExpFactor
        , zc
        , bigV
        , bigmX
        , bigmY
        , bFun
        , rootFindBrent
        , swapSchedule
        , blackPrice
        , logLikelihood
        , testMathMod
        , testG2ppUtil
        ) where

import Constants
import Date

import Debug.Trace

--let (p, a1, a2, a3, a4, a5) = 
--            (0.3275911, 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429)

---------------------------------------------------------------------------
--- polynomial expansion of the erf() function, with error<=1.5e-7
---  formula 7.1.26 (page 300), Handbook of Mathematical Functions, Abramowitz and Stegun
---  http://people.math.sfu.ca/~cbm/aands/frameindex.htm
---------------------------------------------------------------------------
erffPolyOnly :: Double -> Double
erffPolyOnly x = 
    let (p, a1, a2, a3, a4, a5) = 
            (0.3275911, 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429)
        t  = 1.0/(1.0+p*x)
        t2 = t  * t 
        t3 = t  * t2
        t4 = t2 * t2
        t5 = t2 * t3
    in  (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5)

erff1 :: Double -> Double
erff1 x = 
    let poly = erffPolyOnly x
    in  ( 1.0 - poly * exp(0.0-(x*x)) )

-- Cumulative Distribution Function for a standard normal distribution
uGaussianP :: Double -> Double
uGaussianP x = 
    let u = x / sqrt(2.0)
        e = if u < 0.0
            then 0.0 - (erff1 (0.0-u))
            else erff1 u
    in  ( 0.5 * (1.0 + e) )

uGaussianPwithExpFactor :: Double -> Double -> Double
uGaussianPwithExpFactor x exp_factor = 
    let u   = abs ( x / (sqrt 2.0) )
        e   = erffPolyOnly u
        res = 0.5 * e * ( exp (exp_factor - u*u) )
    in  if x >= 0.0
        then (exp exp_factor) - res
        else res

---------------------------------------------------------------------------
-- THE FUNction-PARAMeter to rootFinding_Brent
-- if fid = 33 then for testing to_solve(real x) = (x+3.0)*(x-1.0)*(x-1.0)
-- Otherwise follows the real implementation
-- N is the size of scales (bbi)
---------------------------------------------------------------------------
toSolve :: [Double] -> [Double] -> Double -> Double
toSolve scales bbs yhat = 
    let res_map = zipWith ( \scale bb -> scale * (exp (0.0 - bb * yhat)) ) scales bbs 
        res_red = reduce (+) 0.0 res_map
    in  res_red - 1.0

toSolveTest :: [Double] -> [Double] -> Double -> Double
toSolveTest _ _ yhat = 
    (yhat+3.0)*(yhat-1.0)*(yhat-1.0)

max_num_brent_iters :: Int
max_num_brent_iters = 10000

rootFindBrentGen ::   ([Double] -> [Double] -> Double -> Double)         -- toSolve function
                   -> [Double] -> [Double] -> Double -> Double -> Double -- scales bbs lb ub tolerance 
                   -> (Double, Int, Double)                              -- Result: (root, # needed iter, fb)
rootFindBrentGen toSolveFun scales bbs lb ub toll = 
    let tol = if toll <= 0.0 then 1.0e-9 else toll
        (a, b) = (lb, ub)
        fa  = toSolveFun scales bbs a
        fb  = toSolveFun scales bbs b
  
    in  if fa*fb >= 0.0
        then (0.0, 0, if a >= 0.0 then infty else 0.0-infty)
        else 
        let (a', fa', b', fb') = if (abs fa) < (abs fb) 
                                 then (b, fb, a, fa) 
                                 else (a, fa, b, fb)            
        in  brentWhileLoop (toSolveFun scales bbs) 
                           tol 0 True 0.0 a' fa' b' fb' a' fa' 

rootFindBrent ::   [Double] -> [Double]       -- scales bbs 
                -> Double -> Double -> Double -- lb ub tolerance 
                -> (Double, Int, Double)
rootFindBrent = rootFindBrentGen toSolve

brentWhileLoop ::   (Double -> Double) -- toSolve function
                 -> Double -> Int -> Bool -> Double
                 -> Double -> Double -> Double -> Double -> Double -> Double
                 -> (Double, Int, Double)
brentWhileLoop toSolveFun tol i mflag d a fa b fb c fc =
    if (fb == 0.0) || (abs (b-a)) < tol || i >= max_num_brent_iters
    then (b, i, fb)
    else let s = if (fa == fc) || (fb == fc) 
                 then b - fb * (b - a) / (fb - fa)
                 else let s1 = ( a*fb*fc ) / ( (fa-fb)*(fa-fc) )
                          s2 = ( b*fa*fc ) / ( (fb-fa)*(fb-fc) )
                          s3 = ( c*fa*fb ) / ( (fc-fa)*(fc-fb) )
                      in  s1 + s2 + s3
             cond = (    (3.0 * a + b) / 4.0 > s || s > b )            ||
                    (     mflag  && (abs (b-c)) / 2.0 <= (abs (s-b)) ) ||
                    ( not mflag  && (abs (c-d)) / 2.0 <= (abs (s-b)) ) ||
                    (     mflag  && (abs (b-c))       <= (abs tol)   ) ||
                    ( not mflag  && (abs (c-d))       <= (abs tol)   )
             (mflag', s') = if cond 
                            then ( True , (a+b)/2.0 ) 
                            else ( False, s         )

             fs' = toSolveFun s'
            
             (d', c', fc') = (c, b, fb)
             (a', fa', b', fb') = if fa*fs' < 0.0
                                  then (a, fa, s', fs')
                                  else (s', fs', b, fb)
             (a'', fa'', b'', fb'') = 
                                  if (abs fa') < (abs fb')
                                  then (b', fb', a', fa')
                                  else (a', fa', b', fb')

         in  brentWhileLoop toSolveFun tol (i+1) mflag' d' a'' fa'' b'' fb'' c' fc'

testMathMod :: Int -> Int
testMathMod arg0 = 
    let (scales, bbs) = ([ 0.0 ], [ 0.0 ])

        -- Rootfinder.brent (fun x -> (x+.3.)*.(x-.1.)**2.) (-4.) (4./.3.) 1e-4 == -3
        (tmp1,_,_) = trace ("# Brent test: -3.0 => ") $ 
                        rootFindBrentGen toSolveTest scales bbs (0.0-4.0) (4.0/3.0) 0.0
        arg1 = if equalEps tmp1 (0.0 - 3.0)
               then trace ("... SUCCESS !") arg0
               else trace ("... FAILS   !") arg0
        
        -- erff 0. == 0. ;; 100. *. erff (1./.sqrt 2.)
        tmp2 = trace ("# Erf test => ") $
                    100.0 * ( erff1 (1.0 / (sqrt 2.0)) )
        arg2 = if (equalEps (erff1 0.0) 0.0) && (equalEps tmp2 68.26894723352726)
               then trace ("... SUCCESS !") arg1
               else trace ("... FAILS   !") arg1

        -- ugaussian_P 0. ;; ugaussian_P 1. +. ugaussian_P (-1.)
        tmp3 = trace ("# Gaussian test => ") $
                    (uGaussianP (0.0-1.0)) + (uGaussianP 1.0)
        arg3 = if (equalEps (uGaussianP 0.0) 0.5) && (equalEps tmp3 1.0)
               then trace ("... SUCCESS !") arg2
               else trace ("... FAILS   !") arg3
    in  arg3

--------------------------------------
--------------------------------------
rrr :: Double
rrr =  0.03

zc :: Double -> Double
zc t = exp (0.0 - rrr * (dateAct365 t today) )

bFun :: Double -> Double -> Double
bFun z0 tau = ( 1.0 - (exp (0.0 - z0*tau)) ) / z0

tFun :: Double -> Double -> Double -> Double
tFun sigma x0 tau = 
    let expxtau  = exp ( 0.0 - x0*tau )
        exp2xtau = expxtau*expxtau
    in  sigma*sigma/(x0*x0) * ( tau + 2.0/x0*expxtau-1.0/(2.0*x0)*exp2xtau-3.0/(2.0*x0) )

-- the first parameter `genome' is the five-genes genome used in
--     the genetic algorithms that models the interest rate
-- the second parameter is the time
-- the result is V in Brigo and Mercurio's book page 148.
--     \var(\int_t^T [x(u)+y(u)]du)
bigV ::    Genome -> Double 
        -> (Double,Double,Double)  
bigV (a,b,rho,nu,sigma) tau = 
    let bai = bFun a       tau
        bbi = bFun b       tau
        v   = tFun sigma a tau
        o1  = tFun nu    b tau
        o2  = 2.0 * rho * nu * sigma / (a * b)*
             ( tau - bai - bbi + (bFun (a+b) tau) )
    in  (v + o1 + o2, bai, bbi)


-- the first parameter `genome' is the five-genes genome used in
--     the genetic algorithms that models the interest rate
-- the other parameter are times: today, maturity, and the
--      lower and upper bound of the considered time interval
--
-- the result is: x drift term in tmat-forward measure

bigmX ::   Genome -> Double -> Double -> Double -> Double
        -> Double
bigmX (a,b,rho,nu,sigma) iday tmat s t = 
    let ts    = dateAct365 t    s
        tmatt = dateAct365 tmat t

        tmat0 = dateAct365 tmat iday
        tmats = dateAct365 tmat s
        
        t0    = dateAct365 t    iday
        s0    = dateAct365 s    iday

        tmp0  = (sigma*sigma)/(a*a)+(sigma*rho*nu)/(a*b)
        tmp1  = 1.0 - (exp (0.0 - a*ts))

        tmp2  = sigma * sigma / (2.0 * a * a)
        tmp3  = (exp (0.0 - a*tmatt)) - (exp (0.0 - a*(tmats+ts)))

        tmp4  = rho * sigma * nu / (b * (a + b))
        tmp5  = (exp (0.0 - b*tmatt)) - (exp (0.0 -b*tmat0 - a*t0 + (a+b)*s0))

    in  tmp0*tmp1 - tmp2*tmp3 - tmp4*tmp5


-- the first parameter `genome' is the five-genes genome used in
--     the genetic algorithms that models the interest rate
-- the other parameter are times: today, maturity, and the
--      lower and upper bound of the considered time interval
--
-- the result is: y drift term in tmat-forward measure
bigmY ::   Genome -> Double -> Double -> Double -> Double
        -> Double
bigmY (a,b,rho,nu,sigma) iday tmat s t = 
    let ts    = dateAct365 t    s
        tmatt = dateAct365 tmat t
        tmat0 = dateAct365 tmat iday
        tmats = dateAct365 tmat s
        t0    = dateAct365 t    iday
        s0    = dateAct365 s    iday

        tmp0  = nu*nu/(b*b)+sigma*rho*nu/(a*b)
        tmp1  = 1.0 - (exp (0.0 - b*ts))

        tmp2  = nu * nu / (2.0 * b * b)
        tmp3  = (exp (0.0 -b*tmatt)) - (exp (0.0 - b*(tmats + ts)))

        tmp4  = sigma * rho * nu / (a * (a + b))
        tmp5  = (exp (0.0 - a*tmatt)) - (exp (0.0 - a*tmat0 - b*t0 + (a+b)*s0))

    in  tmp0*tmp1 - tmp2*tmp3 - tmp4*tmp5



swapSchedule ::    Int -> Double -> Double 
                -> (Double, Double, Double)
swapSchedule nschedule maturity freq = 
    let lvlt0tns = 
            map (\j -> let  i   = fromIntegral j
                            a1  = addMonths maturity (freq*i)
                            a2  = addMonths a1        freq
                            lvl = (zc a2) * (dateAct365 a2 a1)
                       in   (lvl, a1, a2)
                ) [0..nschedule-1]
        -- Reduction( lvl: +, t0 : min, tn : max )
    in  reduce  (\ (l,ti,tf) (lvl,tini, tfin) ->
                    (l+lvl, min ti tini, max tf tfin) 
                ) (0.0, max_date, min_date) lvlt0tns


blackPrice ::   Double -> Swaption -> Double
blackPrice iday (sw_mat, sw_freq, sw_ty, sw_vol) = 
    let maturity    = addYears iday sw_mat
        sqrtt       = dateAct365 maturity iday
        n_schedule  = truncate $ 12.0 * sw_ty / sw_freq
        (lvl,t0,tn) = swapSchedule n_schedule maturity sw_freq
        
        strike = ( (zc t0) - (zc tn) ) / lvl
        d1     = 0.5 * sw_vol * sqrtt
        d2     = 0.0 - d1
    in  ( lvl * ( strike * (uGaussianP d1) - strike * (uGaussianP d2) ) )


testG2ppUtil :: Int -> Int
testG2ppUtil arg0 = 
    let (iday, tmat, s, t) = (9000.0, 18000.0, 400000.0, 9000000.0)
        
        res_b_fun = trace ("# bFun test = 0.30490117 => ") $ 
                        bFun 3.24 1.362
        arg1 = if equalEps res_b_fun 0.30490117
               then trace ("... SUCCESS !") arg0
               else trace ("... FAILS   !") arg0

        res_bigmX = trace ("# bigmX test = -0.2356067470979 => ") $ 
                        bigmX (0.02, 0.02, 0.0, 0.01, 0.04) iday tmat s t
        arg2 = if equalEps res_bigmX (-0.2356067470979)
               then trace ("... SUCCESS !") arg1
               else trace ("... FAILS   !") arg1

        res_bigmY = trace ("# bigmY test = -0.01472542169362 => ") $ 
                        bigmY (0.02, 0.02, 0.0, 0.01, 0.04) iday tmat s t
        arg3 = if equalEps res_bigmY (-0.01472542169362)
               then trace ("... SUCCESS !") arg2
               else trace ("... FAILS   !") arg2

        (res_bigv1, res_bigv2, res_bigv3) = 
                trace ("# bigv test  = { 7.8288965347e-4, 1.107549139, 1.107549139 } => ") $ 
                        bigV (0.02, 0.02, 0.0, 0.01, 0.04) 1.12
        arg4 = if   (equalEps res_bigv1 7.8288965347e-4) &&
                    (equalEps res_bigv2 1.107549139)     && 
                    (equalEps res_bigv3 1.107549139)       
               then trace ("... SUCCESS !") arg3
               else trace ("... FAILS   !") arg3

        (sw_mat,freq,sw_ty) = (10.0, 6.0, 4.0) 
                    
        (maturity, strike)  = (22094640.0, 0.030226283149239714)
        swap_schedule1      = [ 22094640.0, 22355280.0, 22620240.0, 22880880.0, 
                                23145840.0, 23407920.0, 23672880.0, 23933520.0 ]
        swap_schedule2      = [ 22355280.0, 22620240.0, 22880880.0, 23145840.0, 
                                23407920.0, 23672880.0, 23933520.0, 24198480.0 ]

        vol = 0.2454
        black_price_res = 
                trace ("# Testing Black Price = 654.1429648454 => ") $ 
                        (blackPrice today (sw_mat, freq, sw_ty, vol)) * 10000.0
        arg5 = if   (equalEps black_price_res 654.1429648454) 
               then trace ("... SUCCESS !") arg4
               else trace ("... FAILS   !") arg4
    in  arg5


--------------------------------------
------- Likelihood Measures ----------
--------------------------------------

sqrtTwoPi :: Double
sqrtTwoPi = sqrt (2*pi)

data Likelihood_Type = CAUCHY | NORMAL

llhood_type :: Likelihood_Type
llhood_type = CAUCHY -- NORMAL;
 
llhood_cauchy_offs :: Double
llhood_cauchy_offs = 5.0

llhood_normal_offs :: Double
llhood_normal_offs = 1.0



normalPdf :: Double -> Double -> Double -> Double
normalPdf z mu sigma = 
    let sigmap = abs sigma
        res    = 1.0 / (sigmap * sqrtTwoPi)
        ecf    = (z - mu) * (z - mu) / (2.0 * sigmap * sigmap)
    in  res * (exp (0.0 - ecf))

cauchyPdf :: Double -> Double -> Double -> Double
cauchyPdf z mu gamma = 
    let x = (z - mu) / gamma
    in  1.0 / ( pi * gamma * ( 1.0 + x*x ) )

logLikelihoodNormal :: Double -> Double -> Double
logLikelihoodNormal y_ref y =
    let sigma = (y_ref / 50.0) * llhood_normal_offs
        pdfs  = normalPdf y y_ref sigma
    in  log (1.0e-20 + pdfs)

logLikelihoodCauchy :: Double -> Double -> Double
logLikelihoodCauchy y_ref y =
    let gamma = ( (abs y_ref) / 50.0 ) * llhood_cauchy_offs + 0.01
        pdfs  = cauchyPdf y y_ref gamma
    in  log (1.0e-20 + pdfs)


logLikelihood :: Double -> Double -> Double
logLikelihood y_ref y = 
    case llhood_type of
        NORMAL -> logLikelihoodNormal y_ref y
        CAUCHY -> logLikelihoodCauchy y_ref y

