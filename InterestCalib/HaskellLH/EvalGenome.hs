module EvalGenome 
        (
          evalGenome
        , testSwaptionPricer
        ) where

import Data.List
import Data.Vector.Unboxed(Vector) 
import qualified Data.Vector.Unboxed as V
import Control.DeepSeq

import Constants
import Date
import MathMod

import Debug.Trace



----------------------------------
--- Main Export of this Module ---
----------------------------------
evalGenome ::   [Double] -> [Double] -> [Swaption] -> Vector Double
             -> ( Double, Vector Double, Vector Double )   -- (logLik,calib_prices, black_prices) [NUM_SWAP_QUOTES]  
evalGenome hermite_coeffs hermite_weights swaption_quotes genome =
    let genome_tup = -- (a, b, rho, nu, sigma) = 
            (genome V.! 0, genome V.! 1, genome V.! 2, genome V.! 3, genome V.! 4)

        (black_prices, swap_prices, swap_logLiks) = unzip3 $
            map (pricerOfSwaption hermite_coeffs hermite_weights genome_tup) swaption_quotes

        swap_logLik    = reduce (+) 0.0 swap_logLiks
        swap_price_vct = V.fromList swap_prices
        black_price_vct= V.fromList black_prices
        
    in  swap_logLik `deepseq` swap_price_vct `deepseq` black_price_vct `deepseq`
        ( swap_logLik, swap_price_vct, black_price_vct )
        -- ( swap_logLik, V.fromList swap_prices, V.fromList black_prices )


pricerOfSwaption ::  [Double] -> [Double] -> Genome -> Swaption 
                  -> (Double, Double, Double)        -- (Black Price, Swapt Price, LogLikelihood)
pricerOfSwaption hermite_coeffs hermite_weights      -- Hermite coefficients & weights
                genome (sw_mat, sw_freq, sw_ty, sw_quote) = -- swaption

    let (a, b, rho, nu, sigma) = genome
        maturity   = addYears today sw_mat
        n_schedi   = truncate $ 12.0 * sw_ty / sw_freq
        tmat0      = dateAct365 maturity today

        -- 1. compute quote, i.e., black price inlined
        (lvl,t0,tn) = swapSchedule n_schedi maturity sw_freq
        strike = ( (zc t0) - (zc tn) ) / lvl
        d1     = 0.5 * sw_quote * tmat0
        d2     = 0.0 - d1
        black_price = ( lvl * ( strike * (uGaussianP d1) - strike * (uGaussianP d2) ) )
        
        -- 2. estimate the swaption price for this particular genome,
        --    and also compute its likelihood
        (v0_mat,_,_) = bigV genome tmat0
        mux     = 0.0 - (bigmX genome today maturity today maturity)
        muy     = 0.0 - (bigmY genome today maturity today maturity)
        zc_mat  = zc maturity
        sqrt_bfun_a = sqrt $ bFun (2.0*a) tmat0
        sqrt_bfun_b = sqrt $ bFun (2.0*b) tmat0
        rhoxy   = rho * (bFun (a+b) tmat0) / (sqrt_bfun_a * sqrt_bfun_b)
        sigmax  = sigma * sqrt_bfun_a
        sigmay  = nu    * sqrt_bfun_b
        rhoxyc  = 1.0 - rhoxy * rhoxy
        rhoxycs = sqrt rhoxyc
        sigmay_rhoxycs = sigmay * rhoxycs
        t4      = (rhoxy * sigmay) / sigmax
        
        (cs, t1_cs, bas, bbs, aics, log_aics, scales) = unzip7 $
            map (\ i -> 
                    let beg_date = addMonths maturity (sw_freq*(fromIntegral i))
                        end_date = addMonths beg_date sw_freq
                        res      = (dateAct365 end_date beg_date) * strike

                        cii = if i == n_schedi-1 then 1.0 + res else res

                        date_tod1    = dateAct365 end_date today
                        (v0_end,_,_) = bigV genome date_tod1

                        date_tod2        = dateAct365 end_date maturity
                        (vt_end,bai,bbi) = bigV genome date_tod2
                        
                        expo_aici = 0.5 * (vt_end - v0_end + v0_mat)
                        fact_aici = cii * zc(end_date) / zc_mat

                        t1_c = bbi * (mux * t4 - (muy - 0.5*rhoxyc*sigmay*sigmay*bbi) ) + expo_aici
                        -- sanity = ! ( isinf(aici[i]) || isnan(aici[i]) );
                    in  (   fact_aici, t1_c, bai, bbi, fact_aici * (exp expo_aici), 
                            expo_aici + (log fact_aici), 0.0 - (bai + bbi * t4) )
                ) [0..n_schedi-1]
        
        exactYhatSpec = exactYhat   n_schedi b sigmax sigmay 
                                    rhoxy rhoxyc rhoxycs mux muy
                                    bas bbs aics log_aics
        eps = 0.5 * sigmax;
        f   = exactYhatSpec  mux
        g   = exactYhatSpec (mux + eps)
        h   = exactYhatSpec (mux - eps)
        df  = 0.5 * ( g - h ) / eps;

        sqrt2sigmax    = (sqrt 2.0) * sigmax
        t2             = rhoxy / (sigmax*rhoxycs)

        accums = 
            zipWith (\x_quad w_quad ->
                        let x      = sqrt2sigmax * x_quad + mux
                            yhat_x = f + df*(x - mux)
                            h1     = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux )

                            accs1= zipWith4 (\bbi t1_csi ci scalei ->
                                                let h2  = h1 + bbi * sigmay_rhoxycs
                                                    expo_aici = t1_csi + scalei*x
                                                    fact_aici = ci
                                                    --fact_aici * (exp expo_aici) * (uGaussianP (-h2))
                                                    expo_part = uGaussianPwithExpFactor (-h2) expo_aici
                                                in  fact_aici * expo_part
                                            ) bbs t1_cs cs scales
                            accum= reduce (+) 0.0 accs1
                               
                            tmp  = (sqrt 2.0) * x_quad; --(x - mux) / sigmax;
                            t1   = exp $ 0.0 - 0.5*tmp*tmp

                        in  w_quad * t1 * ( (uGaussianP (-h1)) - accum )   

                    ) hermite_coeffs hermite_weights 

        accum  = reduce (+) 0.0 accums

        swap_price = zc_mat * ( accum / (sqrt pi) )
        logLik     = logLikelihood black_price swap_price
    in  (black_price, swap_price, logLik)
--        tmp       = (swap_price - black_price) / black_price 
--        in  (black_price, swap_price, tmp * tmp)
    

--kkk :: Double
--kkk = -3.71901648545568  

exactYhat ::   Int ->  Double  ->  Double  ->  Double  ->  Double 
                   ->  Double  ->  Double  ->  Double  ->  Double 
                   -> [Double] -> [Double] -> [Double] -> [Double] 
                   ->  Double  ->  Double
exactYhat n_schedi b sigmax sigmay rhoxy rhoxyc rhoxycs mux muy 
            bai bbi aici log_aici x = 
    let kkk = -3.71901648545568
        (scales, los) = unzip $ zipWith4
                        (\ baii bbii aicii log_aicii -> 
                                let baix    = baii * x
                                    up_term = aicii * (exp (-baix))
                                    lo_term = ( log_aicii - baix ) / bbii
                                in  (up_term, lo_term)
                        ) bai bbi aici log_aici
        (up0,lo) = reduce   (\(u,l) (cu,cl) -> (u+cu, max l cl) ) 
                            (0.0, -infty) (zip scales los)

    in  if n_schedi == 1 then lo
        else 
            let log_s = log up0
                tmp   = log_s / (last bbi) -- [n_schedi-1];
                up    = if   tmp <= 0.0
                        then tmp
                        else let tmp2 = log_s / (head bbi)
                             in  if 0.0 <= tmp2
                                 then tmp2
                                 else (-infty)
                yl = lo - eps
                yu = up + eps

                y0 = sigmay * ( rhoxy * (x-mux) / sigmax + kkk * rhoxycs ) - rhoxyc/b + muy
                y1 = sigmay * ( rhoxy * (x-mux) / sigmax - kkk * rhoxycs )            + muy
            in  if      ( y1 <= yl ) then y1 + 1.0
                else if ( yu <= y0 ) then y0 - 1.0
                else let root_lb = max yl y0
                         root_ub = min yu y1
                         (root,_,err) = rootFindBrent scales bbi root_lb root_ub 1.0e-4 
                     in  if      err == (-infty) then y0 - 1.0
                         else if err ==   infty  then y1 + 1.0
                         else root


testSwaptionPricer :: [Double] -> [Double] -> Int -> Int
testSwaptionPricer hermite_coeffs hermite_weights arg0 = 
    let genome    = (0.02453, 0.98376, -0.82400, 0.11830, 0.02398)

        swaption1 = (10.0, 6.0, 4.0, 0.2454)
        (_,price1,_) = pricerOfSwaption hermite_coeffs hermite_weights genome swaption1
        price1' = trace ("# Pricer_of_swaption test I = 657.8215929143189 => ") $ 
                        1.0e4 * price1
        arg1 = if equalEps price1' 657.8215929143189 -- 657.82158867845
               then trace ("... SUCCESS !") arg0
               else trace ("... FAILS   !"++show price1') arg0

        swaption2 = (30.0, 6.0, 30.0, 0.16355)
        (_,price2,_) = pricerOfSwaption hermite_coeffs hermite_weights genome swaption2
        price2' = trace ("# Pricer_of_swaption test II = 1902.97628191498 => ") $ 
                        1.0e4 * price2
        arg2 = if equalEps price2' 1902.97628191498
               then trace ("... SUCCESS !") arg1
               else trace ("... FAILS   !") arg1

        swaption3 = (30.0, 6.0, 25.0, 0.1686)
        (_,price3,_) = pricerOfSwaption hermite_coeffs hermite_weights genome swaption3
        price3' = trace ("# Pricer_of_swaption test III = 1840.859126408099 => ") $ 
                        1.0e4 * price3
        arg3 = if equalEps price3' 1840.859126408099
               then trace ("... SUCCESS !") arg2
               else trace ("... FAILS   !") arg2
    in  arg3

