module Main where

import Control.Applicative
import Control.Monad
import Data.Maybe

------------------------------------
--- Requires the installation of ---
--- Parsec package, i.e.,        ---
---     cabal install parsec     ---
------------------------------------

import Text.Parsec hiding (optional, (<|>), token)
import Text.Parsec.String

import Data.Bits
import Data.List hiding (tail)
import Prelude hiding (tail)

import Payoff3Cond
import Debug.Trace

------------------------------
--- Parser-related Helpers ---
------------------------------
signed :: Parser String -> Parser String
signed p = (char '-' >> (('-':) <$> p)) <|> p

whitespace :: Parser ()
whitespace = spaces <* optional (string "//" >> manyTill anyChar (char '\n') >> whitespace)

lexeme :: Parser a -> Parser a
lexeme p = p <* whitespace

token :: String -> Parser ()
token = void . lexeme . string

readInt :: Parser Int
readInt = lexeme $ read <$> signed (many1 digit)

readDouble :: Parser Double
readDouble = lexeme $ read <$> signed (s2 <|> s1)
  where s1 = do bef <- many1 digit
                aft <- fromMaybe "" <$> optional ((:) <$> char '.' <*> many1 digit)
                return $ bef ++ aft
        s2 = (++) <$> (char '.' >> pure "0.") <*> many1 digit

readArray :: Parser a -> Parser [a]
readArray p = lexeme $ between (token "[") (token "]") (p `sepBy` token ",")

---------------------------------------------
--- Definition of lish-homomorphic reduce ---
---    requires a binary associative op   ---
---------------------------------------------

reduce :: (a -> a -> a) -> a -> [a] -> a
reduce bop ne lst = foldl bop ne lst

--zip5 :: [a] -> [b] -> [c] -> [d] -> [e] -> [(a,b,c,d,e)]
--zip5 [] [] [] [] [] = []
--zip5 (a:as) (b:bs) (c:cs) (d:ds) (e:es) = 
--    (a,b,c,d,e) :: (zip5 as bs cs ds es)
--zip5 as bs cs ds es = []

---------------------------------------------
--- Chunking Helpers                      ---
---------------------------------------------
chunk :: Int -> [a] -> [[a]]
chunk _ []  = []
chunk n arr = take n arr : chunk n (drop n arr)

---------------------------------------------
--- Sobol Quasi-random vector generation  ---
---------------------------------------------
grayCode :: Int -> Int
grayCode x = (x `shiftR` 1) `xor` x

xorInds :: Int -> Int -> [Int] -> Int
xorInds bits_num n dir_vs = foldl xor 0 $ map (dir_vs!!) is
  where bits = [0..bits_num-1]
        is = filter (testBit $ grayCode n) bits

sobolIndI :: Int -> [[Int]] -> Int -> [Int]
sobolIndI bits_num dir_vs n =
  map (xorInds bits_num n) dir_vs

sobolIndR :: Int -> [[Int]] -> Int -> [Double]
sobolIndR bits_num dir_vs n =
  map ((/divisor) . fromIntegral) arri
  where arri = sobolIndI bits_num dir_vs n
        divisor = 2.0 ** fromIntegral bits_num


------------------------------------------
--- To-Gaussian Distribution Transform.---
------------------------------------------
polyAppr :: Double -> Double -> Double -> Double -> Double -> Double
         -> Double -> Double -> Double -> Double -> Double -> Double
         -> Double -> Double -> Double -> Double -> Double -> Double
polyAppr x a0 a1 a2 a3 a4 a5 a6 a7 b0 b1 b2 b3 b4 b5 b6 b7 =
  (x*(x*(x*(x*(x*(x*(x*a7+a6)+a5)+a4)+a3)+a2)+a1)+a0) /
  (x*(x*(x*(x*(x*(x*(x*b7+b6)+b5)+b4)+b3)+b2)+b1)+b0)

smallcase :: Double -> Double
smallcase q =
  q *
  polyAppr
  ( 0.180625 - q * q)
  3.387132872796366608
  133.14166789178437745
  1971.5909503065514427
  13731.693765509461125
  45921.953931549871457
  67265.770927008700853
  33430.575583588128105
  2509.0809287301226727

  1.0
  42.313330701600911252
  687.1870074920579083
  5394.1960214247511077
  21213.794301586595867
  39307.89580009271061
  28729.085735721942674
  5226.495278852854561

intermediate :: Double -> Double
intermediate r =
  polyAppr (r - 1.6)
  1.42343711074968357734
  4.6303378461565452959
  5.7694972214606914055
  3.64784832476320460504
  1.27045825245236838258
  0.24178072517745061177
  0.0227238449892691845833
  7.7454501427834140764e-4

  1.0
  2.05319162663775882187
  1.6763848301838038494
  0.68976733498510000455
  0.14810397642748007459
  0.0151986665636164571966
  5.475938084995344946e-4
  1.05075007164441684324e-9

tail :: Double -> Double
tail r =
  polyAppr (r - 5.0)
  6.6579046435011037772
  5.4637849111641143699
  1.7848265399172913358
  0.29656057182850489123
  0.026532189526576123093
  0.0012426609473880784386
  2.71155556874348757815e-5
  2.01033439929228813265e-7

  1.0
  0.59983220655588793769
  0.13692988092273580531
  0.0148753612908506148525
  7.868691311456132591e-4
  1.8463183175100546818e-5
  1.4215117583164458887e-7
  2.04426310338993978564e-5

ugaussianEl :: Double -> Double
ugaussianEl p =
  if ((dp < 0.0) && (0.0 - dp <= 0.425)) || ((0.0 <= dp) && (dp <= 0.425))
  then smallcase dp
  else let pp = if dp < 0.0 then dp + 0.5 else 0.5 - dp
           r  = sqrt $ -log pp
           x = if r <= 5.0 then intermediate r else tail r
       in if dp < 0.0 then 0.0 - x else x
  where dp = p - 0.5

ugaussian :: [Double] -> [Double]
ugaussian = map ugaussianEl

------------------------------------------
--- Brownian-Bridge Implementation     ---
------------------------------------------
brownianBridgeDates :: Int -> [[Int]] -> [[Double]] -> [Double] -> [Double]
brownianBridgeDates num_dates bb_inds bb_data gauss =
  let bbrow = update (replicate num_dates 0.0)
              (bi!!0) (sd!!0 * gauss!!0)
      bbrow' = foldl f bbrow $ drop 1 $
               zip3 (zip3 bi li ri) (zip3 sd lw rw) gauss
  in zipWith (-) bbrow' (0.0:bbrow')
  where [bi,li,ri] = map (map $ subtract 1) bb_inds
        [sd,lw,rw] = bb_data
        f bbrow ((l,j,k),(x,y,z),zi) =
          let wk = bbrow !! k
              tmp = z * wk + x * zi
          in update bbrow l $ if (j + 1) == 0
                              then tmp
                              else tmp + y * (bbrow!!j)
        update a i v = take i a ++ [v] ++ drop (i+1) a

brownianBridge :: Int -> Int -> [[Int]] -> [[Double]] -> [Double] -> [[Double]]
brownianBridge num_paths num_dates bb_inds bb_data gaussian_arr =
  transpose $ map (brownianBridgeDates num_dates bb_inds bb_data) gauss2dT
  where gauss2d = chunk num_paths gaussian_arr
        gauss2dT = transpose gauss2d

------------------------------------------
--- Black-Scholes Implementation       ---
------------------------------------------
fftmp :: Int -> [[Double]] -> [Double] -> [Double]
fftmp num_paths md_c zi = map f [0..num_paths-1]
  where f j = sum $ zipWith (*) (take (j+1) zi) (take (j+1) $ md_c!!j)

correlateDeltas :: Int -> [[Double]] -> [[Double]] -> [[Double]]
correlateDeltas num_paths md_c = map $ fftmp num_paths md_c

combineVs :: [Double] -> [Double] -> [Double] -> [Double]
combineVs n_row vol_row dr_row =
  zipWith (+) dr_row $ zipWith (*) n_row vol_row

mkPrices :: [Double] -> [[Double]] -> [[Double]] -> [[Double]] -> [[Double]]
mkPrices md_starts md_vols md_drifts noises =
  drop 1 $ scanl (zipWith (*)) md_starts e_rows
  where e_rows = map (map exp) $ zipWith3 combineVs noises md_vols md_drifts

blackScholes :: Int -> [[Double]] -> [[Double]] -> [[Double]]
             -> [Double] -> [[Double]] -> [[Double]]
blackScholes num_paths md_c md_vols md_drifts md_starts bb_arr =
  mkPrices md_starts md_vols md_drifts noises
  where noises = correlateDeltas num_paths md_c bb_arr


------------------------------------------
--- Payoff functions                   ---
------------------------------------------
payoff :: Int -> [Double] -> [Double] -> [[Double]] -> Double
payoff contract md_disc md_detval xss = 
    if      contract == 1 then payoff1 md_disc md_detval xss
    else if contract == 2 then payoff2 md_disc xss
    else if contract == 3 then payoff3 md_disc xss
    else 0.0                

payoff1 :: [Double] -> [Double] -> [[Double]] -> Double
payoff1 md_disc md_detvals xss = 
    let detval = head md_detvals
        amount = ( (head (head xss)) - 4000.0 ) * detval
        amount'= max 0.0 amount
    in  trajInner amount' 0 md_disc

payoff2 :: [Double] -> [[Double]] -> Double
payoff2 md_disc xss
  | 1.0 <= mins !! 0 = trajInner 1150.0 0 md_disc
  | 1.0 <= mins !! 1 = trajInner 1300.0 1 md_disc
  | 1.0 <= mins !! 2 = trajInner 1450.0 2 md_disc
  | 1.0 <= mins !! 3 = trajInner 1600.0 3 md_disc
  | 1.0 <= mins !! 4 = trajInner 1750.0 4 md_disc
  | 0.75 < mins !! 4 = trajInner 1000.0 4 md_disc
  | otherwise = trajInner (1000 * (mins!!4)) 4 md_disc
  where divs = [ 1.0/3758.05, 1.0/11840.0, 1.0/1200.0 ]
        xss_div = map (zipWith (*) divs) xss
        mins = map minimum xss_div


payoff3 :: [Double] -> [[Double]] -> Double
payoff3 md_disc xss =
    let underlyings (i,j) = (xss !! i) !! j

        goto_50 = payoff3Cond xss
        x3309   = goto_50 || ( underlyings(366,1) <= 8288.0 )
        goto_40 = x3309 && 
                  ( (underlyings(366,0) < 3758.05) || 
                    (underlyings(366,2) < 1200.0 ) ||
                    (underlyings(366,1) < 11840.0)  )

        price1  = trajInner 100.0 0 md_disc

        price2  = if goto_40
                  then let m1 = min ((underlyings(366,2) / 1200.0 ) - 1.0) 
                                    ((underlyings(366,0) / 3758.05) - 1.0)
                           m2 = min ((underlyings(366,1) / 11840.0) - 1.0) m1
                           amount = (1000.0 * (1.0 + m2 ) ) 
                       in  trajInner amount 1 md_disc
                  else     trajInner 1000.0 1 md_disc
    in  price1 + price2
        


trajInner :: Double -> Int -> [Double] -> Double
trajInner amount i disc = amount * disc !! i


-------------------------------------------
--- Main Entry Point: price computation ---
-------------------------------------------
compute :: Int -> Int -> Int -> Int -> Int -> Int -> [[Int]] 
        -> [[[Double]]] -> [[[Double]]] -> [[[Double]]] 
        ->  [[Double]]  ->  [[Double]]  ->  [[Double]]
        -> [[Int]] -> [[Double]] 
        -> [Double]
compute contract  num_mc_it num_dates num_under num_models num_bits 
        dir_vs    md_cs     md_vols   md_drifts md_starts  md_detvals 
        md_discts bb_inds   bb_data =

  let num_det_pricers = length $ head md_detvals
      num_cash_flows  = length $ head md_discts

      sobol_mat = map    (sobolIndR num_bits dir_vs) [1..num_mc_it]
      gauss_mat = map     ugaussian sobol_mat
      bb_mat    = map    (brownianBridge num_under num_dates bb_inds bb_data) gauss_mat

      md_cvdsmat= zip5 md_cs md_vols md_drifts md_starts (replicate num_models bb_mat)
      bs_mats   = map    (\cvdsmat -> let (c, vol, drift, start, bb_mat) = cvdsmat
                                      in  map (blackScholes num_under c vol drift start) bb_mat 
                         )
                         md_cvdsmat 

--      bs_mats'  = trace (show bs_mats ++ show md_discts ++ show md_detvals ++ show contract) bs_mats 
--      bs_mats'  = trace (show bs_mats) bs_mats

      md_discmat= zip3 md_discts md_detvals bs_mats
      payoffs   = map    (\disc_bs -> let (disc, detvals, bsmat) = disc_bs
                                      in  map (payoff contract disc detvals) bsmat
                         )
                         md_discmat

--      payoffs'  = trace (show payoffs) payoffs 

      sum_prices= map (\payoff -> reduce (+) 0.0 payoff) payoffs
      prices    = map (\sp -> sp / fromIntegral num_mc_it) sum_prices
  in  prices

----------------------------------------------
--- Formatting the Output of the Benchmark ---
----------------------------------------------
validate :: [Double] -> [Double] -> (Int,Int,Int,Int) -> [String]
validate prices_ref prices info=
    let errs = map abs $ zipWith (-) prices_ref prices
        err  = reduce max 0.0 errs
        (contract,num_mc_it, num_dates, num_under) = info
    in ["// Generic Pricing Haskell Benchmark (List-Homomorphism Style):",
        "// Contract: " ++ show contract ++ ", MC Its#: " ++ show num_mc_it ++
            ", #Underlyings: " ++ show num_under ++ ", #Path Dates: " ++ 
            show num_dates,
        ( if err <= 0.0005 
          then "1\t\t// VALID   Result,"
          else "0\t\t// INVALID Result," ), 
        "0\t\t// Runtime in microseconds,",
        "1\t\t// CPU Threads,",
        show prices ++ "\t// Generic Pricing Result."]
        

-----------------------------------------
--- Entry point for Generic Pricing   ---
--- The only place where using Monads,---
--- e.g., parsing Dataset from StdIn  ---
-----------------------------------------
main :: IO ()
main = do s <- getContents 
          case parse run "input" s of
            Left  e        -> error $ show e
            Right (p,pr,i) -> do let msg_lst = validate pr p i
                                 putStrLn (msg_lst !! 0)  
                                 putStrLn (msg_lst !! 1)
                                 putStrLn (msg_lst !! 2)
                                 putStrLn (msg_lst !! 3)
                                 putStrLn (msg_lst !! 4)
                                 putStrLn (msg_lst !! 5)
                                 -- print $ validate pr p i 
--          either (error . show) print $ parse run "input" s
  where run = do contract <- readInt
                 num_mc_it<- readInt
                 num_dates<- readInt
                 num_under<- readInt

                 v <- compute contract num_mc_it num_dates num_under <$>
                        readInt <*> readInt <*> readInt2d <*>
                        readDouble3d <*> readDouble3d <*> readDouble3d <*>
                        readDouble2d <*> readDouble2d <*> readDouble2d <*>
                        readInt2d    <*> readDouble2d

                 r <- readDouble1d
                 return (v, r, (contract,num_mc_it, num_dates, num_under))

        readInt2d = readArray $ readArray readInt
        readDouble1d = readArray readDouble
        readDouble2d = readArray $ readArray readDouble
        readDouble3d = readArray $ readArray $ readArray readDouble

-- ghc -O2 -msse2 -rtsopts  PricingLexiFi.hs
-- ./PricingLexiFi +RTS -K128m -RTS < ../Data/Medium/input.data
