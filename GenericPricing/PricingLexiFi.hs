module Main where

import Control.Applicative
import Control.Monad
import Data.Maybe

import Text.Parsec hiding (optional, (<|>), token)
import Text.Parsec.String

import Data.Bits
import Data.List hiding (tail)
import Prelude hiding (tail)

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
                aft <- fromMaybe "" <$> optional s2
                return $ bef ++ aft
        s2 = (++) <$> (char '.' >> pure "0.") <*> many1 digit

readArray :: Parser a -> Parser [a]
readArray p = lexeme $ between (token "[") (token "]") (p `sepBy` token ",")

chunk :: Int -> [a] -> [[a]]
chunk _ []  = []
chunk n arr = take n arr : chunk n (drop n arr)

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

brownianBridgeDates :: Int -> [[Int]] -> [[Double]] -> [Double] -> [Double]
brownianBridgeDates num_dates bb_inds bb_data gauss =
  let bbrow = update (replicate num_dates 0.0)
              (bi!!0) ((sd!!0) * (gauss!!0))
      bbrow' = foldl f bbrow (zip3 (zip3 bi li ri) (zip3 sd lw rw) gauss)
  in zipWith (-) bbrow' (0:bbrow')
  where [bi,li,ri] = map (drop 1 . map (subtract 1)) bb_inds
        [sd,lw,rw] = map (drop 1) bb_data
        f bbrow ((l,j,k),(x,y,z),zi) =
          let wk = bbrow !! k
              tmp = z * wk + x * zi
          in update bbrow l $ if j + 1 == 0
                              then tmp
                              else tmp + y * (bbrow!!j)
        update a i v = take i a ++ [v] ++ drop (i+1) a

brownianBridge :: Int -> Int -> [[Int]] -> [[Double]] -> [Double] -> [[Double]]
brownianBridge num_paths num_dates bb_inds bb_data gaussian_arr =
  transpose $ map (brownianBridgeDates num_dates bb_inds bb_data) gauss2dT
  where gauss2d = chunk num_paths gaussian_arr
        gauss2dT = transpose gauss2d

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

trajInner :: Double -> Int -> [Double] -> Double
trajInner amount i disc = amount * disc !! i

compute :: Int -> Int -> Int -> Int -> [[Int]]
        -> [[Double]] -> [[Double]] -> [[Double]]
        -> [Double] -> [Double] -> [Double]
        -> [[Int]] -> [[Double]] -> Double
compute num_mc_it num_dates num_und num_bits dir_vs
  md_c md_vols md_drifts md_st _ md_disc bb_inds bb_data =
  payoff / fromIntegral num_mc_it
  where sobol_mat = map (sobolIndR num_bits dir_vs) [1..num_mc_it]
        gauss_mat = map ugaussian sobol_mat
        bb_mat = map (brownianBridge num_und num_dates bb_inds bb_data) gauss_mat
        bs_mat = map (blackScholes num_und md_c md_vols md_drifts md_st) bb_mat
        payoff = sum $ map (payoff2 md_disc) bs_mat

main :: IO ()
main = do s <- getContents
          either (error . show) print $ parse run "input" s
  where run = compute <$>
              readInt <*> readInt <*> readInt <*> readInt <*> readInt2d <*>
              readDouble2d <*> readDouble2d <*> readDouble2d <*>
              readDouble1d <*> readDouble1d <*> readDouble1d <*>
              readInt2d <*> readDouble2d
        readInt2d = readArray $ readArray readInt
        readDouble1d = readArray readDouble
        readDouble2d = readArray $ readArray readDouble
