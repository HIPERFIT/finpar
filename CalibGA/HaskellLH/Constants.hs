module Constants 
        (
          eps
        , eps0
        , infty
--        , rrr
        , sobolNum
--        , pop_size
--        , num_mcmc_its
        , genome_dim
        , genome_scale
        , MOVE_TYPE(DIMS_ALL,DIMS_ONE,DEMCMC)
        , selectMoveType
        , equalEps
        , reduce
        , Swaption
        , Genome
        ) where


import Data.Bits
import Data.Vector.Unboxed(Vector) 
import qualified Data.Vector.Unboxed as V

import Debug.Trace

data MOVE_TYPE  = DIMS_ALL 
                | DIMS_ONE 
                | DEMCMC
                    deriving(Show)

type Swaption = (Double, Double, Double, Double)
type Genome   = (Double, Double, Double, Double, Double)

equalEps :: Double -> Double -> Bool
equalEps x1 x2 = abs (x1 - x2) <= 1.0e-8

selectMoveType :: Double -> MOVE_TYPE
selectMoveType move_selected = 
    if      move_selected <= 0.2 then DIMS_ALL
    else if move_selected <= 0.5 then DIMS_ONE
    else {- move_selected <= 1.0 -}   DEMCMC  


eps0 :: Double   
eps0 = 1.0e-3

eps :: Double
eps  = 1.0e-5

infty :: Double
infty = 1.0e49

--pop_size :: Int
--pop_size = 128
--
--num_mcmc_its :: Int
--num_mcmc_its = truncate $ (2.0 + (log (fromIntegral genome_dim))) * 
--                          (1.0 + (log genome_scale)) * 100.0

genome_scale :: Double
genome_scale = 1.0

genome_dim :: Int
genome_dim = 5

---------------------------------------------
--- Definition of lish-homomorphic reduce ---
---    requires a binary associative op   ---
---------------------------------------------

reduce :: (a -> a -> a) -> a -> [a] -> a
reduce bop ne lst = foldl bop ne lst

---------------------------------
--- Sobol Numbers Computation ---
---------------------------------

num_sobol_bits :: Int
num_sobol_bits = 30

sobol_norm :: Double
sobol_norm = 1.0 / (1.0 + (2.0 ** (fromIntegral num_sobol_bits)))

grayCode :: Int -> Int
grayCode x = (x `shiftR` 1) `xor` x

sobolNum  :: Vector Int -> Int -> Double
sobolNum dir_vct n = 
    let n_gray = grayCode n
        ires = V.ifoldl (\res i dir -> 
                            if testBit n_gray i
                            then res `xor` dir
                            else res
                        ) 0 dir_vct
        rres = sobol_norm * (fromIntegral ires)
    in  rres -- trace ("Sobol Num "++show n++" is: "++show rres) rres

--sobolNums :: Vector int -> Int -> Int -> Vector Double
--sobolNums dir_vct beg_offs end_offs = 
    

-- ifoldl :: Unbox b => (a -> Int -> b -> a) -> a -> Vector b -> a
