module Genome (
                g_mins
              , g_maxs
              , mutateDimsALL
              , mutateDimsONE
              , mcmcDE
              ) where

import Constants

import Data.Vector.Unboxed(Vector) 
import qualified Data.Vector.Unboxed as V

import Debug.Trace

g_mins :: Vector Double 
g_mins = V.fromList [ eps0,      eps0,      -1.0+eps0, eps0, eps0 ]

g_maxs :: Vector Double
g_maxs = V.fromList [ 1.0-eps0,  1.0-eps0,   1.0-eps0, 0.2,  0.2  ]

g_inis :: Vector Double
g_inis = V.fromList [ 0.02,      0.02,       0.0,      0.01, 0.04 ]
--                  [ g_a,       g_b,        g_rho,    g_nu, g_sigma ]

unif_ampl_ratio :: Double
unif_ampl_ratio = 0.005 / genome_scale



-- rand_vct[5], i.e., needs 5 random numbers
mutateDimsALL ::   Vector Double -> Vector Double -> Vector Double 
                -> (Vector Double, Double)
mutateDimsALL genome genome_prop rand_vct = -- 5 rand nums
    let (genome_res, bf_facts) = V.unzip $
            V.zipWith5  (mutateHelper unif_ampl_ratio)
                        rand_vct g_mins g_maxs genome genome_prop
        bf_res = V.foldl (*) 1.0 bf_facts
    in  (genome_res, bf_res)


-- needs one random number
mutateDimsONE ::   Int -> Vector Double -> Vector Double -> Vector Double 
                -> (Vector Double, Double)
mutateDimsONE dim_j genome genome_prop rand_vct = -- r01 1 rand num
    let (gene, gene_prop)    = (genome V.! dim_j, genome_prop V.! dim_j)
        (gene_min, gene_max) = (g_mins V.! dim_j, g_maxs      V.! dim_j)
        r01 = rand_vct V.! dim_j        

        (gene_new, bf_res) = mutateHelper unif_ampl_ratio r01 
                                          gene_min gene_max gene gene_prop
        genome_res = V.imap (\i g -> if i == dim_j then gene_new else g) genome
    in  (genome_res, bf_res)


-- rand_vct[8], i.e., needs 8 random numbers
mcmcDE ::   Int -> [Vector Double] -> Vector Double -> Int
         -> (Vector Double, Double)
mcmcDE pop_size genomes rand_vct j = -- 8 rand nums
    let gamma_avg  = 2.38 / (sqrt (2.0*(fromIntegral genome_dim)))
        ampl_ratio = 0.1  * unif_ampl_ratio
        
        k0 = truncate $ (fromIntegral (pop_size-1))*(rand_vct V.! 0)
        (k, cand_UB) = if ( k0 == j ) 
                       then (pop_size - 1, pop_size - 2)
                       else (k0,           pop_size - 1)

        l0 = truncate $ (fromIntegral cand_UB)*(rand_vct V.! 1)
        l = -- trace ("(k,j,l0,cand_UB) is: "++show k++" "++show j++" "++show l0++" "++show cand_UB) $
            if l0 == j || l0 == k
            then cand_UB 
            else l0 

        -- proposal gamma: integrated out from the adviced 
        --   Multivariate Gaussian with Gaussian target (Braak, 2006)
        gamma1 = gamma_avg - 0.5 + (rand_vct V.! 2)  -- * instead of + ?

        gene_rands = V.drop 3 rand_vct
        (genome_j, genome_k, genome_l) = (genomes !! j, genomes !! k, genomes !! l) 
        -- g_a     [j + POP_SIZE] = constrain_dim1( 0, perturbation( g_a    [j], g_a    [k], g_a    [l], 0, gamma1, ampl_ratio ) );
        genome_res = V.zipWith6 (perturbation gamma1 ampl_ratio)
                                gene_rands g_mins g_maxs genome_j genome_k genome_l
    in  (genome_res, 1.0)

--perturbation gamma1 amplitude_ratio r01 gene_min gene_max gene gene_k gene_l

----------------------------------------
----------------------------------------
---------- Helper Functions ------------
----------------------------------------
----------------------------------------

mutateHelper ::    Double -> Double -> Double -> Double -> Double -> Double 
                -> (Double, Double)
mutateHelper amplitude_ratio r01 gene_min gene_max gene gene_prop = 
    let amplitude     = abs $ (gene_max - gene_min) * amplitude_ratio
        semiamplitude = amplitude / 2.0

        tmp_min_max_f = min gene_max (gene + semiamplitude)
        tmp_max_min_f = max gene_min (gene - semiamplitude)
        forward_range = tmp_min_max_f - tmp_max_min_f

        bf_fact = if forward_range > 0.0 && amplitude > 0.0
                  then let tmp_min_max_b  = min gene_max (gene_prop + semiamplitude)
                           tmp_max_min_b  = max gene_min (gene_prop - semiamplitude)
                           backward_range = tmp_min_max_b - tmp_max_min_b
                       in  backward_range / forward_range
                  else 1.0
        diff = amplitude * r01 - semiamplitude
        new_gene = constrainDim gene_min gene_max (gene+diff)
    in  (new_gene, bf_fact)

constrainDim :: Double -> Double -> Double -> Double
constrainDim gene_min gene_max gene = 
    max gene_min $ min gene_max gene


perturbation ::    Double -> Double -> Double -> Double 
                -> Double -> Double -> Double -> Double 
                -> Double
perturbation gamma1 amplitude_ratio r01 gene_min gene_max gene gene_k gene_l = 
    let amplitude     = abs $ (gene_max - gene_min) * amplitude_ratio
        semiamplitude = amplitude / 2.0
        perturb       = amplitude * r01 - semiamplitude
        new_gene = constrainDim gene_min gene_max $ 
                                gene + perturb + gamma1 * ( gene_k - gene_l )
    in  new_gene

