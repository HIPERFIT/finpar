module Main where

import Control.DeepSeq
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
import Data.List 
import Prelude 

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

---------------------------------------------
---------------------------------------------

initGrid :: Double -> Double -> Double -> Double
        ->  Int    -> Int    -> Int 
        -> ( Int, Int, [Double], [Double], [Double], 
             [[Double]], [[Double]], [[Double]], [[Double]] )
initGrid s0 alpha nu t num_x num_y num_t = 
    let logAlpha = log alpha
        myTimeline = map  (\i -> t * (fromIntegral i) / ((fromIntegral num_t) - 1.0)) [0..num_t-1]
        stdX = 20.0 * alpha * s0 * (sqrt t)
        stdY = 10.0 * nu         * (sqrt t)
        (dx, dy) = (stdX / fromIntegral num_x, stdY / fromIntegral num_y)
        (myXindex, myYindex) = (truncate (s0 / dx), num_y `div` 2)
        myX = map (\i -> (fromIntegral i) * dx - (fromIntegral myXindex) * dx + s0      ) [0..num_x-1]
        myY = map (\i -> (fromIntegral i) * dy - (fromIntegral myYindex) * dy + logAlpha) [0..num_y-1]
        xXy = replicate num_x (replicate num_y 0.0)
        (myMuX, myVarX, myMuY, myVarY) = (xXy, xXy, xXy, xXy)
    in  (myXindex, myYindex, myX, myY, myTimeline, myMuX, myVarX, myMuY, myVarY)

---------------------------------------------
---------------------------------------------

---------------------------------------------
---------------------------------------------

initOperator :: [Double] 
            -> ( [(Double,Double,Double)], [(Double,Double,Double)] )
initOperator x = 
    let n      = length x
        mid_x  = zip3 x (tail x) (tail (tail x))

        dxu    = (x !! 1) - (x !! 0)
        dxl    = 0.0
        dxlow  = (0.0, -1.0 / dxu, 1.0 / dxu)
        dxxlow = (0.0, 0.0, 0.0)

        --- Implements a parallel loop such as:   ---
        --- doall i = 2, N-1, 1                   ---
        ---   (a,b,c)  = f(x[i-1], x[i], x[i+1])  ---
        ---   res[i-1] = (a,b,c)                  ---
        dxmid  = map (\(xim1, xi, xip1) -> 
                        let (dxl, dxu) = (xi - xim1, xip1 - xi)
                            dxlpxu     = dxl + dxu
                        in ( (-dxu/dxl )          / dxlpxu,
                             ( dxu/dxl - dxl/dxu) / dxlpxu,
                             ( dxl/dxu )          / dxlpxu
                           ) 
                     ) mid_x 
        dxxmid = map (\(xim1, xi, xip1) -> 
                        let (dxl, dxu) = (xi - xim1, xip1 - xi)
                            dxlpxu     = dxl + dxu
                        in ( ( 2.0/dxl )            / dxlpxu,
                             -2.0*(1.0/dxl+1.0/dxu) / dxlpxu,
                             ( 2.0/dxu )            / dxlpxu
                           ) 
                     ) mid_x
        dxll   = (x !! (n-1)) - (x !! (n-2))
        dxul   = 0.0
        dxhigh = [(-1.0 / dxll, 1.0 / dxll, 0.0)]
        dxxhigh= [(0.0, 0.0, 0.0)]
    in ( dxlow : (dxmid++dxhigh), dxxlow : (dxxmid++dxxhigh) ) 

-----------------------------------------------
--- A tail-recursive function modeling a    ---
--- do loop that goes from i = 0 to bound-1 ---
-----------------------------------------------

doLoop :: Int -> Int
       -> ( [Double], [(Double,Double,Double)], [(Double,Double,Double)], 
            [Double], [(Double,Double,Double)], [(Double,Double,Double)], 
            [Double], Double, Double, Double )
       -> ([[Double]], [[Double]], [[Double]], [[Double]], [[Double]])
       -> ([[Double]], [[Double]], [[Double]], [[Double]], [[Double]])
doLoop i bound loop_ros loop_variants  = 
    if  i == bound 
    then loop_variants
    else let (myResult, myMuX, myVarX, myMuY, myVarY) = loop_variants

             ( myX, myDx, myDxx, myY, myDy, myDyy,
               myTimeline, alpha, beta, nu )          = loop_ros

             j = (bound-1-i)

             (myMuX', myVarX', myMuY', myVarY')       =
                    updateParams myX myY myTimeline j alpha beta nu

             

             myResult'= rollback j  myX  myY   myTimeline  myResult
                                 myMuX'  myDx  myDxx       myVarX'
                                 myMuY'  myDy  myDyy       myVarY' 

             loop_variants' = (myResult', myMuX', myVarX', myMuY', myVarY')

         in  -- Hack to avoid space leak
          myResult' `deepseq` myMuX' `deepseq` myVarX' `deepseq` myMuY' `deepseq` myVarY' `deepseq`
          doLoop (i+1) bound loop_ros loop_variants'

---------------------------------------------
---------------------------------------------

updateParams :: [Double] -> [Double] -> [Double] -> Int -> Double -> Double -> Double 
             -> ( [[Double]], [[Double]], [[Double]], [[Double]] )
updateParams myX myY myTimeline g alpha beta nu =
    unzip4 $ map (\ xi -> unzip4 $
                          map (\ yj -> let b = beta * log(xi) + yj
                                           c = 0.5 * nu * nu * (myTimeline !! g)
                                       in  ( 0.0, exp (2.0 * (b - c) ), 0.0, nu * nu )
                              ) myY
                 ) myX
---------------------------------------------
---  Helpers for the rollback function    --- 
---    ToDo: make them nice!!!            ---
---------------------------------------------

-- explicitX = explicitXY  dt_inv 0.5  (transpose my_result)  (transpose myMuX)  (transpose myVarX)  myDx  myDxx
-- explicitY = explicitXY  0.0    1.0  my_result  myMuY  myVarY  myDy  myDyy
explicitXY ::   Double   ->   Double
           -> [[Double]] -> [[Double]] -> [[Double]]
           -> [ (Double,Double,Double) ] 
           -> [ (Double,Double,Double) ]
           -> [[Double]]
explicitXY  dt_inv  fact  my_result  myMuY  myVarY  myDy  myDyy =
        map (\ tup_i -> 
                let (res_row, myMuY_row, myVarY_row) = tup_i
                    res_row3 = zip3 (0.0:res_row) res_row ((tail res_row) ++ [0.0])

                in  map (\ ((myDy_i, myDyy_i, myMuY_ij, myVarY_ij),(r_jm1, r_j, r_jp1)) ->
                            let ((dy0, dy1, dy2), (dyy0, dyy1, dyy2)) = (myDy_i, myDyy_i)
                                res0 =  dt_inv * r_j
                                res1 =  fact * (myMuY_ij*dy0 + 0.5*myVarY_ij*dyy0) * r_jm1
                                res2 =  fact * (myMuY_ij*dy2 + 0.5*myVarY_ij*dyy2) * r_jp1
                                res3 =  fact * (myMuY_ij*dy1 + 0.5*myVarY_ij*dyy1) * r_j 
                            in  res0 + res1 + res2 + res3

                        ) (zip (zip4 myDy myDyy myMuY_row myVarY_row) res_row3)

            ) (zip3 my_result myMuY myVarY)

------------------------------------------
--- loop-like definition for explicitX ---
---  inneficient due to list traversal ---
---  to find every index.              ---
--- Same is possible for explicitY0    ---
------------------------------------------
explicitX0 :: Int -> Int -> Double -> [[Double]] -> [[Double]] 
           -> [ (Double,Double,Double) ] -> [ (Double,Double,Double) ] -> [[Double]]
           -> [[Double]]
explicitX0    num_x  num_y  dt_inv  my_result  
              myMuX  myDx   myDxx   myVarX    = 
    map (\ j -> 
                map (\ i -> let res0 = dt_inv * ((my_result !! i) !! j) -- my_result[i,j]
                                ((dx0, dx1, dx2), (dxx0, dxx1, dxx2)) = (myDx !! i, myDxx !! i)

                                myMuX_ij  = ( myMuX !! i) !! j --  myMuX[i,j]
                                myVarX_ij = (myVarX !! i) !! j -- myVarX[i,j]

                                res1 =  if i == 0 then 0.0
                                        else 0.5 * (myMuX_ij*dx0 + 0.5*myVarX_ij*dxx0) * 
                                             ((my_result !! (i-1)) !! j) -- my_result[i-1, j] 
                                res2 =  if i == num_x - 1 then 0.0
                                        else 0.5 * (myMuX_ij*dx2 + 0.5*myVarX_ij*dxx2) * 
                                             ((my_result !! (i+1)) !! j) -- my_result[i+1,j]

                                res3 =       0.5 * (myMuX_ij*dx1 + 0.5*myVarX_ij*dxx1) * 
                                             ((my_result !! i) !! j)     -- my_result[i  ,j]

                            in  res0 + res1 + res2 + res3 
                                       

                    ) [0..num_x-1]

        ) [0..num_y-1]

tridag :: [Double] -> [Double] -> [Double] -> [Double] 
       -> ( [Double], [Double] )
tridag a b c r = 
    let u0  = head r
        uu0 = head b
 
        -- scanl's binary operator is NOT associative
        uu  = scanl (\ uuim1 (ai, bi, cim1) -> 
                        let beta = ai / uuim1
                        in  bi - beta*cim1
                    ) uu0 (zip3 (tail a) (tail b) c)
 
        -- scanl's binary operator is NOT associative
        u   = scanl (\ uim1 (ai, ri, uuim1) ->
                        let beta = ai / uuim1
                        in  ri - beta*uim1
                    ) u0  (zip3 (tail a) (tail r) uu)

        ur  = reverse u
        uur = reverse uu
        ur0'= (head ur) / (head uur)
        ur' = scanl (\ uip1 (uri, uuri, cri) -> 
                            (uri - cri*uip1) / uuri
                    ) ur0' (zip3 (tail ur) (tail uur) (tail $ reverse c))
    in  (reverse ur', uu)


tridagPar :: [Double] -> [Double] -> [Double] -> [Double] 
          -> ( [Double], [Double] )
tridagPar a b c r = 
    let -- u0  = head r
        -- uu0 = head b
        u0  = head b
        uu0 = head r 

        -- creating the 2x2 matrices 
        mats    = map (\(ai,bi,cim1)-> (bi, -ai*cim1, 1.0, 0.0)) (zip3 (tail a) (tail b) c)
        -- scan with 2x2 matrix multiplication
        scanmat = scanl (\(x1,y1,z1,w1) (x2,y2,z2,w2) -> 
                            let dv = 1.0/(x1*x2)
                            in  ( (x1*x2+y1*z2)*dv, 
                                  (x1*y2+y1*w2)*dv,
                                  (z1*x2+w1*z2)*dv,
                                  (z1*y2+w1*w2)*dv )
                        ) (1.0, 0.0, 0.0, 1.0) mats

        -- compute the first recurrence result
        uu = map (\(x,y,z,w) -> (x*u0 + y) / (z*u0 + w) ) scanmat

        -- compute the second recurrence
        pairs = map (\(ai,ri,uuim1)->(ri, -ai/uuim1)) (zip3 (tail a) (tail r) uu) 

        scanpairs = scanl (\ (x1,y1) (x2,y2) -> (x2+y2*x1, y1*y2) ) (0.0,1.0) pairs 

        u = map (\(x,y) -> x + u0*y) scanpairs

        -- backwards recurrence
        ur  = reverse u
        uur = reverse uu
        ur0'= (head ur) / (head uur)

        pairsr = map (\(uri, uuri, cri)->(uri/uuri, -cri/uuri)) (zip3 (tail ur) (tail uur) (tail $ reverse c)) 
        scanpairsr = scanl (\ (x1,y1) (x2,y2) -> (x2+y2*x1, y1*y2) ) (0.0,1.0) pairsr
        ur' = map (\(x,y) -> x + ur0'*y) scanpairsr

    in  (reverse ur', uu)

---------------------------------------------
--- rollback: the brain of the program    ---
---------------------------------------------

rollback ::  Int -> [Double] -> [Double] -> [Double] -> [[Double]]
         -> [[Double]] -> [(Double,Double,Double)] -> [(Double,Double,Double)] -> [[Double]]
         -> [[Double]] -> [(Double,Double,Double)] -> [(Double,Double,Double)] -> [[Double]]
         -> [[Double]] 

rollback g  myX  myY  myTimeline  myResult
         myMuX  myDx  myDxx  myVarX
         myMuY  myDy  myDyy  myVarY = 

    let (numX, numY) = (length myX, length myY)
        numZ  = max numX numY
        dtInv = 1.0 / ( (myTimeline !! (g+1)) - (myTimeline !! g) )

        u0= explicitXY dtInv  0.5  (transpose myResult)  (transpose myMuX)  (transpose myVarX)  myDx  myDxx
        v0= explicitXY 0.0    1.0  myResult  myMuY  myVarY  myDy  myDyy
 
        u1= map (\ (us, vs) -> zipWith (+) us vs )
                (zip u0 (transpose v0))

        u2= map (\ t -> let (uj, myDx, myDxx, myMuX, myVarX) = t
                            (a,b,c) = unzip3 $
                                        map (\ tt -> let (myDx, myDxx, myMuX, myVarX) = tt
                                                         (dx0,  dx1,  dx2 ) = myDx  
                                                         (dxx0, dxx1, dxx2) = myDxx
                                                     in  (  -0.5*(myMuX*dx0 + 0.5*myVarX*dxx0),
                                                            dtInv - 0.5*(myMuX*dx1 + 0.5*myVarX*dxx1),
                                                            -0.5*(myMuX*dx2+0.5*myVarX*dxx2)          ) )
                                            (zip4 myDx myDxx myMuX myVarX)

                            (uj', yy) = tridag a b c uj  -- tridagPar a b c uj 
                        in  uj' )
                (zip5 u1 (replicate numY myDx) (replicate numY myDxx) (transpose myMuX) (transpose myVarX))

    in  map (\ t -> let (ui, vi, myDy, myDyy, myMuY, myVarY) = t
                        (a,b,c) = unzip3 $
                                    map (\ tt -> let (myDy, myDyy, myMuY, myVarY) = tt
                                                     (dy0,  dy1,  dy2 ) = myDy
                                                     (dyy0, dyy1, dyy2) = myDyy
                                                 in  (  -0.5*(myMuY*dy0+0.5*myVarY*dyy0),
                                                        dtInv - 0.5*(myMuY*dy1+0.5*myVarY*dyy1),
                                                        -0.5*(myMuY*dy2+0.5*myVarY*dyy2)        ) )
                                        (zip4 myDy myDyy myMuY myVarY)

                        y = map (\ (u,v) -> dtInv * u - 0.5 * v) (zip ui vi)

                        (ri, yy) = tridag a b c y -- tridagPar a b c y 
                    in  ri )
            (zip6 (transpose u2) v0 (replicate numX myDy) (replicate numX myDyy) myMuY myVarY)

---------------------------------------------
---------------------------------------------

-------------------------------------------------
--- Cosmin ToDo: make MyX & MyY & myTimeline  ---
---                  vectors instead of lists ---
-------------------------------------------------

value ::    (Double, Double, Double, Double, Double)
        ->  (Int,    Int,    Int,    Int) -> Double
        ->  Double
value  params it_spaces strike = 
    let (outer, num_x, num_y, num_t) = it_spaces
        (s0, t, alpha, nu, beta)     = params

        (myXindex, myYindex, myX, myY, myTimeline, myMuX, myVarX, myMuY, myVarY) =
            initGrid s0 alpha nu t num_x num_y num_t
        (myDx, myDxx) = initOperator myX
        (myDy, myDyy) = initOperator myY

        my_result_ini = map (\ xi -> replicate num_y (max (xi - strike) 0.0 )) myX
        
        loop_variants = (my_result_ini, myMuX, myVarX, myMuY, myVarY)
        loop_ros      = (myX, myDx, myDxx, myY, myDy, myDyy, myTimeline, alpha, beta, nu)

        (my_result,_,_,_,_) = doLoop 0 (num_t-1) loop_ros loop_variants 

    in  (my_result !! myXindex) !! myYindex  -- my_result[myXindex, myYindex]

------------------------------------------------
--- Main Entry Point: volatility calibration ---
------------------------------------------------
compute :: (Int,    Int,    Int,    Int) 
        -> (Double, Double, Double, Double, Double)
        -> [Double]
compute    it_spaces params =
    let (s0, t, alpha, nu, beta) = params
        (outer, num_x, num_y, num_t) = it_spaces
        strikes = map (\i -> 0.001 * (fromIntegral i)) [0..outer-1]
        res     = map (value params it_spaces) strikes 
    in  res

----------------------------------------------
--- Formatting the Output of the Benchmark ---
----------------------------------------------
validate :: [Double] -> [Double] -> (Int,Int,Int,Int) -> [String]
validate res_ref res info=
    let errs = map abs $ zipWith (-) res_ref res
        err  = reduce max 0.0 errs
        (outer, num_x, num_y, num_t) = info
    in ["// Generic Pricing Haskell Benchmark (List-Homomorphism Style):",
        "// OUTER: " ++ show outer ++ ", NUM_X: " ++ show num_x ++
         ", NUM_Y: " ++ show num_y ++ ", NUM_T: " ++ show num_t,
        ( if err <= 0.00001 -- 0.0005 
          then "1\t\t// VALID   Result,"
          else "0\t\t// INVALID Result," ), 
        "0\t\t// Runtime in microseconds,",
        "1\t\t// CPU Threads,",
        show res ++ "\t// Volatility Calibration Result."]
        

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
                                 _ <- mapM putStrLn msg_lst  
                                 putStrLn ""
  where run = do outer <- readInt
                 num_x <- readInt
                 num_y <- readInt
                 num_t <- readInt

                 let params = (0.03, 5.0, 0.2, 0.6, 0.5) --(s0, t, alpha, nu, beta)

                 let v = compute (outer, num_x, num_y, num_t) params
                 r <- readDouble1d
                 return (v, r, (outer, num_x, num_y, num_t))

        readInt2d    = readArray $ readArray readInt
        readDouble1d = readArray readDouble
        readDouble2d = readArray $ readArray readDouble
        readDouble3d = readArray $ readArray $ readArray readDouble

-- ghc -O2 -msse2 -rtsopts  PricingLexiFi.hs
-- ./PricingLexiFi +RTS -K128m -RTS < ../Data/Medium/input.data
