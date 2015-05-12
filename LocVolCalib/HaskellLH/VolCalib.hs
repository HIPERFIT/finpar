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

import Data.List 
import Prelude 

--import Debug.Trace

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
        -> ( Int, Int, [Double], [Double], [Double])
initGrid s0 alpha nu t num_x num_y num_t = 
    let logAlpha = log alpha
        myTimeline = map  (\i -> t * (fromIntegral i) / ((fromIntegral num_t) - 1.0)) [0..num_t-1]
        stdX = 20.0 * alpha * s0 * (sqrt t)
        stdY = 10.0 * nu         * (sqrt t)
        (dx, dy) = (stdX / fromIntegral num_x, stdY / fromIntegral num_y)
        (myXindex, myYindex) = (truncate (s0 / dx), num_y `div` 2)
        myX = map (\i -> (fromIntegral i) * dx - (fromIntegral myXindex) * dx + s0      ) [0..num_x-1]
        myY = map (\i -> (fromIntegral i) * dy - (fromIntegral myYindex) * dy + logAlpha) [0..num_y-1]
    in  (myXindex, myYindex, myX, myY, myTimeline)

---------------------------------------------
---------------------------------------------

---------------------------------------------
---------------------------------------------

initOperator :: [Double] 
            -> ( [(Double,Double,Double)], [(Double,Double,Double)] )
initOperator x = 
    let n      = length x
        mid_x  = zip3 x (tail x) (tail (tail x))

        dxu0   = (x !! 1) - (x !! 0)
        dxlow  = (0.0, -1.0 / dxu0, 1.0 / dxu0)
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
       -> [[Double]]
       -> [[Double]]
doLoop i bound loop_ros myResult  = 
    if  i == bound 
    then myResult
    else let ( myX, myDx, myDxx, myY, myDy, myDyy,
               myTimeline, alpha, beta, nu )          = loop_ros

             j = (bound-1-i)

             (myMuX, myVarX, myMuY, myVarY)       =
                    updateParams myX myY myTimeline j alpha beta nu

             

             myResult'= rollback j  myTimeline  myResult
                                 myMuX  myDx  myDxx  myVarX
                                 myMuY  myDy  myDyy  myVarY 

         in  -- Hack to avoid space leak
          myResult' `deepseq` myMuX `deepseq` myVarX `deepseq` myMuY `deepseq` myVarY `deepseq`
          doLoop (i+1) bound loop_ros myResult'

---------------------------------------------
---------------------------------------------

updateParams :: [Double] -> [Double] -> [Double] -> Int -> Double -> Double -> Double 
             -> ( [[Double]], [[Double]], [[Double]], [[Double]] )
updateParams myX myY myTimeline g alpha beta nu =
    let ( numX, numY ) = ( length myX, length myY )    
        myMuY  = replicate numX (replicate numY (0.0*alpha))
        myVarY = replicate numX (replicate numY (nu*nu)    )
        myMuX  = replicate numY (replicate numX 0.0        )
        myVarX = map (\ yj -> 
                        map (\ xi -> 
                                let b = beta * log(xi) + yj
                                    c = 0.5 * nu * nu * (myTimeline !! g)
                                in  exp (2.0 * (b - c))
			                ) myX
		             ) myY
    in  ( myMuX, myVarX, myMuY, myVarY )

------------------------------------------------------
--- TRIDAG: tridiagonal solver.                    ---
---         sequential and parallel implementation ---
------------------------------------------------------

tridagSeq :: [Double] -> [Double] -> [Double] -> [Double] 
          -> [Double]
tridagSeq a b c r = 
    let u0  = head r
        uu0 = head b
 
        -- sequential scanl: binary operator is NOT associative
        uu  = scanl (\ uuim1 (ai, bi, cim1) -> 
                        let beta = ai / uuim1
                        in  bi - beta*cim1
                    ) uu0 (zip3 (tail a) (tail b) c)
 
        -- sequential scanl: binary operator is NOT associative
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
    in  reverse ur'


tridagPar :: [Double] -> [Double] -> [Double] -> [Double] 
          -> [Double]
tridagPar a b c y =
    ----------------------------------------------------
    -- Recurrence 1: b[i] = b[i] - a[i]*c[i-1]/b[i-1] --
    --   solved by scan with 2x2 matrix mult operator --
    ----------------------------------------------------
    let bfst = head b
        mats = map  (\ (bi,ai,cim1) -> (bi, 0.0-ai*cim1, 1.0, 0.0) )  
                    (zip3 (tail b) (tail a) c)

        scmt = -- parallel scan with 2x2 mat mult op
               scanl(\ (a0,a1,a2,a3) (b0,b1,b2,b3) -> 
		                let val = 1.0/(a0*b0) in
			            ( (b0*a0 + b1*a2)*val,
			              (b0*a1 + b1*a3)*val,
			              (b2*a0 + b3*a2)*val,
			              (b2*a1 + b3*a3)*val ) ) 
		            (1.0,0.0,0.0,1.0) mats

        b'   = map  (\ (t0,t1,t2,t3) -> (t0*bfst + t1) / (t2*bfst + t3) )
		            scmt

    ------------------------------------------------------
    -- Recurrence 2: y[i] = y[i] - (a[i]/b[i-1])*y[i-1] --
    --   solved by scan with linear func comp operator  --
    ------------------------------------------------------
        y0   = head y
        lffns= map (\ (yi,ai,bim1) -> (yi, 0.0-ai/bim1) )
                   (zip3 (tail y) (tail a) b')

        cffns= -- parallel scan with lin fun composition op
               scanl (\ (a0,a1) (b0,b1) -> (b0 + b1*a0, a1*b1) )
		             (0.0, 1.0)  lffns

        yf   = map  (\ (t1,t2) -> t1 + t2*y0 ) cffns

    ------------------------------------------------------
    -- Recurrence 3: backward recurrence solved via     --
    --             scan with linear func comp operator  --
    ------------------------------------------------------
        yn   = (last yf) / (last b')
        lbfn0= map (\ (yi,bi,ci) -> (yi/bi, 0.0-ci/bi) )
		           (tail $ zip3 (reverse yf) (reverse b') (reverse c))
        lbfns= (0.0, 1.0) : lbfn0
        cbfns= tail $ -- parallel scan with lin fun composition op
               scanl (\ (a0,a1) (b0,b1) -> (b0 + b1*a0, a1*b1) )
		             (0.0, 1.0)  lbfns

        yb   = map  (\ (t1,t2) -> t1 + t2*yn ) cbfns

    in  (reverse yb)

---------------------------------------------
---  Implicit Method (parameterization)   --- 
---------------------------------------------

implicitMethod :: Double
               -> [(Double,Double,Double)] 
               -> [(Double,Double,Double)] 
               -> [[Double]] -> [[Double]]
               -> [[Double]] -> [[Double]]
implicitMethod dtInv myD myDD myMu myVar u =
    zipWith3(\ mu_row var_row u_row ->
                let (a, b, c) = unzip3 $
                        zipWith4 (\ mu var (d0,d1,d2) (dd0,dd1,dd2) ->
                                    ( 0.0   - 0.5*(mu*d0 + 0.5*var*dd0)
			                        , dtInv - 0.5*(mu*d1 + 0.5*var*dd1)
			                        , 0.0   - 0.5*(mu*d2 + 0.5*var*dd2)
			                        )
                                 )  mu_row var_row myD myDD 
                in  if (1::Integer) == (1::Integer)
                    then tridagSeq a b c u_row
                    else tridagPar a b c u_row
            ) myMu myVar u

---------------------------------------------
---  Explicit Method (parameterization)   --- 
---------------------------------------------

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
--- UNUSED                             ---
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


---------------------------------------------
--- rollback: the brain of the program    ---
---------------------------------------------

rollback ::  Int -> [Double] -> [[Double]]
         -> [[Double]] -> [(Double,Double,Double)] -> [(Double,Double,Double)] -> [[Double]]
         -> [[Double]] -> [(Double,Double,Double)] -> [(Double,Double,Double)] -> [[Double]]
         -> [[Double]] 

rollback g  myTimeline  myResult
         myMuX  myDx  myDxx  myVarX
         myMuY  myDy  myDyy  myVarY = 

    let dtInv = 1.0 / ( (myTimeline !! (g+1)) - (myTimeline !! g) )

        u0= explicitXY dtInv  0.5  (transpose myResult)  myMuX  myVarX  myDx  myDxx
        v0= explicitXY 0.0    1.0  myResult  myMuY  myVarY  myDy  myDyy
 
        u1= map (\ (us, vs) -> zipWith (+) us vs )
                (zip u0 (transpose v0))

        u2= implicitMethod dtInv myDx myDxx myMuX myVarX u1
        y = zipWith (\ u_row v_row ->
			            zipWith (\ u_el v_el -> dtInv*u_el - 0.5*v_el) u_row v_row
			        ) (transpose u2) v0

    in  implicitMethod dtInv myDy myDyy myMuY myVarY y

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
    let (_, num_x, num_y, num_t) = it_spaces
        (s0, t, alpha, nu, beta) = params

        (myXindex, myYindex, myX, myY, myTimeline) =
            initGrid s0 alpha nu t num_x num_y num_t
        (myDx, myDxx) = initOperator myX
        (myDy, myDyy) = initOperator myY

        my_result_ini = map (\ xi -> replicate num_y (max (xi - strike) 0.0 )) myX
        
        loop_ros      = (myX, myDx, myDxx, myY, myDy, myDyy, myTimeline, alpha, beta, nu)

        my_result     = doLoop 0 (num_t-1) loop_ros my_result_ini 

    in  (my_result !! myXindex) !! myYindex  -- my_result[myXindex, myYindex]

------------------------------------------------
--- Main Entry Point: volatility calibration ---
------------------------------------------------
compute :: (Int,    Int,    Int,    Int) 
        -> (Double, Double, Double, Double, Double)
        -> [Double]
compute    it_spaces params =
    let -- (s0, t, alpha, nu, beta) = params
        (outer, _, _, _) = it_spaces -- was (outer, num_x, num_y, num_t)
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
                 s0    <- readDouble
                 t     <- readDouble
                 alpha <- readDouble
                 nu    <- readDouble
                 beta  <- readDouble

                 let params = (s0, t, alpha, nu, beta)

                 let v = compute (outer, num_x, num_y, num_t) params
                 r <- readDouble1d
                 return (v, r, (outer, num_x, num_y, num_t))

        readInt2d    = readArray $ readArray readInt
        readDouble1d = readArray readDouble
        readDouble2d = readArray $ readArray readDouble
        readDouble3d = readArray $ readArray $ readArray readDouble

-- ghc -O2 -msse2 -rtsopts  PricingLexiFi.hs
-- ./PricingLexiFi +RTS -K128m -RTS < ../Data/Medium/input.data
