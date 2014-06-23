{-# LANGUAGE CPP,RecordWildCards #-}

-- running finpar benchmarks 
-- Usage: runhaskell benchmarks.hs

import System.Environment
import System.IO
import System.Process -- run commands
import System.Exit    -- exit codes
import System.Directory
import System.FilePath
import System.Console.GetOpt
import Control.Exception
import Text.ParserCombinators.ReadP -- (string, readP_to_S)
import Control.Monad
import Control.Applicative((<*), (<$>), (<*>), pure )
import Data.Maybe


------------------------------- options of this script: -------------------------
options =
    [ Option ['v'] ["verbose"] (NoArg (\o -> o{optVerbose = True}))
      "report activities on stderr"
    -- , Option ['g'] ["graph"]  (ReqArg (\f o -> o{optGraph = Just f}) "out.png")
    --   "output a png graph to given filename"
    , Option [] ["all"]  (NoArg (\o -> o{optAll = True}))
      "run all configured benchmarks"
    , Option ['s'] ["small"]   (NoArg (\o -> o{optSize = Small}))
      "run on small data (default)"
    , Option ['m'] ["medium"]  (NoArg (\o -> o{optSize = Medium}))
      "run on medium data"
    , Option ['l'] ["large"]   (NoArg (\o -> o{optSize = Large}))
      "run on large data"
    ]

data Options = Options
    { optVerbose :: Bool
    , optGraph :: Maybe FilePath
    , optAll :: Bool
    , optSize :: BenchSize
    }
instance Show Options where
    show Options{..} = "// Options:"
                       ++ "Running with " ++ show optSize ++ " data size, "
                       ++ (if optVerbose then "verbose, " else "")
                       -- ++ (case optGraph of
                       --      Just f  -> "output to " ++ f ++ ", "
                       --      Nothing -> ""
                       --    )
                       ++ (if optAll then "(all benchmarks)" else "")

defaults = Options { optVerbose = False, optGraph = Nothing
                   , optAll = False, optSize = Small }

data BenchSize = Small | Medium | Large
               deriving (Eq,Ord)
instance Show BenchSize where
    show Small  = "small"
    show Medium = "medium"
    show Large  = "large"
instance Read BenchSize where
    readsPrec _  = readP_to_S $ 
                   (string "small" >> return Small) +++
                   (string "medium" >> return  Medium)  +++
                   (string "large"  >> return Large)

---------------------------------------------------------------------------------

usage = usageInfo "runhaskell benchmarks.hs [OPTIONS] benchmarkname" options

main = do putStrLn "Please use make!"
          args <- getArgs
          let (opts,bs) = case getOpt Permute options args of
                            (o, [], []) -> 
                                (foldl (flip ($)) defaults{optAll=True} o, [])
                                    -- (flip id) :-]
                            (o, n, [])  -> (foldl (flip ($)) defaults o, n)
                            (_,_,errs)  -> error (concat errs ++ '\n':usage)
          benches <- filterM doesDirectoryExist 
                               (if optAll opts then dirnames else bs)
          when (optVerbose opts) $ do hPutStrLn stderr usage
                                      print benches
                                      print opts
          mapM_ (runBench opts) benches


---------------------------------------------------------------------------------
runBench :: Options -> String -> IO ()
runBench Options{..} dir =
    do (e,msg) <- runMake optVerbose dir []
       (ret,out) <- if e /= ExitSuccess then return (e,msg)
                    else runMake optVerbose dir ["run_" ++ show optSize]
       when (optVerbose || (ret /= ExitSuccess))
            (hPutStrLn stderr ("make returned " ++ show ret))
       -- parsing output 
       putStrLn out

runMake :: Bool -> FilePath -> [String] -> IO (ExitCode,String)
runMake verbose dir targets = do
  when verbose (putStrLn ("make " ++ unwords ("-C":dir:targets)))
  dirOK <- doesDirectoryExist dir
  if not dirOK then return (ExitFailure 1, dir ++ ": does not exist")
     else do
       let args = (if verbose then ["-C"] else ["--no-print-directory","-C"]) 
       (code, stdo, stde) <- readProcessWithExitCode "make" (args ++dir:targets) ""
       when verbose $ hPutStrLn stderr ("Output from " ++ dir ++ "\nstdout:\n" 
                                        ++ stdo ++ "\nstderr:\n" ++ stde)
       return (code,stdo)

-- every directory with source code must supply a makefile with
-- default target compiling the program, and targets run_small,
-- run_medium, run_large with test data in a relative path "../Data"

-- names of target directories (only those with code)
dirnames = ["CalibVolDiff" </> "Orig_COpenMP"
           , "CalibVolDiff" </> "Original"
           , "CalibVolDiff" </> "All_COpenMPCL"
           , "CalibVolDiff" </> "Outer_COpenMPCL"
           , "CalibVolDiff" </> "HaskellLH"
--           , "CalibGA/CppAndGPU"
--           , "CalibGA/python"
--           , "CalibGA/OCaml"
           , "GenericPricing" </> "Orig_COpenMp"
           , "GenericPricing" </> "CppOpenCL"
           , "GenericPricing" </> "HaskellLH"
           ]

-- names of different benchmarks
benchnames = [ "CalibVolDiff"
             , "GenericPricing"
             -- , "CalibGA"
             ]

-- code directories inside benchmark directories; kind-of ad-hoc standardised
codenames = ["Orig_COpenMP", "Original", "CppOpenCL"
            , "All_COpenMPCL", "Outer_COpenMPCL"
            , "HaskellLH"
            ] 


---------------------
-- parsing using ReadP (no external library required)

signed p =  (char '-' >> (('-':) <$> p)) +++ p

whitespace :: ReadP ()
whitespace = skipSpaces <*
             optional (string "//" >> manyTill anyChar (char '\n') >> whitespace)
anyChar = satisfy (const True)
digit = satisfy (`elem` ['0'..'9']) 

lexeme :: ReadP a -> ReadP a
lexeme p = p <* whitespace

-- token :: String -> Parser ()
token = return () . lexeme . string

readInt :: ReadP Int
readInt = lexeme $ read <$> signed (many1 digit)

readDouble :: ReadP Double
readDouble = read <$> signed (s2 +++ s1)
  where s1 = do bef <- many1 digit
                aft <- option "" ((:) <$> char '.' <*> many1 digit)
                return $ bef ++ aft
        s2 = (++) <$> (char '.' >> pure "0.") <*> many1 digit

-- not ready:
-- readArray :: ReadP a -> ReadP [a]
-- readArray p = lexeme $ between (token "[") (token "]") (p `sepBy` token ",")

