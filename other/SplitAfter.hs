--
-- Break a string into pieces after the given delimeter
--


simpleDigest      :: (a -> Bool) -> [a] -> [[a]]
simpleDigest p    =  addCleave . digest p

addCleave         :: [[a]] -> [[a]]
addCleave []      =  []
addCleave (x:[])  =  [x]
addCleave (x:y:z) =  x : (x++y) : addCleave (y:z)

digest      :: (a -> Bool) -> [a] -> [[a]]
digest _ [] =  []
digest p xs =  a : digest p b
    where
        (a,b) = splitAfter p xs

splitAfter          :: (a -> Bool) -> [a] -> ([a], [a])
splitAfter _ z@[]   =  (z, z)
splitAfter p (x:xs)
    | p x           =  ([x],xs)
    | otherwise     =  let (ys,zs) = splitAfter p xs in (x:ys,zs)

