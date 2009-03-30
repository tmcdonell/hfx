
import Data.ByteString.Internal (w2c)
import Data.ByteString.Lazy.Internal (ByteString(..))
import qualified Data.ByteString      as B
import qualified Data.ByteString.Lazy as L

digest :: (Char -> Bool) -> ByteString -> [ByteString]
digest _  Empty                 = []
digest p cs'@(Chunk c cs)       =
    case B.findIndex (p . w2c) c of
        Nothing                    -> [cs']
        Just n | n+1 == B.length c -> [cs']
               | otherwise         ->
               let (a,b) = B.splitAt (n+1) c in
                    Chunk a Empty : digest p (Chunk b cs)

addCleave   :: [ByteString] -> [ByteString]
addCleave l =  l ++ (zipper l (tail l))
    where
        zipper (a:as) (b:bs) = L.append a b : zipper as bs
        zipper _      _      = []
