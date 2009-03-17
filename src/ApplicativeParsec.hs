{-
 - Applicative Parsec
 -
 - Define an applicative typeclass for the Parsec module, which represents an
 - applicative functor (??).
 -
 - Taken from Real World Haskell, Bryan O'Sullivan, chapter 16
 -}

 module ApplicativeParsec
    (
      module Control.Applicative
    , module Text.ParserCombinators.Parsec
    ) where

import Control.Applicative
import Control.Monad (MonadPlus(..), ap)

-- Hide a few names supplied by applicative
import Text.ParserCombinators.Parsec hiding (many, optional, (<|>))

-- The Applicative instance for every Monad looks like this
instance Applicative (GenParser s a) where
    pure  = return
    (<*>) = ap

-- The Alternative instance for every MonadPlus looks like this
instance Alternative (GenParser s a) where
    empty = mzero
    (<|>) = mplus

