'''
Created on Mar 1, 2010

@author: patrick
'''

from pyparsing import *

cvtInt = lambda toks: int(toks[0])
cvtReal = lambda toks: float(toks[0])
cvtTuple = lambda toks : tuple(toks.asList())
nameJoin = lambda toks : "".join([tok.replace("#","") for tok in toks[0]])
#lambda toks: " ".join([str(t) for t in toks[0]])

# define punctuation as suppressed literals
lparen,rparen,lbrack,rbrack,lbrace,rbrace,colon = map(Suppress,"()[]{}:")

integer = Combine(Optional(oneOf("+ -")) + Word(nums))\
    .setName("integer")\
    .setParseAction( cvtInt )
real = Combine(Optional(oneOf("+ -")) + Optional(Word(nums)) + "." +
               Optional(Word(nums)) +
               Optional(oneOf("e E")+Optional(oneOf("+ -")) +Word(nums))).setName("real").setParseAction( cvtReal )
file = Word(alphanums+"._-/")

class Feature(object):
    
    def __init__(self, **kwargs):
        for key in kwargs:
            self.__setattr__(key,kwargs[key])
            
    def __str__(self):
        vars = [str(key)+":"+str(self.__dict__[key]) for key in self.__dict__ if not key.startswith("_")]
        return "<Feature %s>" % (" ".join(vars))
            

class Consume(SkipTo):
    """
    Wrapper for SkipTo(*args,include=True)
    Advances the pointer and includes the expression.
    """
    def __init__(self, *args, **kwargs):
        kwargs['include']=True
        SkipTo.__init__(self, *args, **kwargs)
        