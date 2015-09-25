'''
Created on Feb 26, 2010

@author: Patrick
'''

from pyparsing import *
from hpf.parsing import *

class PositiveSelectionParser(object):

    MODELS = "models"
    POSITIVES = "positives"
    COLUMN = "column"
    AA = "aa"
    PROBABILITY = "Pr(w>1)"
    POST_MEAN = "post mean"
    PLUS_MINUS = "+-"
    STANDARD_ERROR = "SE for w"
    STAR = "*"
    MODEL = "Positively Selected Sites : Model" 




#    # The section title name
#    selection_title = (Literal("Bayes Empirical Bayes (BEB)")+
#                       Consume(lparen+
#                              Literal("amino acids refer to 1st sequence:")+
#                              Consume(rparen)).suppress()
#                              )
    
    # The column header line
    selection_header = (Literal(PROBABILITY)+
                        Literal(POST_MEAN)+
                        Literal(PLUS_MINUS)+
                        Literal(STANDARD_ERROR)+
                        lineEnd+lineEnd
                        )
    # Line defining sites of positive selection
    selection_line = Group(integer.setResultsName(COLUMN)+
                      Word(alphas+"*",exact=1).setResultsName(AA)+
                      real.setResultsName(PROBABILITY)+ZeroOrMore(Literal(STAR)).setResultsName(STAR)+
                      real.setResultsName(POST_MEAN)+
                      (
                       (Literal(PLUS_MINUS)+real.setResultsName(STANDARD_ERROR)+lineEnd)
                       | lineEnd
                      )
                      )
    
    
        # One model of selected sites
    model = Group(Consume(Literal(MODEL)).suppress()+integer.setResultsName(MODEL)+
                  Consume(selection_header).suppress()+
                  ZeroOrMore(selection_line).setResultsName(POSITIVES)
                  )
    # file is defined as a set of models
    expression = ZeroOrMore(model).setResultsName(MODELS)
    
#    expression = (Consume(selection_title).suppress()+
#                  Consume(selection_header.suppress())+
#                  ZeroOrMore(selection_line).setResultsName(POSITIVES)
#                  )
#    
#    sd = (Literal("Bayes Empirical Bayes (BEB)")+
#                       SkipTo(Literal("sorghum#15448)")))
#    test = SkipTo(selection_line).suppress()
        
    def parse(self, input, debug=False):
        """
        CodeML objects have the 'positive_selection' attribute which is a list
        of PositiveSelection objects.
        @type input: open file handle type object.
        @return generator of CodeML objects.
        """
        from hpf.hddb.db import PositiveSelection, CodeML
        ps = PositiveSelectionParser
        ps.expression.setDebug(debug)
        result = ps.expression.parseFile(input)
        result = result.asDict()
        if result.has_key(ps.MODELS):
            for m in result[ps.MODELS]:
                m = m.asDict()
                model = CodeML(model=m[ps.MODEL])
                model.ps = []
                if m.has_key(ps.POSITIVES):
                    for site in m.get(ps.POSITIVES):
                        vars = {'column':site.get(ps.COLUMN),
                                'aa':site.get(ps.COLUMN),
                                'probability':site.get(ps.PROBABILITY),
                                'post_mean':site.get(ps.POST_MEAN),
                                'stderr':site.get(ps.STANDARD_ERROR)
                                }
                        # This stderr measure seems to be optional for
                        # some models.
                        model.ps.append(PositiveSelection(**vars))
                yield model
            
if __name__=="__main__":
    import sys
    parser = PositiveSelectionParser()
    for file in sys.argv[1:]:
        print "parsing",file
        with open(file) as handle:
            for vars in parser.parse(handle):
                print "\t",vars
