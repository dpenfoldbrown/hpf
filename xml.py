'''
Some ideas taken from
http://mail.python.org/pipermail/xml-sig/1998-October/000423.html
'''
from __future__ import absolute_import
import sys
from xml.sax import saxutils

import string
class DefaultXMLGenerator(saxutils.XMLGenerator):
    def __init__(self, out, pretty=False):
        saxutils.XMLGenerator.__init__(self,out)
        self.pretty = pretty
        self.indent = 0
        self.out = out
        self.char = False

    def characters(self,s):
        self.char = True
        if not isinstance(s,basestring):
            s = str(s)
        return saxutils.XMLGenerator.characters(self,s)

    def startElement(self,name,attrs=None):

        if self.pretty:
            self.out.write("\n")
            for i in xrange(self.indent):
                self.out.write("  ")

        saxutils.XMLGenerator.startElement(self, name,attrs if attrs else Attributes())
        self.indent=self.indent+1

    def endElement(self,name):
        self.indent=self.indent-1
        if self.pretty and not self.char:
            self.out.write("\n")
            for i in xrange(self.indent):
                self.out.write("  ")

        saxutils.XMLGenerator.endElement(self, name)
        self.char=False

class Attributes(dict):
    
    def __init__(self, **kwargs):
        for key in kwargs:
            kwargs[key]=str(kwargs[key])
        super(Attributes,self).__init__(**kwargs)

    def getLength(self):
        return len(self.keys)
    
    def getNames(self):
        return self.keys()
    
    def getType(self, name):
        return 'CDATA'

    def getValue(self, name):
        return self[name]

class Element(object):

    def __init__(self, handler, name, parent):
        self.__handler = handler
        self.__name = name
        self.__started = False
        self.__ended = False
        self.__children = []

    def __start__(self,**attrs):
        if self.__started:
            return
        self.__handler.startElement(self.__name,Attributes(**attrs))
        self.__started = True

    def __content__(self,*cdata):
        content = saxutils.escape(string.joinfields(cdata, ''))
        self.__handler.characters(content)

    def __getattr__(self, name):
        self.__start__()
        if name[:2] == '__':
            #return object.__getattribute__(self, name)
            raise AttributeError, name

        # Before creating and starting a new element...
        while(len(self.__children)>0):
            child = self.__children.pop()
            child.__end__()        

        el = Element(self.__handler,name, self)
        self.__children.append(el)
        return el

    def __call__(self, *cdata, **attrs):
        self.__start__(**attrs)
        self.__content__(*cdata)
        self.__end__()
        return self

    def __end__(self):
        self.__start__()
        if self.__ended:
            return

        # We close all children, and remove references as we do
        while(len(self.__children)>0):
            child = self.__children.pop()
            child.__end__()

        self.__handler.endElement(self.__name)
        self.__parent.__remove__(self)
        self.__ended = True

    def __remove__(self, child):
        if self.__children.count(child)>0:
            self.__children.remove(child)

    def __str__(self):
        s = '<' + self.__name
        for name, value in self.__attrs.items():
            s = s + ' ' + name + '="' + saxutils.escape(str(value)) + '"'
        if self.__cdata or self.__children:
            s = s + '>' + saxutils.escape(string.joinfields(self.__cdata, ''))
            for child in self.__children:
                s = s + str(child)
            s = s + '</' + self.__name + '>'
        else:
            s = s + '/>'
        return s

class SAXFactory:
    def __init__(self, handler):
        """
        @param hadler: ContentHandler for SAX events. Most likely an XMLGenerator
        to produce output as the document is constructed.
        @type handler: xml.sax.ContentHandler
        """
        self.__handler = handler
        self.__handler.startDocument()

    def __getattr__(self, name):
        if name[:2] == '__':
            raise AttributeError, name
        self.root = Element(self.__handler, name, self)
        return self.root

    def __remove__(self, child):
        pass

    def __end__(self):
        self.root.__end__()
        self.__handler.endDocument()

# class Intermediate:
#     def __init__(self, handler, parent, child, generator):
#         self.__handler = handler
#         self.__parent = parent
#         self.__child = child

#     def __getattr__(self, name):
#         inter = getattr(self.__child, name)
#         return Intermediate(self.__parent, inter.__child)

#     def __call__(self, *cdata, **attrs):
#         apply(self.__child, cdata, attrs)
#         return self

#     def __getitem__(self, items):
#         self.__child[items]
#         return self

#     def __str__(self):
#         return str(self.__parent)


def test():
    import sys
    handler = saxutils.XMLGenerator(sys.stdout)
    f = SAXFactory(handler)
    html = f.html
    head = html.head
    head.title("BLAH",id="2")
    head.other("test")
    ul = html.body.div.ul
    for i in range(10):
        li = ul.li(id=str(i))
    f.__end__()

