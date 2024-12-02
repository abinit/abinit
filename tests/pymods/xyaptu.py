"XYAPTU: Lightweight XML/HTML Document Template Engine for Python. Taken from http://code.activestate.com/recipes/162292/"
from __future__ import print_function, division, absolute_import #, unicode_literals

__version__ = '1.0.0'
__author__= [
  'Alex Martelli (aleax@aleax.it)',
  'Mario Ruggier (mario@ruggier.org)'
]
__copyright__ = '(c) Python Style Copyright. All Rights Reserved. No Warranty.'
__dependencies__ = ['YAPTU 1.2, http://aspn.activestate.com/ASPN/Python/Cookbook/Recipe/52305']
__history__= {
  '1.0.0' : '2002/11/13: First Released Version',
}

####################################################

import sys, re, string

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


from .yaptu import copier

class xcopier(copier):
    ' xcopier class, inherits from yaptu.copier '

    def __init__(self, dns, rExpr=None, rOpen=None, rClose=None, rClause=None,
                 ouf=sys.stdout, dbg=0, dbgOuf=sys.stdout):
        ' set default regular expressions required by yaptu.copier '

        # Default regexps for yaptu delimeters (what xyaptu tags are first converted to)
        # These must be in sync with what is output in self._x2y_translate
        _reExpression = re.compile('_:@([^:@]+)@:_')
        _reOpen       = re.compile(r'\++yaptu ')
        _reClose      = re.compile('--yaptu')
        _reClause     = re.compile('==yaptu ')

        rExpr         = rExpr  or _reExpression
        rOpen         = rOpen  or _reOpen
        rClose        = rClose or _reClose
        rClause       = rClause or _reClause

        # Debugging
        self.dbg = dbg
        self.dbgOuf = dbgOuf
        _preproc = self._preProcess
        if dbg: _preproc = self._preProcessDbg

        # Call super init
        copier.__init__(self, rExpr, dns, rOpen, rClose, rClause,
                        preproc=_preproc, handle=self._handleBadExps, ouf=ouf)


    def xcopy(self, input=None):
        '''
        Converts the value of the input stream (or contents of input filename)
        from xyaptu format to yaptu format, and invokes yaptu.copy
        '''

        # Read the input
        inf = input
        try:
            inputText = inf.read()
        except AttributeError:
            inf = open(input)
            if inf is None:
                raise ValueError("Can't open file (%s)" % input)
            inputText = inf.read()
        try:
            inf.close()
        except:
            pass

        # Translate (xyaptu) input to (yaptu) input, and call yaptu.copy()

        yinf = StringIO(self._x2y_translate(inputText))
        self.copy(inf=yinf)
        yinf.close()

    def _x2y_translate(self, xStr):
        ' Converts xyaptu markup in input string to yaptu delimeters '

        # Define regexps to match xml elements on.
        # The variations (all except for py-expr, py-close) we look for are:
        # <py-elem code="{python code}" /> |
        # <py-elem code="{python code}">ignored text</py-elem> |
        # <py-elem>{python code}</py-elem>

        # ${py-expr} | $py-expr | <py-expr code="pvkey" />
        reExpr = re.compile(r'''
          \$\{([^}]+)\} |  # ${py-expr}
          \$([_\w]+) | # $py-expr
          <py-expr\s+code\s*=\s*"([^"]*)"\s*/> |
          <py-expr\s+code\s*=\s*"([^"]*)"\s*>[^<]*</py-expr> |
          <py-expr\s*>([^<]*)</py-expr\s*>
        ''', re.VERBOSE)

        # <py-line code="pvkeys=pageVars.keys()"/>
        reLine = re.compile(r'''
          <py-line\s+code\s*=\s*"([^"]*)"\s*/> |
          <py-line\s+code\s*=\s*"([^"]*)"\s*>[^<]*</py-line> |
          <py-line\s*>([^<]*)</py-line\s*>
        ''', re.VERBOSE)

        # <py-open code="for k in pageVars.keys():" />
        reOpen = re.compile(r'''
          <py-open\s+code\s*=\s*"([^"]*)"\s*/> |
          <py-open\s+code\s*=\s*"([^"]*)"\s*>[^<]*</py-open\s*> |
          <py-open\s*>([^<]*)</py-open\s*>
        ''', re.VERBOSE)

        # <py-clause code="else:" />
        reClause = re.compile(r'''
          <py-clause\s+code\s*=\s*"([^"]*)"\s*/> |
          <py-clause\s+code\s*=\s*"([^"]*)"\s*>[^<]*</py-clause\s*> |
          <py-clause\s*>([^<]*)</py-clause\s*>
        ''', re.VERBOSE)

        # <py-close />
        reClose = re.compile(r'''
          <py-close\s*/> |
          <py-close\s*>.*</py-close\s*>
        ''', re.VERBOSE)

        # Call-back functions for re substitutions
        # These must be in sync with what is expected in self.__init__
        def rexpr(match,self=self):
            return '_:@%s@:_' % match.group(match.lastindex)
        def rline(match,self=self):
            return '\n++yaptu %s #\n--yaptu \n' % match.group(match.lastindex)
        def ropen(match,self=self):
            return '\n++yaptu %s \n' % match.group(match.lastindex)
        def rclause(match,self=self):
            return '\n==yaptu %s \n' % match.group(match.lastindex)
        def rclose(match,self=self):
            return '\n--yaptu \n'

        # Substitutions
        xStr = reExpr.sub(rexpr, xStr)
        xStr = reLine.sub(rline, xStr)
        xStr = reOpen.sub(ropen, xStr)
        xStr = reClause.sub(rclause, xStr)
        xStr = reClose.sub(rclose, xStr)

        # When in debug mode, keep a copy of intermediate template format
        if self.dbg:
            _sep = '====================\n'
            self.dbgOuf.write('%sIntermediate YAPTU format:\n%s\n%s' % (_sep, xStr, _sep))

        return xStr

    # Handle expressions that do not evaluate
    def _handleBadExps(self, s):
        ' Handle expressions that do not evaluate '
        if self.dbg:
            self.dbgOuf.write('!!! ERROR: failed to evaluate expression: %s \n' % s)
        return '***! %s !***' % s

    # Preprocess code
    def _preProcess(self, s, why):
        ' Preprocess embedded python statements and expressions '
        return self._xmlDecode(s)
    def _preProcessDbg(self, s, why):
        ' Preprocess embedded python statements and expressions '
        self.dbgOuf.write('!!! DBG: %s %s \n' % (s, why))
        return self._xmlDecode(s)

    # Decode utility for XML/HTML special characters
    _xmlCodes = [
      ['"', '&quot;'],
      ['>', '&gt;'],
      ['<', '&lt;'],
      ['&', '&amp;'],
    ]
    def _xmlDecode(self, s):
        ' Returns the ASCII decoded version of the given HTML string. '
        codes = self._xmlCodes
        for code in codes:
            #s = string.replace(s, code[1], code[0])
            s = s.replace(code[1], code[0])
        return s


####################################################

if __name__=='__main__':

    ##################################################
    # Document Name Space (a dictionary, normally prepared by runtime application,
    # and that serves as the substitution namespace for instantiating a doc template).
    #
    DNS = {
      'pageTitle' : 'Event Log (xyaptu test page)',
      'baseUrl' : 'http://xproject.sourceforge.net/',
      'sid' : 'a1b2c3xyz',
      'session' : 1,
      'userName' : 'mario',
      'startTime' : '12:31:42',
      'AllComputerCaptions' : 'No',
      'ComputerCaption' : 'mymachine01',
      'LogSeverity' : ['Info', 'Warning', 'Error' ],
      'LogFileType' : 'Application',
      'logTimeStamp' : 'Event Log Dump written on 25 May 2001 at 13:55',
      'logHeadings' : ['Type', 'Date', 'Time', 'Source', 'Category', 'Computer', 'Message'] ,
      'logEntries' : [
        ['Info', '14/05/2001', '15:26', 'MsiInstaller', '0', 'PC01', 'winzip80 install ok...'],
        ['Warning', '16/05/2001', '02:43', 'EventSystem', '4', 'PC02', 'COM+ failed...'],
        ['Error', '22/05/2001', '11:35', 'rasctrs', '0', 'PC03', '...', ' ** EXTRA ** ' ],
      ]
    }

    # and a function...
    def my_current_time():
        import time
        return str(time.clock())
    DNS['my_current_time'] = my_current_time

    '''
    # To use functions defined in an external library
    import externalFunctionsLib
    dict['fcn'] = externalFunctionsLib
    # which will therefore permit to call functions with:
    ${fcn.somefun()}
    '''

    ##################################################
    # Sample page template that uses the xyaptu tags and pcdata expressions.
    # Note that:
    #  - source code indentation here is irrelevant for xyaptu
    #  - xyaptu tags may span more than one source line
    #
    templateString = '''<html>
   <head>
    <title>$pageTitle</title>
   </head>
   <body bgcolor="#FFFFFF" text="#000000">

    <py-open code="if session:"/>
     Logged on as $userName, since <py-expr>startTime</py-expr>
     (<a href="$baseUrl?sid=$sid&amp;linkto=Logout">Logout?</a>)
    <py-close/>
    <hr>
    <h1>${pageTitle}</h1>
    <hr>
    <p>${a bad expression}</p>
    <p>
     <b>Filtering Event Log With:</b><br>
     All Computers: $AllComputerCaptions <br>
     Computer Name: $ComputerCaption <br>
     Log Severity:
      <py-open code="for LG in LogSeverity:"/>
        $LG
      <py-close/>
      <br>
     Log File Type: <py-expr code="LogFileType" />
    </p>
    <hr>
    <p>$logTimeStamp</p>

    <table width="100%" border="0" cellspacing="0" cellpadding="2">

     <tr valign="top" align="left">
      <py-open code = "for h in logHeadings:" > code attribute takes precedence
       over this text, which is duly ignored </py-open>
       <th>$h</th>
      <py-close/>
     </tr>

     <py-line code = "numH=len(logHeadings)" />

     <py-open code="for logentry in logEntries:"/>
      <tr valign="top" align="left">
       <py-open>for i in range(0,len(logentry)):</py-open>
        <py-open code="if i &lt; numH:" />
         <td>${logentry[i]}</td>
        <py-clause code="else:" />
         <td bgcolor="#cc0000">Oops! <!-- There's more log entry fields than headings! --></td>
        <py-close/>
       <py-close>### close (this is ignored) </py-close>
      </tr>
     <py-close/>

    </table>
    <hr>
    Current time: ${my_current_time()}
    <hr>
   </body>
  </html>
    '''

    ##################################################
    # Set a filelike object to templateString
    templateStream = StringIO(templateString)

    ##################################################
    # Initialise an xyaptu xcopier, and call xcopy
    xcp = xcopier(DNS)
    xcp.xcopy(templateStream)


    ##################################################
    # Test DBG 1
    # Set dbg ON (writing dbg statements on output stream)
    '''
    xcp = xcopier(DNS, dbg=1)
    xcp.xcopy(templateStream)
    '''

    ##################################################
    # Test DBG 2
    # Write dbg statements to a separate dbg stream
    '''
    dbgStream = StringIO()
    dbgStream.write('DBG info: \n')
    xcp = xcopier(DNS, dbg=1, dbgOuf=dbgStream)
    xcp.xcopy(templateStream)
    print dbgStream.getvalue()
    dbgStream.close()
    '''

####################################################

__doc__ = """
Xyaptu is python-centric, in the sense that the XML tags offered reflect python constructs
(such as python expressions, statements, opening and closing blocks) and not particularly
constructs typically identified in web page templates. The advantage is simplicity, while
still keeping all the options open for the HTML or XML document designer.

The primary requirements of xyaptu are:

(a) expression evaluation, e.g. variable substitutions, function calls (b) loop over, and
format, a python data sequence (c) xyaptu mark-up must pass through an XML parser, to
naturally allow using XML tools of choice, such as XSLT, for generation of page templates
(d) but, since HTML is not XML, page templates need not otherwise be XML compliant. In
some future time xhtml may be the norm for web pages, but as yet eb design tools currently
in wide use do not produce XML compliant output. Thus, non-XML page templates, such as
HTML, must still be considered as valid input. (For the implementation, this implies that
xyaptu tags be matched using regular expressions, and not by parsing HTML or XML.) (e)
simplicity of use, with minimum learning and runtime overhead (f) separation of
presentation and application logic

There are only 5 XML tags, to handle python statements (expression, line, block open,
block continuation clause, block close). Python expressions are also mapped to a 'pcdata'
token, to allow the use of python expressions also in places where tags are not allowed,
i.e. in attribute values. XML special characters (< > & ") must be encoded (< > & ") to be
used in python code. This, unfortunately, is unavoidable.

Xyaptu may be run in debug mode, and debug output may be sent to either the specified
output filelike object, or to any writable filelike object. Please see the module
self-test for sample code. When in debug mode, the intermediate format of the template is
also copied to the debug stream (done in copier._x2y_translate). As a default behaviour,
expressions that do not evaluate are written out (surrounded with '! ' and ' !') to the
specified output stream, and, if in debug mode, an error message is written out to the
debug stream (which defaults to the output stream). To change this behaviour (and that of
debug in general) you would need to override the methods _handleBadExps and _preprocessDbg
in a sublcassed xyaptu.

Mark-up syntax: A template may contain anything acceptable to targeted clients will
accept, plus the following xyaptu tags and 1 expression, to mark-up python expressions and
statements:

<py-expr code="pvkey" /> -- expression <py-line code="pvkeys=pVars.keys()" /> -- line
statement <py-open code="if inSession:"/> -- block open <py-clause code="else:"/> -- block
continuation clause <py-close/> -- block close

$pyvar | ${pyexpr} -- for simple interpolation of python variables, expressions, or
function calls (the second syntax option is for expressions that may contain spaces or
other characters not allowed in variable names)

The advantage of pcdata tokens over XML tags is that they may be used anywhere, including
within XML attribute values, e.g.: <a href="http://host/$language/topic">topic</a>

Alternate mark-up tag syntax: Because most web browsers by default do not display
attribute values, but they do show element values, as a convenience for those who like to
preview page templates in web browsers, an alternate tag syntax is provided. The five tag
examples above therefore become:

<py-expr>pvkey -- expression <py-line>pvkeys=pVars.keys() -- line statement <py-open>if
inSession: -- block open <py-clause>else: -- block continuation clause <py-close># close
-- block close (content is ignored)

Note that, in the case that a "code" attribute is specified, then _that_ code is executed,
and element content, if any, is ignored.

Usage: (1) A runtime application prepares the document namespace in the form of a python
dictionary. (2) An xyaptu xcopier is initialised with this dictionary: from xyaptu import
xcopier xcp = xcopier(DNS) (3) The xcopy method of this xcopier instance may be called
with either the name of a page template file, or a filehandle, to instantiate the page
template within this namespace: xcp.xcopy( templateFileName | templateFileHandle ) For
page templates available as strings in memory, use StringIO: from cStringIO import
StringIO pageTemplateStream = StringIO(pageTemplateString) xcp.xcopy(pageTemplateStream)
(4) Output is by default sent to sys.stdout. A different output stream may be specified at
initialisation, as the value of an 'ouf' parameter: xcp = xcopier(DNS, ouf=myOutputStream)
(5) Debugging may be turned on by initialising xyaptu with a 'dbg=1' switch. Debug output
is sent to the specified output stream unless a separate stream is specified by the dbgOuf
parameter, also at initialisation time.

For a full working example, see the module self-test source.

Enhancements to consider:

Add support for XML namespaces (py:expr, py:line, ...)
For each python statement (open, close, ...), yaptu adds an extra blank line. To be
'faithful' to the template source, it would be better if this is not so.
Do not process xyaptu mark-up when this is inside an XML comment in the source document
template
Add possibility to include xyaptu mark-up as verbatim document content, i.e. to be able to
write out ${pyexpr} as is.
"""
