"Yet Another Python Templating Utility, Version 1.2. Taken from http://code.activestate.com/recipes/52305/"
from __future__ import print_function, division, absolute_import #, unicode_literals

import sys

from .six import exec_

#def exec_(stat, globals, locals):
#    """Comment the last line if you are using python3."""
#    # Py2 version
#    exec stat in globals, locals
#    # Py3 version
#    #exec(stat, globals, locals)


# utility stuff to avoid tests in the mainline code
class _nevermatch:
    "Polymorphic with a regex that never matches"
    def match(self, line):
        return None
_never = _nevermatch()     # one reusable instance of it suffices
def identity(string, why):
    "A do-nothing-special-to-the-input, just-return-it function"
    return string
def nohandle(string):
    "A do-nothing handler that just re-raises the exception"
    raise

# and now the real thing
class copier:
    "Smart-copier (YAPTU) class"
    def copyblock(self, i=0, last=None):
        "Main copy method: process lines [i,last) of block"
        def repl(match, self=self):
            "return the eval of a found expression, for replacement"
            # uncomment for debug: print '!!! replacing',match.group(1)
            expr = self.preproc(match.group(1), 'eval')
            try: return str(eval(expr, self.globals, self.locals))
            except: return str(self.handle(expr))
        block = self.locals['_bl']
        if last is None: last = len(block)
        while i<last:
            line = block[i]
            match = self.restat.match(line)
            if match:   # a statement starts "here" (at line block[i])
                # i is the last line to _not_ process
                stat = match.string[match.end(0):].strip()
                j=i+1   # look for 'finish' from here onwards
                nest=1  # count nesting levels of statements
                while j<last:
                    line = block[j]
                    # first look for nested statements or 'finish' lines
                    if self.restend.match(line):    # found a statement-end
                        nest = nest - 1     # update (decrease) nesting
                        if nest==0: break   # j is first line to _not_ process
                    elif self.restat.match(line):   # found a nested statement
                        nest = nest + 1     # update (increase) nesting
                    elif nest==1:   # look for continuation only at this nesting
                        match = self.recont.match(line)
                        if match:                   # found a contin.-statement
                            nestat = match.string[match.end(0):].strip()
                            stat = '%s _cb(%s,%s)\n%s' % (stat,i+1,j,nestat)
                            i=j     # again, i is the last line to _not_ process
                    j=j+1
                stat = self.preproc(stat, 'exec')
                stat = '%s _cb(%s,%s)' % (stat,i+1,j)
                # for debugging, uncomment...: print "-> Executing: {"+stat+"}"
                exec_(stat, self.globals, self.locals)
                i=j+1
            else:       # normal line, just copy with substitution
                self.ouf.write(self.regex.sub(repl,line))
                i=i+1
    def __init__(self, regex=_never, dict={},
            restat=_never, restend=_never, recont=_never,
            preproc=identity, handle=nohandle, ouf=sys.stdout):
        "Initialize self's attributes"
        self.regex   = regex
        self.globals = dict
        self.locals  = { '_cb':self.copyblock }
        self.restat  = restat
        self.restend = restend
        self.recont  = recont
        self.preproc = preproc
        self.handle  = handle
        self.ouf     = ouf
    def copy(self, block=None, inf=sys.stdin):
        "Entry point: copy-with-processing a file, or a block of lines"
        if block is None: block = inf.readlines()
        self.locals['_bl'] = block
        self.copyblock()

if __name__=='__main__':
    "Test: copy a block of lines, with full processing"
    import re
    rex=re.compile('@([^@]+)@')
    rbe=re.compile(r'\+')
    ren=re.compile('-')
    rco=re.compile('= ')
    x=23 # just a variable to try substitution
    cop = copier(rex, globals(), rbe, ren, rco)
    lines_block = [line+'\n' for line in """
A first, plain line -- it just gets copied.
A second line, with @x@ substitutions.
+ x+=1   # non-block statements MUST end with comments
-
Now the substitutions are @x@.
+ if x>23:
After all, @x@ is rather large!
= else:
After all, @x@ is rather small!
-
+ for i in range(3):
  Also, @i@ times @x@ is @i*x@.
-
One last, plain line at the end.""".split('\n')]
    print("*** input:")
    print(''.join(lines_block))
    print("*** output:")
    cop.copy(lines_block)


__doc__ = """
We may often want to copy some "template" text (normally from an input file) to an output
file-like object, while expanding Python expressions (and possibily executing Python
statements, e.g. for selection or repetition) that may be "embedded" in the template text.

YAPTU is a small but complete Python module for this purpose, suitable for processing most
any kind of structured-text input, since it lets client-code decide which
regular-expressions will denote embedded Python expressions and/or statements (so, such
re's can be selected to avoid conflicting with whatever syntax is needed by the kind of
structured-text that is being processed -- be it HTML, a programming language, RTF, ...).

The compiled-re object that identifies expressions, if not None, is used for a .sub on
each line of the input; for each MatchObject "match" that results, match.group(1) is
eval'd as a Python expression, and the result, transformed to a string, gets substituted
in place; a dictionary to be used as the namespace for the evaluation must also be passed
as an argument. Many such (non-overlapping!) matches per line are possible, but the
resulting text is NOT re-scanned for 'embeddings'.

Python statements can also be embedded; this is mostly intended to be used with
if/elif/else, for, while, and is line-based. Statement-related lines are recognized
through three more regular-expression objects that are passed in, one each for
'statement', 'continuation', 'finish', to be used for regular-expression _match_ (i.e.,
from line start) [again, each can be None if no such statements are to be embedded].

The 'stat' and 'cont' re's are followed by the corresponding statement lines (beginning
statement, and continuation statement, respectively -- the latter normally makes sense
only if it's an 'else' or 'elif'). Statements can nest without limits.

If a statement must be embedded that does NOT end with a colon (e.g., an assignment
statement), then a Python comment MUST terminate its line; conversely, such comments are
NOT allowed on the kind of statements most often embedded (if, else, for, while) --
_their_ lines must terminate with their ':' (optionally followed by whitespace). This
peculiarity is due to the somewhat tricky technique used in YAPTU's implementation:
embedded statements (with their continuations) are exec'd with _recursive calls to yaptu's
copyblock function_ substituted in place of the blocks of template-text they contain,
taking advantage of the fact that such a single-statement "suite" can be placed on the
same line as the controlling statement, right after the colon, and this avoids any
whitespace-issue (yaptu does NOT rely on whitespace to discern embedded-statement
structure, but on the explicit statement/continuation/end markers!).

Net of comments, whitespace, and docstrings, YAPTU is just 50 source-lines of code, but
rather a lot happens within that small compass. Instances of the _nevermatch auxiliary
class are used in lieu of regular-expression objects that are passed in as 'None' -- their
polymorphism with compiled-re objects (regarding the only two methods of the latter that
yaptu uses, .sub and .match) saves quite a few tests in the main body of code, and
simplifies it -- a good general idiom to keep in mind.

An instance of the 'copier' class has a certain amount of state, besides the relevant
compiled-re's (or nevermatch instances) and the output file-like object being used (the
latter need only implement method .write), that is held in two dictionary attributes --
self.globals, the dictionary that was originally passed in for expression-substitution,
and self.locals, another dictionary which is used as the local-namespace for all of
yaptu's exec and eval uses. Two internal-use-only items in self.locals, in particular
(with names starting with _) indicate the block of template-text being 'copied' (a
sequence of lines, each ending in a '\n'), at key '_bl', and the bound-method that
performs the copying, self.copyblock, at key '_cb'.

Holding these two pieces of state in self.locals items is not quaint personal usage --
it's part of the key to yaptu's workings, since self.locals is what is guaranteed to be
made available to the code that yaptu exec's (self.globals, too, but yaptu does NOT dirty
THAT dictionary -- it is owned by its caller!). Since .copyblock must be recursive (the
simplest way to ensure no nesting limitations), it is important that nested recursive
calls be always able to further recurse, if needed, through their exec statements! Access
to _bl is similarly necessary -- .copyblock only takes as arguments the line _indices_
inside _bl that a given recursive call is processing (in the usual form -- index of first
line to process, index of first following line to AVOID processing; i.e.,
lower-bound-included, upper-bound-excluded, as everywhere in Python).

copyblock's 32 SLOCs are the heart of YAPTU. The repl nested function is the one that is
passed to the .sub method of compiled RE objects to get the text to be used for each
expression substitution -- it uses eval on the expression string, and str() on the result
to ensure it, too, is turned back into a string. Most of copyblock is a simple while loop
that examines relevant block lines from the start; when it doesn't match a statement-start
RE, it copies the line out to the output file-object, with substitutions; when it does
match statement-start, it enters a smaller nested-loop looking for statement continuations
and statement-end (with proper accounting for nesting-levels, of course!). As it goes, it
builds up in local variable 'stat' a string, with the original statement [and possibly its
continuations at the same nesting-level] followed by a recursive call to _cb(i,j) right
after each semicolon [with newlines as separators between continuations, if any]. Lastly,
'stat' gets passed to the exec statement; the nested loop terminates, and the main loop
resumes from right after the embedded-statement just processed. Note that the
exec-statement will inevitably invoke copyblock recursively, but that does not disturb the
loops' state [based on local variables unoriginally named i and j, since they are
loop-counters and indices on the _bl list...], thanks to perfectly normal
recursive-invocation mechanisms.

Are the slightly-tricky subtleties in yaptu justified -- or does it violate the Prime
Directive, to "do the simplest thing that can possibly work"? I lean towards the former
opinion -- that there is no _gratuitous_ subtlety in yaptu, that it uses the minimal
amount of trickiness compatible with doing its job sensibly -- its job being to supply a
reusable templating utility, flexible, effective, and small, for a variety of actual and
potential uses; I can be biased, but I can see no _substantial_ simplification that could
be made for the sake of clarity without impairing functionality. However, your comments
(and proposed rewrites!) are welcome, of course!

Late-breaking additions: all expressions and statements may now be "preprocessed" by
passing an optional callable "preproc" when creating the copier -- default is no
preprocessing. Exceptions in expressions (only) may be handled by passing an optional
callable "handle" (default is re-raising the exception, which terminates YAPTU's
processing and propagates outwards).
"""
