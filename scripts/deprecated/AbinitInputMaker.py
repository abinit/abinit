#=================================
# AbinitInputMaker.py
version = '1.0'
#=================================
# last modified : january 17 2006
# written by Benjamin Tardif
# benjamin.tardif@umontreal.ca
#=================================
headline = '\n=====================\n AbinitInputMaker.py\n version %s\n=====================' %version

from Tkinter import *
from xml.sax import make_parser, SAXException
from xml.sax.handler import ContentHandler
import os # TEMP
import sys
from tkSimpleDialog import askstring
from tkMessageBox import showerror, showinfo

from class_ScrolledWidgetList import *

#font=(('Verbatim'),12,'bold')
#font=(('MS','Sans','Serif'),12)
font=(('Arial'),12)

################################################################################################################
# CONTENT HANDLER
################################################################################################################

class MAIN_ContentHandler(ContentHandler):
    def __init__(self):
        if debug.current:print 'MAIN_ContentHandler' # ---------- DEBUG
        locator.found = 0
        locator.tagnumberlist = []
        writeInputFileButton['state']='disabled'

    def startDocument(self):
        pass

    def startElement(self,Element,Attribute):
        locator.tagnumberlist.append(int(Attribute.get('tag')))
        locator.currenttag = [int(Attribute.get('tag')),'opening']
        if debug.prtcurrenttag:print 'locator.currenttag = ',locator.currenttag # ---------- DEBUG

        if locator.currenttag == locator.wantedtag:
            locator.found = 1
            if debug.prtlocatorfound:print 'locator found for currenttag',locator.currenttag # ---------- DEBUG
            locator.element = ''
            locator.description = ''
            locator.readdecision = 1

        if (locator.found==1) & (Element=='DECISION') & (locator.readdecision==0):
            locator.decisionindex += 1
            if debug.decisionindex:print 'locator.decisionindex (begin) ',locator.decisionindex # ---------- DEBUG

        if (locator.found==1) & (locator.readdecision==1):

            if Element == 'SECTION':
                locator.currentsection = str(Attribute.get('sectiontitle'))
                output.container.append('#=================================================\n# %s\n' % Attribute.get('sectiontitle'))
                output.pieceindex += 1
                UpdateOutput()
            
            if Element == 'CHOICE':
                locator.element = 'CHOICE'
                locator.instruction = Attribute.get('instruction')
                locator.wantedtag=locator.currenttag
                raise SAXException('') # stop parsing

            if Element == 'DECISION':
                if not int(Attribute.get('tag')) in locator.readdecisionlist:
                    locator.readdecision = 0
                    locator.decisionindex = 1

            if Element == 'MOUSEENTRY':
                locator.element = 'MOUSEENTRY'
                locator.instruction = Attribute.get('instruction')
                locator.wantedtag=locator.currenttag
                raise SAXException('') # stop parsing

            if Element == 'KEYBOARDENTRY':
                locator.element = 'KEYBOARDENTRY'
                locator.instruction = Attribute.get('instruction')
                locator.variablename = Attribute.get('variablename')
                locator.textlen = Attribute.get('textlen')
                if locator.textlen.split()[0] == 'getvalue':
                    locator.textlen = getvalue(locator.textlen.split()[1])
                    if debug.getvalue:print 'returned value',locator.textlen
                raise SAXException('') # stop parsing

            if Element == 'DIRECTENTRY':
                output.container.append('%s %s\n' % (Attribute.get('variablename'),Attribute.get('variablevalue')))
                output.pieceindex += 1
                UpdateOutput() 
              
    def characters(self,Content):
        pass
        
    def endElement(self,Element):
        locator.currenttag = [locator.tagnumberlist[-1],'closing']
        locator.tagnumberlist.pop(-1)
        if debug.prtcurrenttag:print 'locator.currenttag = ',locator.currenttag # ---------- DEBUG

        if locator.currenttag == locator.wantedtag:
            locator.found=1
            if debug.prtlocatorfound:print 'locator found for currenttag',locator.currenttag # ---------- DEBUG
        
        if locator.found:

            if (Element == 'DECISION') & (locator.readdecision==0):
                locator.decisionindex -= 1
                if debug.decisionindex:print 'locator.decisionindex ( end ) ',locator.decisionindex # ---------- DEBUG
        
            if (Element == 'DECISION') & (locator.decisionindex==0):
                locator.readdecision = 1

            if Element == 'SECTION':
                output.container.append('#=================================================\n\n\n')
                output.pieceindex += 1
                UpdateOutput()
                SECTIONENDERwidget(frameQuestion.container)
                frameQuestion.container.update()
                frameQuestion.refresh()
                
    def endDocument(self):
        pass

class CHOICE_ContentHandler(ContentHandler): # ne fait que construire la decisionlist utilise pour faire le CHOICEwidget
    def __init__(self):
        if debug.current:print 'CHOICE_ContentHandler' # ---------- DEBUG
        locator.found = 0
        locator.decisioncounter = 0
        locator.decisionlist = []
    def startElement(self,Element,Attribute):
        locator.currenttag = [int(Attribute.get('tag')),'opening']
        if locator.currenttag==locator.wantedtag:
            locator.found = 1
        if locator.found:
            if Element == 'DECISION':
                locator.decisioncounter += 1
            if Element == 'DECISION' and locator.decisioncounter == 1:
                locator.decisionlist.append([Attribute.get('tag'),Attribute.get('description')])
    def endElement(self,Element):
        if locator.found:
            if Element == 'DECISION':
                locator.decisioncounter -= 1
            if Element == 'CHOICE' and locator.decisioncounter == 0:
                raise SAXException('') # stop parsing

class MOUSEENTRY_ContentHandler(ContentHandler): # ne fait que construire la optionlist utilise pour faire le MOUSEENTRYwidget
    def __init__(self):
        if debug.current:print 'MOUSEENTRY_ContentHandler' # ---------- DEBUG
        locator.found = 0
        locator.optionlist = []
    def startElement(self,Element,Attribute):
        locator.currenttag = [int(Attribute.get('tag')),'opening']
        if locator.currenttag==locator.wantedtag:
            locator.found = 1
        if locator.found:
            if Element == 'OPTION':
                locator.optionlist.append([Attribute.get('description'),Attribute.get('variablename'),Attribute.get('variablevalue')])
    def endElement(self,Element):
        if locator.found:
            if Element == 'MOUSEENTRY':
                raise SAXException('') # stop parsing


def enablebuttons():
    L = len(sections.buttonlist)
    for i in range(L):
        if sections.buttonlist[i][1] == 'disabled':
            req = 0
            for requirement in sections.buttonlist[i][2]:
                if requirement in sections.completedsectionlist:
                    req += 1
            if req == len(sections.buttonlist[i][2]):
                sections.buttonlist[i][1] = 'enabled'
                sections.buttonlist[i][3]['state'] = 'normal'
    for section in sections.completedsectionlist:
        for i in range(L):
            if sections.buttonlist[i][0] == section:
                sections.buttonlist[i][1] = 'disabled'
                sections.buttonlist[i][3]['state'] = 'disabled'

def endsection(section):
    sections.completedsectionlist.append(section)
    output.permanent = textOutput.get(0.0,END) 
    enablebuttons()
    clear()
    UpdateOutput()
    writeInputFileButton['state'] = 'normal'


def parsefile(filename):
    clear()
    XMLfile.filename = filename
    MAINparser()
    #sections.completedsectionlist.append(filename)
    #enablebuttons()

class SECTIONS_ContentHandler(ContentHandler):
    def __init__(self):
        pass
    def startElement(self,Element,Attribute):
        if Element == 'SECTION':
            sections.buttonlist.append([str(Attribute.get('sectionname')),'disabled',str(Attribute.get('requirements')).split(),\
                Button(frameSection,state='disabled',text=Attribute.get('sectionname'),font=font,command=lambda:parsefile(str(Attribute.get('associatedfile'))))])
            sections.buttonlist[-1][-1].pack(side=TOP,padx=5,pady=5)
    def endElement(self,Element):
        pass

def MAINparser():
    if debug.current:print 'MAINparser' # ---------- DEBUG
    MAINparser = make_parser()
    MAINhandler = MAIN_ContentHandler()
    MAINparser.setContentHandler(MAINhandler)
    try:
        MAINparser.parse(open(sys.path[0]+'/'+XMLfile.filename))
    except SAXException:
        if locator.element == 'CHOICE':
            CHOICEparser()
        if locator.element == 'MOUSEENTRY':
            MOUSEENTRYparser()
        if locator.element == 'KEYBOARDENTRY':
            KEYBOARDENTRYwidget(frameQuestion.container,instruction=locator.instruction,variablename=locator.variablename,textlen=locator.textlen)
            frameQuestion.container.update()
            frameQuestion.refresh()

def CHOICEparser():
    if debug.current:print 'CHOICEparser' # ---------- DEBUG
    CHOICEparser = make_parser()
    CHOICEhandler = CHOICE_ContentHandler()
    CHOICEparser.setContentHandler(CHOICEhandler)
    try:
        CHOICEparser.parse(open(sys.path[0]+'/'+XMLfile.filename))
    except SAXException:
        CHOICEwidget(frameQuestion.container,instruction=locator.instruction,decisionlist=locator.decisionlist)
        frameQuestion.container.update()
        frameQuestion.refresh()

def MOUSEENTRYparser():
    if debug.current:print 'MOUSEENTRYparser' # ---------- DEBUG
    MOUSEENTRYparser = make_parser()
    MOUSEENTRYhandler = MOUSEENTRY_ContentHandler()
    MOUSEENTRYparser.setContentHandler(MOUSEENTRYhandler)
    try:
        MOUSEENTRYparser.parse(open(sys.path[0]+'/'+XMLfile.filename))
    except SAXException:
        MOUSEENTRYwidget(frameQuestion.container,instruction=locator.instruction,optionlist=locator.optionlist)
        frameQuestion.container.update()
        frameQuestion.refresh()

def getvalue(variable):
    for outputpiece in output.container:
        if outputpiece.split()[0] == variable:
            return str(abs(int(outputpiece.split()[1])))

################################################################################################################
# DEFINITIONS
################################################################################################################

class VariableContainer:
    pass

debug = VariableContainer()
debug.current = 0 # imprime le nom de chaque programme lors de son utilisation
debug.widget = 0  # imprime le W.index et le self.outputpiecebefore de chaque widget quand il est cree
debug.getvalue = 0 # imprime la valeur retournee par la fonction getvalue
debug.decisionindex = 0 # imprime le numero de l indentation des <DECISION> lors de la phase locator.readdecision = 0
debug.prtcurrenttag = 0 # imprime le locator.currenttag chaquefois qu'on frappe un element
debug.prtlocatorfound = 0 # imprime un message chaque fois que le locator est found dans le MAIN_ContentHandler

XMLfile = VariableContainer()
XMLfile.filename = 'PracticalExample.xml'
#XMLfile.filename = 'AbinitInputMaker8.xml'
#XMLfile.filename = 'keyboardtest.xml'

locator = VariableContainer()
locator.found = 0           # used in MAIN_ContentHandler, CHOICE_ContentHandler, MOUSEENTRY_ContentHandler
locator.currentsection = None # indique le nom de la section presente
locator.element = 'SECTION' # fait pour que le locator devienne TRUE des qu'il rencontre le premier widget <SECTION>
locator.description = None  # fait pour que le locator devienne TRUE des qu'il rencontre le premier widget <SECTION>
locator.instruction = None

locator.currenttag = None # current tag
locator.wantedtag = [1,'opening'] # wanted tag, on arrange le tout pour que ca devienne locator.found au premier element
locator.tagnumberlist = [] # list qui contient la liste de tout les element ouvert en ce moment

locator.decisioncounter = 0   # used in CHOICE_ContentHandler
locator.decisionlist = []     # used in CHOICEwidget, CHOICE_ContentHandler, CHOICEparser
locator.readdecisionlist = [] # used in MAINparser, contient la liste des readdecisionpiece attache a chacun des CHOICEwidget
locator.readdecision = 1      # used in MAINparser
locator.decisionindex = 0     # indique dans quelle indentation de decision nous sommes chaque fois qu on croise une widget DECISION

locator.currentmouseentry = None # used in MAIN_ContentHandler
locator.optionlist = []          # used in MOUSEENTRYwidget, MOUSEENTRY_ContentHandler, MOUSEENTRYparser

locator.currentkeyboardentry = None # used in MAIN_ContentHandler
locator.variablename = None         # used in KEYBOARDENTRYwidget, MAIN_ContentHandler
locator.textlen = None              # used in KEYBOARDENTRYwidget, MAIN_ContentHandler

W = VariableContainer()
W.index = 0 # index des widgets dans le framequestion.container (le premier a un index de 1)

output = VariableContainer()
output.container = [] # contient tous les outputpiece
output.pieceindex = 0 # index des outputpiece dans le output.container (la premiere a un pieceindex de 1 )
output.permanent = ''
output.inputfilename = None

sections = VariableContainer()
sections.buttonlist = [] # liste qui contient les widgets des boutons de chaque section
sections.completedsectionlist = ['None'] # liste qui contient toutes les sections qui sont terminees



colors = VariableContainer()
colors.framequestion = '#555555' # gris fonce
colors.widget_bg = '#cccccc' #gris pale
colors.widget_relief = 'raised'
colors.widget_bd = 8

colors.instruction_bg = 'blue'
colors.instruction_fg = 'white'
colors.instruction_relief = 'raised'
colors.instruction_bd = 2
colors.instruction_ipadx = 6
colors.instruction_ipady = 3
colors.decisionframe = '#dddddd' #gris pale
colors.optionframe = '#dddddd' #gris pale
colors.keyboardframe = '#dddddd' #gris pale

colors.option_fg = 'black'
colors.option_bg = 'white'
colors.option_active_fg = 'blue'
#colors.option_active_bg = 'black'
#colors.option_selected_fg = 'red'
colors.option_selected_bg = 'red'

class CHOICEwidget(Frame):    
    def __init__(self,master=None,instruction='<default>',decisionlist=['<default1>','<default2>','<default3>']): # default values, not used
       Frame.__init__(self,master,bg=colors.widget_bg,relief=colors.widget_relief,bd=colors.widget_bd)
       if debug.current:print 'CHOICEwidget' # ---------- DEBUG

       W.index += 1
       self.readdecisionpiece = None
       self.outputpiecebefore = output.pieceindex # indique le nombre de outputpiece present dans le output.container avant ce widget
       if debug.widget:print 'CHOICEwidget : W.index = %s , self.outputpiecebefore = %s' %(W.index,self.outputpiecebefore) # ---------- DEBUG

       self.instructionlabel = Label(self,text=instruction,font=font,bg=colors.instruction_bg,fg=colors.instruction_fg,relief=colors.instruction_relief,bd=colors.instruction_bd,wraplen=colors.wraplen)
       self.instructionlabel.pack(side=TOP,padx=10,pady=10,ipadx=colors.instruction_ipadx,ipady=colors.instruction_ipady)
       
       self.decisionframe = Frame(self,bg=colors.decisionframe,relief='solid',bd=1)
       self.decisionframe.pack(side=TOP,padx=10,pady=10)

       for i in range(len(decisionlist)):
           Label(self.decisionframe,text=decisionlist[i][1],font=font,anchor='w',relief='solid',bd=1,bg=colors.option_bg,fg=colors.option_fg,wraplen=colors.wraplen).pack(side=TOP,fill='x',padx=10,pady=5)
           self.decisionframe.winfo_children()[i].tag = int(decisionlist[i][0]) # attache le tag sur chaque label (seulement le numero du tag)
           self.decisionframe.winfo_children()[i].index = W.index # attache le W.index index sur chaque label

       for label in self.decisionframe.winfo_children():
           label.bind('<Enter>',self.labelenter)
           label.bind('<Leave>',self.labelleave)
           label.bind('<ButtonPress-1>',self.labelpress)

    def labelenter(self,event):
        event.widget.configure(fg=colors.option_active_fg)

    def labelleave(self,event):
        event.widget.configure(fg=colors.option_fg)

    def labelpress(self,event):
        for label in self.decisionframe.winfo_children():
            label.configure(bg=colors.option_bg) # change les couleurs
        event.widget.configure(bg=colors.option_selected_bg)

        while len(frameQuestion.container.winfo_children()) != event.widget.index:
            frameQuestion.container.winfo_children()[-1].destroy() # efface les widgets en dessous
            frameQuestion.container.update()
            frameQuestion.refresh()
            W.index -= 1

        while len(output.container) != frameQuestion.container.winfo_children()[event.widget.index-1].outputpiecebefore:
            output.container.remove(output.container[-1]) # efface les outputpiece en dessous
            output.pieceindex -= 1
            UpdateOutput()

        frameQuestion.container.winfo_children()[event.widget.index-1].readdecisionpiece = event.widget.tag # attache le tag du label selectionne sur le .readdecisionpiece de ce widget
        UpdateReadDecisionList()

        locator.wantedtag = [event.widget.tag,'opening']
        MAINparser()

class MOUSEENTRYwidget(Frame):
    def __init__(self,master=None,instruction='<empty>',optionlist=['<empty1>','<empty2>','<empty3>']): # default values, not used
       Frame.__init__(self,master,bg=colors.widget_bg,relief=colors.widget_relief,bd=colors.widget_bd)
       if debug.current:print 'MOUSEENTRYwidget' # ---------- DEBUG

       W.index += 1
       self.outputpiecebefore = output.pieceindex # indique le nombre de outputpiece present dans le output.container avant ce widget
       self.optionlist = locator.optionlist
       self.instruction = locator.instruction
       if debug.widget:print 'MOUSEENTRYwidget : W.index = %s , self.outputpiecebefore = %s' %(W.index,self.outputpiecebefore) # ---------- DEBUG

       self.instructionlabel = Label(self,text=instruction,font=font,bg=colors.instruction_bg,fg=colors.instruction_fg,relief=colors.instruction_relief,bd=colors.instruction_bd,wraplen=colors.wraplen)
       self.instructionlabel.pack(side=TOP,padx=10,pady=10,ipadx=colors.instruction_ipadx,ipady=colors.instruction_ipady)
       
       self.optionframe = Frame(self,bg=colors.optionframe,relief='solid',bd=1)
       self.optionframe.pack(side=TOP,padx=10,pady=10)

       for i in range(len(self.optionlist)):
           Label(self.optionframe,text=self.optionlist[i][0],font=font,anchor='w',relief='solid',bd=1,bg=colors.option_bg,fg=colors.option_fg,wraplen=colors.wraplen).pack(side=TOP,fill='x',padx=10,pady=5)

       for label in self.optionframe.winfo_children():
           label.index = W.index # attache le W.index index sur chaque label
           label.tag = locator.currenttag[0] # attache le tag sur le label (seulement le numero du tag)
           label.bind('<Enter>',self.labelenter)
           label.bind('<Leave>',self.labelleave)
           label.bind('<ButtonPress-1>',self.labelpress)

    def labelenter(self,event):
        event.widget.configure(fg=colors.option_active_fg)

    def labelleave(self,event):
        event.widget.configure(fg=colors.option_fg)

    def labelpress(self,event):
        for label in self.optionframe.winfo_children():
            label.configure(bg=colors.option_bg)
        event.widget.configure(bg=colors.option_selected_bg)

        while len(frameQuestion.container.winfo_children()) != event.widget.index:
            frameQuestion.container.winfo_children()[-1].destroy() # efface les widgets en dessous
            frameQuestion.container.update()
            frameQuestion.refresh()
            W.index -= 1

        while len(output.container) != frameQuestion.container.winfo_children()[event.widget.index-1].outputpiecebefore:
            output.container.remove(output.container[-1]) # efface les outputpiece en dessous
            output.pieceindex -= 1
            UpdateOutput()

        for (description,variablename,variablevalue) in self.optionlist:
            if event.widget['text'] == description:
                output.container.append('%s %s\n' %(variablename,variablevalue))
                output.pieceindex += 1 # imprime la variable dans le output
                UpdateOutput()

        locator.wantedtag = [event.widget.tag,'closing']
        MAINparser()


class KEYBOARDENTRYwidget(Frame):
    def __init__(self,master=None,instruction='<empty>',variablename='<empty>',textlen=1): # default values
       Frame.__init__(self,master,bg=colors.widget_bg,relief=colors.widget_relief,bd=colors.widget_bd)
       if debug.current:print 'KEYBOARDENTRYwidget' # ---------- DEBUG

       W.index += 1
       self.outputpiecebefore = output.pieceindex # indique le nombre de outputpiece present dans le output.container avant ce widget
       self.instruction = locator.instruction
       self.variablename = locator.variablename
       self.textlen = locator.textlen
       if debug.widget:print 'KEYBOARDENTRYwidget : W.index = %s , self.outputpiecebefore = %s' %(W.index,self.outputpiecebefore) # ---------- DEBUG


       self.instructionlabel = Label(self,text=instruction,font=font,bg=colors.instruction_bg,fg=colors.instruction_fg,relief=colors.instruction_relief,bd=colors.instruction_bd,wraplen=colors.wraplen)
       self.instructionlabel.pack(side=TOP,padx=10,pady=10,ipadx=colors.instruction_ipadx,ipady=colors.instruction_ipady)

       self.keyboardframe = Label(self,bg=colors.keyboardframe,relief='solid',bd=1)
       self.keyboardframe.pack(side=TOP,padx=10,pady=10)

       self.text = Text(self.keyboardframe,width=30,height=locator.textlen,bg=colors.option_bg,fg=colors.option_fg)
       self.text.pack(side=TOP,padx=10,pady=10)

       self.button = Button(self.keyboardframe,bg=colors.instruction_bg,fg=colors.instruction_fg,activebackground=colors.instruction_bg,activeforeground=colors.instruction_fg,text='Validate',font=font,command=self.validate)
       self.button.index = W.index # attache le W.index sur chaque le bouton
       self.button.pack(side=TOP,padx=10,pady=10)
       self.button.tag = locator.currenttag[0] # attache le tag sur le bouton (seulement le numero du tag)

    def validate(self):
        while len(frameQuestion.container.winfo_children()) != self.button.index:
            frameQuestion.container.winfo_children()[-1].destroy() # efface les widgets en dessous
            frameQuestion.container.update()
            frameQuestion.refresh()
            W.index -= 1
        
        while len(output.container) != frameQuestion.container.winfo_children()[self.button.index-1].outputpiecebefore:
            output.container.remove(output.container[-1]) # efface les outputpiece en dessous
            output.pieceindex -= 1
            UpdateOutput()

        if self.textlen == '1':
            output.container.append('%s %s' % (self.variablename,self.text.get(0.0,END)))
            output.pieceindex += 1
            UpdateOutput()
        else:
            output.container.append('%s\n%s' % (self.variablename,self.text.get(0.0,END)))
            output.pieceindex += 1
            UpdateOutput()

        locator.wantedtag = [self.button.tag,'closing']
        MAINparser()

class SECTIONENDERwidget(Frame):    
    def __init__(self,master=None):
       Frame.__init__(self,master,bg=colors.widget_bg,relief=colors.widget_relief,bd=colors.widget_bd)
       if debug.current:print 'SECTIONENDERwidget' # ---------- DEBUG

       W.index += 1
       self.outputpiecebefore = output.pieceindex # indique le nombre de outputpiece present dans le output.container avant ce widget
       if debug.widget:print 'SECTIONENDERwidget : W.index = %s , self.outputpiecebefore = %s' %(W.index,self.outputpiecebefore) # ---------- DEBUG

       self.endbutton = Button(self,text='end this section',font=font,bg=colors.instruction_bg,fg=colors.instruction_fg,activebackground=colors.instruction_bg,activeforeground=colors.instruction_fg,command=lambda: endsection(locator.currentsection))
       self.endbutton.pack(side=TOP,padx=10,pady=10)


def UpdateOutput():
    if debug.current:print 'UpdateOutput' # ---------- DEBUG
    textOutput.configure(state=NORMAL)
    textOutput.delete(0.0,END)
    textOutput.insert(END,output.permanent)
    for outputpiece in output.container:
        textOutput.insert(END,outputpiece)
    textOutput.configure(state=DISABLED)
    textOutput.yview('moveto',1)
    
    

def UpdateReadDecisionList():
    locator.readdecisionlist = []
    for widget in frameQuestion.container.winfo_children():
        try:
            locator.readdecisionlist.append(widget.readdecisionpiece) # parcourt les widgets et ajoute les elements a la readdecision liste pour les elements opportuns            
        except:
            pass


################################################################################################################
# MAIN CLASS
################################################################################################################
print headline
print '\nloadind graphical interface ...\n'

#root===========================================================================================================
root = Tk()
root.title('AbinitInputMaker.py version %s' %version)
xmax = int(0.85*int(str(root.wm_maxsize()).replace('(','').replace(')','').replace(',','').split()[0]))
ymax = int(0.85*int(str(root.wm_maxsize()).replace('(','').replace(')','').replace(',','').split()[1]))
root.geometry('+0+0')
#===============================================================================================================

Frame(root,height=20,width=20).grid(row=0,column=0) # only for aesthetics

#frameSection==================================================================================================
frameSection = Frame(root,bg='black',relief='solid',bd=1)
Label(frameSection,text='SECTIONS',bg='black',fg='white',font=font).pack(padx=5,pady=5)
#===============================================================================================================

#compute sections.buttonlist and create buttons
SECTIONSparser = make_parser()
SECTIONShandler = SECTIONS_ContentHandler()
SECTIONSparser.setContentHandler(SECTIONShandler)
SECTIONSparser.parse(open(sys.path[0]+'/Sections.xml'))

#equalize all buttons width
width = []
for button in sections.buttonlist:
    width.append(len(button[0]))
width = max(width)
for button in sections.buttonlist:
    button[3]['width']=width

#frameOutput====================================================================================================
frameOutput = Frame(root,relief='solid',bd=1)
#===============================================================================================================

#textOutput=====================================================================================================
textOutput = Text(frameOutput,width=51,font=font,height=5,bg='white',relief='groove',wrap='none',state=DISABLED)
textOutput.grid(row=0,column=0)

textOutputScrollbarX = Scrollbar(frameOutput,orient=HORIZONTAL,command=textOutput.xview,width=16 )
textOutputScrollbarX.grid(row=1,column=0,sticky=EW)

textOutputScrollbarY = Scrollbar(frameOutput,orient=VERTICAL,command=textOutput.yview,width=16)
textOutputScrollbarY.grid(row=0,column=1,sticky=NS)

textOutput["xscrollcommand"]  =  textOutputScrollbarX.set
textOutput["yscrollcommand"]  =  textOutputScrollbarY.set
#===============================================================================================================

#compute frameQuestion width and height
frameSection.grid(row=1,column=0,padx=20,pady=20)
frameOutput.grid(row=1,column=2,padx=20,pady=0)
root.update()

width=xmax-80-frameSection.winfo_reqwidth()-frameOutput.winfo_reqwidth()
height=ymax-80

#frameQuestion==================================================================================================
frameQuestion = ScrolledWidgetList(root,width=width,height=height,bg=colors.framequestion,pady=10,bottomspace=0)
#===============================================================================================================

frameQuestion.grid(row=1,column=1,padx=0,pady=0)
colors.wraplen = width-100

#ClearButton====================================================================================================
clearButton = Button(root,text='Clear',font=font,command=lambda:clear())
#===============================================================================================================

#writeInputFileButton===========================================================================================
writeInputFileButton = Button(root,text='Write Input File',font=font,state='disabled',command=lambda:writeinputfile())
#===============================================================================================================

#clear==========================================================================================================
def clear():
    while len(frameQuestion.container.winfo_children()) != 0:
        frameQuestion.container.winfo_children()[-1].destroy() # efface les widgets en dessous
        frameQuestion.container.update()
        try:
            frameQuestion.refresh()
        except:
            pass
        W.index -= 1
        
    while len(output.container) != 0:
        output.container.remove(output.container[-1]) # efface les outputpiece en dessous
        output.pieceindex -= 1
        UpdateOutput()

    locator.wantedtag = [1,'opening'] # wanted tag, on arrange le tout pour que ca devienne locator.found au premier element
    if output.permanent != '':
        writeInputFileButton['state']='normal'
    
#===============================================================================================================

#writeinputfile=================================================================================================
def writeinputfile():
    output.inputfilename = askstring('AbinitInputMaker.py','Enter the name of the input file you wish to produce')
    if output.inputfilename != None:
        if output.inputfilename[-3:] != '.in':
            showerror('AbinitInputMaker.py','The input file must be a *.in file')
            output.inputfilename = None
        else:
            writer = open(output.inputfilename,'w')
            writer.write('%s' %textOutput.get(0.0,END))
            writer.close
            showinfo('AbinitInputMaker.py','input file produced successfully:\n--> %s <--' %output.inputfilename)
#===============================================================================================================

clearButton.grid(row=2,column=1,padx=20,pady=20)
writeInputFileButton.grid(row=2,column=2,padx=20,pady=20)

#adjust textOutput height
height1 = frameOutput.winfo_reqheight()
textOutput['height']=int(textOutput['height'])+1
root.update()
height2 = frameOutput.winfo_reqheight()
dh = height2-height1

count=0
while height2 < height:
    count+=1
    height2 += dh

textOutput['height'] = int(textOutput['height'])+count-1

enablebuttons()

root.mainloop()
