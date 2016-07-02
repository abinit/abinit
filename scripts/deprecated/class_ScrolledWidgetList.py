from Tkinter import *

class ScrolledWidgetList(Frame):
    
    def __init__(self,master=None,width=100,height=100,bg='white',pady=10,topspace=0,bottomspace=0):
        
        self.width = width
        self.height = height
        self.pady = pady
        self.topspace = topspace
        self.bottomspace = bottomspace

        Frame.__init__(self,master)

        self.border = Frame(self,relief='solid',bd=1)
        self.border.grid(row=0,column=0)

        self.container =  Canvas(self.border,width=width,height=height,bg=bg,scrollregion=(0,0,0,height),highlightthickness=0)
        self.container.pack()

        self.scrollbarY = Scrollbar(self,orient=VERTICAL,command=self.scrollY)
        self.scrollbarY.grid(row=0,column=1,sticky=NS)
        self.container['yscrollcommand'] = self.scrollbarY.set
   
    def scrollY(self,mode=None,value=None,units=None):
        self.container.yview(mode,value,units)
        self.container.update()
        self.placewidget()        
            
    def refresh(self):
        try:
            if len(self.container.winfo_children()) > len(self.Y):  # a widget has been added
                self.Y.append([self.Y[-1][0]+self.Y[-1][1]+self.pady,self.container.winfo_children()[-1].winfo_reqheight()])
            elif len(self.container.winfo_children()) < len(self.Y): # one or more widgets has been deleted
                while len(self.container.winfo_children()) < len(self.Y):
                    self.Y.remove(self.Y[-1])
        except: # this is the very first widget to be added
            self.Y = []
            self.Y.append([self.pady+self.topspace,self.container.winfo_children()[0].winfo_reqheight()])
        
        self.container.configure(scrollregion=(0,0,0,self.Y[-1][0]+self.Y[-1][1]+self.pady+self.bottomspace))

        self.container.yview('moveto',1)
        self.container.update()
        self.placewidget()
 
    def placewidget(self):
        self.container.yview('moveto',self.scrollbarY.get()[0])
        Ymodificator = self.scrollbarY.get()[0]*float(self.container['scrollregion'].split()[3])
        for i in range(len(self.container.winfo_children())):
            self.container.winfo_children()[i].place(relx=0.5,y=self.Y[i][0]-Ymodificator,anchor=N)
