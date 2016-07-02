#!/usr/bin/env python
# Author: Damien Caliste <damien.caliste@cea.fr>
# Copyright: CEA 2008
#  Read an ABINIT input file and display the attributes
#  values per dataset. When clicking on one attribute name,
#  the contains of the web help is displayed.

import ctypes
ctypes.PyDLL("./ab7_invars.so", mode=ctypes.RTLD_GLOBAL)
import ab7_invars as abinit
import numpy
import pygtk, gtk, gobject
import re, sys, os
import xml.etree.ElementTree as ElementTree

def delete_event(widget, event, data=None):
    return False

def destroy(widget, data=None):
    gtk.main_quit()

def loadData(widget, response_id, window):
  if (response_id != gtk.RESPONSE_ACCEPT):
    return
  if (widget is not window.get_data("fileSel")):
    widget.hide()

  # Remove old columns.
  view = window.get_data("treeview")
  for col in view.get_columns()[1:]:
    view.remove_column(col)

  wd = window.get_data("hboxResults")
  wd.set_sensitive(False)
  
  window.set_data("dtsets", None)
  # We parse the file.
  try:
      toto = abinit.Dtsets(widget.get_filename())
  except SyntaxError, err:
      msg = "While parsing '%s', ABINIT reports:\n<span font_desc=\"monospace\">" % widget.get_filename() + err.__str__() + "</span>"
  except Exception, err:
      msg = err.__str__()
  else:
      msg = None
  if (msg is not None):
      window.get_data("fileSel").unselect_all()
      wd = gtk.MessageDialog(window, gtk.DIALOG_MODAL, \
                             gtk.MESSAGE_ERROR, gtk.BUTTONS_OK, "Parse error")
      wd.format_secondary_markup(msg)
      wd.set_title("gtkParse error")
      wd.set_icon(window.get_icon())
      wd.run()
      wd.destroy()
      return
  window.set_data("dtsets", toto)

  wd = window.get_data("hboxResults")
  wd.set_sensitive(True)
  
  # Create a new model.
  args = "gtk.ListStore(gobject.TYPE_STRING"
  for i in range(toto.get_ndtset()):
    args += ", gobject.TYPE_STRING, gobject.TYPE_STRING"
  args += ")"
  model = eval(args)
  modelFilter = model.filter_new()
  window.set_data("model", model)
  modelFilter.set_visible_func(hideCycleIter, window)
  modelSort = gtk.TreeModelSort(modelFilter)
  modelSort.set_sort_column_id(0, gtk.SORT_ASCENDING)
  view.set_model(modelSort)

  for strId in abinit.get_ids().iterkeys():
    iter = model.append()
    model.set(iter, 0, strId)
    default = toto.get(strId, 0)
    for i in range(toto.get_ndtset()):
        val = toto.get(strId, i + 1)
        if (numpy.all(val == default)):
            misc = "#505050"
        else:
            misc = "blue"
        if (val is None):
            val = "<span size=\"smaller\"><i>Unsupported type</i></span>"
        model.set(iter, 2 * i + 1, val.__str__(), 2 * i + 2, misc)
  # Add the column to the tree view.
  for i in range(1, toto.get_ndtset() + 1):
    col = view.insert_column_with_attributes(-1, "Value dtset %d" % i, \
                                             gtk.CellRendererText(), \
                                             markup = 2 * i - 1, foreground = 2 * i)
    col.set_resizable(True)
    col.set_sort_column_id(i)

def hideCycleIter(model, iter, window):
    (strId,) = model.get(iter, 0)
    if (strId is None):
      return False
    if (window.get_data("hideDefault").get_active()):
      dtsets = window.get_data("dtsets")
      hide = True
      for i in range(dtsets.get_ndtset()):
        (color,) = model.get(iter, 2 * i + 2)
        hide = hide and (color != "blue")
    else:
      hide = False
    return (not(hide) and re.match(window.get_data("filter"), strId) is not None)

def applyFilter(widget, window, entry):
    patStr = entry.get_text()
    if (patStr[-1] != "$"):
        patStr += "$"
    pat = re.compile(patStr)
    if (pat is not None):
        window.set_data("filter", pat)
        window.get_data("treeview").get_model().get_model().refilter()

def applyHide(widget, window):
  window.get_data("treeview").get_model().get_model().refilter()

def withArg(argv, window):
    fileSel = window.get_data("fileSel")
    fileSel.set_filename(os.path.abspath(sys.argv[1]))
    while (gtk.events_pending()):
        gtk.main_iteration()
    if (fileSel.get_filename() is not None):
        loadData(fileSel, gtk.RESPONSE_ACCEPT, window)
    else:
        raise ValueError("No file '%s'." % os.path.abspath(sys.argv[1]))
    return False

def XMLGetText(var):
  if (var is None):
    return ""

  if (var.text is not None):
    text = var.text.strip().replace("\n", " ") + " "
  else:
    text = ""
  for child in var.getchildren():
    if (child.tag == "p"):
        text += XMLGetText(child) + "\n"
    elif (child.tag == "ul"):
        text += "\n" + XMLGetText(child)
    elif (child.tag == "li"):
        text += "  * " + XMLGetText(child) + "\n"
    elif (child.tag == "br"):
        text += "\n"
    else:
        text += XMLGetText(child)
  if (var.tail is not None):
    text += var.tail.strip().replace("\n", " ") + " "
  return re.sub("  *", " ", text)
  
def selectAtt(selection, window):
    (model, iter) = selection.get_selected()
    if (iter is None):
        return
    iSets = []
    (att,) = model.get(iter, 0)
    window.get_data("lblKeyword").set_markup("<b>" + att + "</b>")
    window.get_data("lblMnemonic").set_markup("")
    window.get_data("lblCategory").set_markup("")
    window.get_data("lblType").set_markup("")
    window.get_data("lblDefault").set_markup("")
    window.get_data("description").set_text("")
    varlist = window.get_data("varlist")
    for var in varlist:
      if (var.find("keyword").text == att):
        window.get_data("lblMnemonic").set_markup(XMLGetText(var.find("mnemonics")))
        window.get_data("lblCategory").set_markup(XMLGetText(var.find("categories")))
        window.get_data("lblType").set_markup(XMLGetText(var.find("type")))
        window.get_data("lblDefault").set_markup(XMLGetText(var.find("default")))
        window.get_data("description").set_text(XMLGetText(var.find("description")))
        break

if __name__ == "__main__":
  window = gtk.Window(gtk.WINDOW_TOPLEVEL)
  window.set_border_width(5)
  window.set_icon_from_file("../../../extras/logos/abinit-logo-64x64.png")
  window.connect("delete_event", delete_event)
  window.connect("destroy", destroy)

  vbox = gtk.VBox()
  window.add(vbox)

  hbox = gtk.HBox()
  vbox.pack_start(hbox, False, False, 5)

  wd = gtk.Label("<b>ABINIT input file:</b>")
  wd.set_use_markup(True)
  wd.set_alignment(0., 0.5)
  wd.set_padding(10, 0)
  hbox.pack_start(wd, False, False, 0)

  wd = gtk.FileChooserDialog(title = "Choose an ABINIT data file", \
                             buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_REJECT,
                                        gtk.STOCK_OPEN, gtk.RESPONSE_ACCEPT))
  wd.connect("response", loadData, window)
  fileSel = gtk.FileChooserButton(wd)
  window.set_data("fileSel", fileSel)
  hbox.pack_start(fileSel, True, True, 0)

  hbox = gtk.HPaned()
  hbox.set_sensitive(False)
  window.set_data("hboxResults", hbox)
  vbox.pack_start(hbox, True, True, 0)

  wd = gtk.VBox()
  hbox.pack1(wd, True, False)

  scroll = gtk.ScrolledWindow()
  scroll.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
  scroll.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
  wd.pack_start(scroll, True, True, 0)

  view = gtk.TreeView()
  wd.set_size_request(-1, 400)
  view.set_rules_hint(True)
  view.set_grid_lines(gtk.TREE_VIEW_GRID_LINES_VERTICAL)
  view.get_selection().set_mode(gtk.SELECTION_SINGLE)
  view.get_selection().connect("changed", selectAtt, window)
  window.set_data("treeview", view)
  scroll.add(view)

  col = view.insert_column_with_attributes(-1, "Attribute", \
                                           gtk.CellRendererText(), \
                                           markup=0)
  col.set_sort_column_id(0)

  hbox2 = gtk.HBox()
  wd.pack_start(hbox2, False, False, 0)
  entry = gtk.Entry()
  entry.set_text(".*")
  window.set_data("filter", re.compile(".*"))
  entry.connect("activate", applyFilter, window, entry)
  hbox2.pack_start(entry, True, True, 0)
  bt = gtk.Button(label="Filter")
  bt.connect("clicked", applyFilter, window, entry)
  hbox2.pack_start(bt, False, False, 0)
  bt = gtk.ToggleButton(label="Hide defaults")
  bt.connect("toggled", applyHide, window)
  window.set_data("hideDefault", bt)
  hbox2.pack_start(bt, False, False, 0)

  table = gtk.Table(5, 4)
  table.set_row_spacings(10);
  hbox.pack2(table, True, False)

  wd = gtk.Label("<i>Keyword:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 0, 1, 0, 1, xoptions = gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)
  wd = gtk.Label("")
  wd.set_alignment(0., 0.5)
  table.attach(wd, 1, 2, 0, 1, yoptions = gtk.SHRINK, xpadding = 5)
  window.set_data("lblKeyword", wd)

  wd = gtk.Label("<i>Mnemonic:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 2, 3, 0, 1, xoptions = gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)
  wd = gtk.Label("")
  wd.set_alignment(0., 0.5)
  wd.set_line_wrap(True)
  table.attach(wd, 3, 4, 0, 1, yoptions = gtk.SHRINK, xpadding = 5)
  window.set_data("lblMnemonic", wd)

  wd = gtk.Label("<i>Found in:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 0, 1, 1, 2, xoptions = gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)
  wd = gtk.Label("")
  wd.set_alignment(0., 0.5)
  table.attach(wd, 1, 2, 1, 2, yoptions = gtk.SHRINK, xpadding = 5)
  window.set_data("lblFound", wd)

  wd = gtk.Label("<i>Category:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 2, 3, 1, 2, xoptions = gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)
  wd = gtk.Label("")
  wd.set_alignment(0., 0.5)
  table.attach(wd, 3, 4, 1, 2, yoptions = gtk.SHRINK, xpadding = 5)
  window.set_data("lblCategory", wd)

  wd = gtk.Label("<i>Type:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 0, 1, 2, 3, xoptions = gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)
  wd = gtk.Label("")
  wd.set_alignment(0., 0.5)
  wd.set_line_wrap(True)
  table.attach(wd, 1, 2, 2, 3, yoptions = gtk.SHRINK, xpadding = 5)
  window.set_data("lblType", wd)

  wd = gtk.Label("<i>Default:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 2, 3, 2, 3, xoptions = gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)
  wd = gtk.Label("")
  wd.set_alignment(0., 0.5)
  wd.set_line_wrap(True)
  table.attach(wd, 3, 4, 2, 3, yoptions = gtk.SHRINK, xpadding = 5)
  window.set_data("lblDefault", wd)

  wd = gtk.Label("<i>Description:</i>")
  wd.set_use_markup(True)
  wd.set_alignment(1., 0.5)
  table.attach(wd, 0, 1, 3, 4, gtk.SHRINK | gtk.FILL, yoptions = gtk.SHRINK)

  scroll = gtk.ScrolledWindow()
  scroll.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
  scroll.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
  table.attach(scroll, 0, 4, 4, 5)
  buf = gtk.TextBuffer()
  window.set_data("description", buf)
  wd = gtk.TextView(buf)
  wd.set_editable(False)
  wd.set_wrap_mode(gtk.WRAP_WORD)
  wd.set_size_request(400, -1)
  scroll.add(wd)

  window.show_all()

  if (len(sys.argv) > 1):
    gobject.idle_add(withArg, sys.argv, window)

  # We parse the variable defs.
  try:
    varlist = ElementTree.parse("varlists.xml").findall("list/variable")
  except:
    varlist  =[]
  window.set_data("varlist", varlist)

  gtk.main()
