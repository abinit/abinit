$ = jQuery

$.extend true,$.ui.dialogExtend.prototype,
  modes:
    "minimize":
      option:"minimizable"
      state:"minimized"
  options:
    "minimizable" : false
    "minimizeLocation" : "left"
    "icons" :
      "minimize" : "ui-icon-minus"
    # callback
    "beforeMinimize" : null
    "minimize" : null
  
  minimize:()->
    # trigger custom event
    @_trigger "beforeMinimize"
    unless @_state is "normal"
      @_restore()
    # caculate new dimension
    newWidth = 200
    # create container for (multiple) minimized dialogs (when necessary)
    if $("#dialog-extend-fixed-container").length
      fixedContainer = $("#dialog-extend-fixed-container")
    else
      fixedContainer = $('<div id="dialog-extend-fixed-container"></div>').appendTo("body")
      fixedContainer.css
        "position" : "fixed"
        "bottom" : 1
        "left" : 1
        "right" : 1
        "z-index" : 9999
    # prepare dialog buttons for new state
    @_toggleButtons("minimized")
    dialogcontrols = $(@element[0]).dialog("widget").clone().children().remove().end()
    $(@element[0]).dialog("widget").find('.ui-dialog-titlebar').clone(true,true).appendTo(dialogcontrols)
    dialogcontrols.css
      # float is essential for stacking dialog when there are many many minimized dialogs
      "float" : @options.minimizeLocation,
      "margin" : 1
    fixedContainer.append(dialogcontrols)
    $(@element[0]).data("dialog-extend-minimize-controls",dialogcontrols)
    # disable draggable-handle (for <titlebar=none> only)
    if $(@element[0]).dialog("option","draggable")
      dialogcontrols.removeClass("ui-draggable")
    # modify dialogcontrols
    dialogcontrols.css
      "height": "auto"
      "width": newWidth
      "position": "static"
    # restore dialog before close
    $(@element[0]).on('dialogbeforeclose',@_minimize_restoreOnClose)
    # hide original dialog
    .dialog("widget").hide()
    # mark new state
    @_setState "minimized"
    # trigger custom event
    @_trigger "minimize"

  _restore_minimized:()->
    # restore dialog
    $(@element[0]).dialog("widget").show()
    # disable close handler
    $(@element[0]).off('dialogbeforeclose',@_minimize_restoreOnClose)
    # remove dialogcontrols
    $(@element[0]).data("dialog-extend-minimize-controls").remove()
    $(@element[0]).removeData("dialog-extend-minimize-controls")

  _initStyles_minimize:()->
    if not $(".dialog-extend-minimize-css").length
      style = ''
      style += '<style class="dialog-extend-minimize-css" type="text/css">'
      style += '.ui-dialog .ui-dialog-titlebar-minimize { width: 19px; height: 18px; }'
      style += '.ui-dialog .ui-dialog-titlebar-minimize span { display: block; margin: 1px; }'
      style += '.ui-dialog .ui-dialog-titlebar-minimize:hover,'
      style += '.ui-dialog .ui-dialog-titlebar-minimize:focus { padding: 0; }'
      style += '</style>'
      $(style).appendTo("body")
  
  _verifyOptions_minimize:()->
    if not @options.minimizeLocation or @options.minimizeLocation not in ['left','right']
      $.error( "jQuery.dialogExtend Error : Invalid <minimizeLocation> value '" + @options.minimizeLocation + "'" )
      @options.minimizeLocation = "left"
  
  _minimize_restoreOnClose:()->
    $(@).dialogExtend("restore")
