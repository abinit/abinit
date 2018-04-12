$ = jQuery

$.extend true,$.ui.dialogExtend.prototype,
  modes:
    "collapse":
      option:"collapsable"
      state:"collapsed"
  options:
    "collapsable" : false
    "icons" :
      "collapse": "ui-icon-triangle-1-s"
    # callbacks
    "beforeCollapse" : null
    "collapse" : null

  collapse:()->
    newHeight = $(@element[0]).dialog("widget").find(".ui-dialog-titlebar").height()+15;
    # start!
    # trigger custom event
    @_trigger "beforeCollapse"
    # restore to normal state first (when necessary)
    unless @_state is "normal"
      @_restore()
    # remember original state
    @_saveSnapshot()
    pos = $(@element[0]).dialog("widget").position()
    $(@element[0])
      # modify dialog size (after hiding content)
      .dialog("option",
        "resizable" : false
        "height" : newHeight
        "maxHeight" : newHeight
        "position" : [pos.left - $(document).scrollLeft(),pos.top - $(document).scrollTop()]
      )
      .on('dialogclose',@_collapse_restore)
      # hide content
      # hide button-pane
      # make title-bar no-wrap
      .hide()
      .dialog("widget")
        .find(".ui-dialog-buttonpane:visible").hide().end()
        .find(".ui-dialog-titlebar").css("white-space", "nowrap").end()
      .find(".ui-dialog-content")
      # mark new state
      @_setState "collapsed"
      # modify dialog buttons according to new state
      @_toggleButtons()
      # trigger custom event
      @_trigger "collapse"

  _restore_collapsed:()->
    original = @_loadSnapshot()
    # restore dialog
    $(@element[0])
      # show content
      # show button-pane
      # fix title-bar wrap
      .show()
      .dialog("widget")
        .find(".ui-dialog-buttonpane:hidden").show().end()
        .find(".ui-dialog-titlebar").css("white-space", original.titlebar.wrap).end()
      .find(".ui-dialog-content")
      # restore config & size
      .dialog("option",
        "resizable" : original.config.resizable
        "height" : original.size.height
        "maxHeight" : original.size.maxHeight
      )
      .off('dialogclose',@_collapse_restore)

  _initStyles_collapse:()->
    if not $(".dialog-extend-collapse-css").length
      style = ''
      style += '<style class="dialog-extend-collapse-css" type="text/css">'
      style += '.ui-dialog .ui-dialog-titlebar-collapse { width: 19px; height: 18px; }'
      style += '.ui-dialog .ui-dialog-titlebar-collapse span { display: block; margin: 1px; }'
      style += '.ui-dialog .ui-dialog-titlebar-collapse:hover,'
      style += '.ui-dialog .ui-dialog-titlebar-collapse:focus { padding: 0; }'
      style += '</style>'
      $(style).appendTo("body")
  
  _collapse_restore:()->
    $(@).dialogExtend("restore")
