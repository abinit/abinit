$ = jQuery

$.extend true,$.ui.dialogExtend.prototype,
  modes:
    "maximize":
      option:"maximizable"
      state:"maximized"
  options:
    "maximizable" : false
    "icons" :
      "maximize" : "ui-icon-extlink"
    # callbacks
    "beforeMaximize" : null
    "maximize" : null

  maximize:()->
    newHeight = $(window).height()-11;
    newWidth = $(window).width()-11;
    # start!
    # trigger custom event
    @_trigger "beforeMaximize"
    # restore to normal state first (when necessary)
    unless @_state is "normal"
      @_restore()
    # remember original state
    @_saveSnapshot()
    # disable draggable-handle (for <titlebar=none> only)
    if $(@element[0]).dialog("option","draggable")
      $(@element[0])
      .dialog("widget")
        .draggable("option", "handle", null)
        .find(".ui-dialog-draggable-handle").css("cursor", "text").end()
    $(@element[0])
      # fix dialog from scrolling
      .dialog("widget")
        .css("position", "fixed")
      .find(".ui-dialog-content")
      # show content
      # show button-pane (when minimized/collapsed)
      .show()
      .dialog("widget")
        .find(".ui-dialog-buttonpane").show().end()
      .find(".ui-dialog-content")
      # modify dialog with new config
      .dialog("option",
        "resizable" : false
        "draggable" : false
        "height" : newHeight
        "width" : newWidth
        "position" :
            my: "left top"
            at: "left top"
            of: window
      )
      # mark new state
      @_setState "maximized"
      # modify dialog buttons according to new state
      @_toggleButtons()
      # trigger custom event
      @_trigger "maximize"
  _restore_maximized:()->
    original = @_loadSnapshot()
    # restore dialog
    $(@element[0])
      # free dialog from scrolling
      # fix title-bar wrap (if dialog was minimized/collapsed)
      .dialog("widget")
        .css("position", original.position.mode)
        .find(".ui-dialog-titlebar").css("white-space", original.titlebar.wrap).end()
      .find(".ui-dialog-content")
      # restore config & size
      .dialog("option",
        "resizable" : original.config.resizable
        "draggable" : original.config.draggable
        "height" : original.size.height
        "width" : original.size.width
        "maxHeight" : original.size.maxHeight
        "position" :
          my: "left top"
          at: "left+"+original.position.left+" top+"+original.position.top
          of: window
      )
      # restore draggable-handle (for <titlebar=none> only)
      if $(@element[0]).dialog("option","draggable")
        $(@element[0])
        .dialog("widget")
          .draggable("option", "handle", if $(@element[0]).dialog("widget").find(".ui-dialog-draggable-handle").length then $(@element[0]).dialog("widget").find(".ui-dialog-draggable-handle") else".ui-dialog-titlebar")
          .find(".ui-dialog-draggable-handle")
          .css("cursor", "move");

  _initStyles_maximize:()->
    if not $(".dialog-extend-maximize-css").length
      style = ''
      style += '<style class="dialog-extend-maximize-css" type="text/css">'
      style += '.ui-dialog .ui-dialog-titlebar-maximize { width: 19px; height: 18px; }'
      style += '.ui-dialog .ui-dialog-titlebar-maximize span { display: block; margin: 1px; }'
      style += '.ui-dialog .ui-dialog-titlebar-maximize:hover,'
      style += '.ui-dialog .ui-dialog-titlebar-maximize:focus { padding: 0; }'
      style += '</style>'
      $(style).appendTo("body")