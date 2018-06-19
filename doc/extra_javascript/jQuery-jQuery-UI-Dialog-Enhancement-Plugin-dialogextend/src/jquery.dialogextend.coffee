$ = jQuery

$.widget "ui.dialogExtend",
  version: "2.0.0"
  
  modes:{}
  options:
    "closable" : true
    "dblclick" : false
    "titlebar" : false
    "icons" :
      "close" : "ui-icon-closethick"
      "restore" : "ui-icon-newwin"
    # callbacks
    "load" : null
    "beforeRestore" : null
    "restore" : null

  _create: ()->
    @_state = "normal"
    if not $(@element[0]).data "ui-dialog"
      $.error "jQuery.dialogExtend Error : Only jQuery UI Dialog element is accepted" 
    @_verifyOptions()
    @_initStyles()
    @_initButtons()
    @_initTitleBar()
    @_setState "normal"
    @_on "load",(e)->
      console.log "test",e
    @_trigger "load"
  
  _setState: (state)->
    $(@element[0])
    .removeClass("ui-dialog-"+@_state)
    .addClass("ui-dialog-"+state)
    @_state = state
  
  _verifyOptions: ()->
    # check <dblckick> option
    if @options.dblclick and @options.dblclick not of @modes
      $.error "jQuery.dialogExtend Error : Invalid <dblclick> value '" + @options.dblclick + "'"
      @options.dblclick = false
    # check <titlebar> option
    if @options.titlebar and @options.titlebar not in ["none","transparent"]
      $.error "jQuery.dialogExtend Error : Invalid <titlebar> value '" + @options.titlebar + "'"
      @options.titlebar = false;
    # check modules options
    for name of @modes
      if @["_verifyOptions_"+name] then @["_verifyOptions_"+name]()
  
  _initStyles:()->
    if not $(".dialog-extend-css").length
      style = ''
      style += '<style class="dialog-extend-css" type="text/css">'
      style += '.ui-dialog .ui-dialog-titlebar-buttonpane>a { float: right; }'
      style += '.ui-dialog .ui-dialog-titlebar-restore { width: 19px; height: 18px; }'
      style += '.ui-dialog .ui-dialog-titlebar-restore span { display: block; margin: 1px; }'
      style += '.ui-dialog .ui-dialog-titlebar-restore:hover,'
      style += '.ui-dialog .ui-dialog-titlebar-restore:focus { padding: 0; }'
      style += '.ui-dialog .ui-dialog-titlebar ::selection { background-color: transparent; }'
      style += '</style>'
      $(style).appendTo("body")
    for name of @modes
      @["_initStyles_"+name]()

  _initButtons:()->
    # start operation on titlebar
    titlebar = $(@element[0]).dialog("widget").find ".ui-dialog-titlebar"
    # create container for buttons
    buttonPane = $('<div class="ui-dialog-titlebar-buttonpane"></div>').appendTo titlebar
    buttonPane.css
      "position" : "absolute"
      "top" : "50%"
      "right" : "0.3em"
      "margin-top" : "-10px"
      "height" : "18px"
    # move 'close' button to button-pane
    titlebar
      .find(".ui-dialog-titlebar-close")
        # override some unwanted jquery-ui styles
        .css
          "position" : "relative"
          "float" : "right"
          "top" : "auto"
          "right" : "auto"
          "margin" : 0
        # change icon
        .find(".ui-icon").removeClass("ui-icon-closethick").addClass(@options.icons.close).end()
        # move to button-pane
        .appendTo(buttonPane)
      .end()
    # append restore button to button-pane
    buttonPane
      .append('<a class="ui-dialog-titlebar-restore ui-corner-all ui-state-default" href="#"><span class="ui-icon '+@options.icons.restore+'" title="restore">restore</span></a>')
      # add effect to button
      .find('.ui-dialog-titlebar-restore')
        .attr("role", "button")
        .mouseover(()-> $(@).addClass("ui-state-hover"))
        .mouseout(()-> $(@).removeClass("ui-state-hover"))
        .focus(()-> $(@).addClass("ui-state-focus"))
        .blur(()-> $(@).removeClass("ui-state-focus"))
      .end()
      # default show buttons
      # set button positions
      # on-click-button
      .find(".ui-dialog-titlebar-close")
        .toggle(@options.closable)
      .end()
      .find(".ui-dialog-titlebar-restore")
        .hide()
        .click((e)=>
          e.preventDefault()
          @restore()
        )
      .end();
      # add buttons from modules
      for name,mode of @modes
        @_initModuleButton name,mode

    # other titlebar behaviors
    titlebar
      # on-dblclick-titlebar
      .dblclick((evt)=>
        if @options.dblclick
          if @_state != "normal"
            @restore()
          else
            @[@options.dblclick]()
      )
      # avoid text-highlight when double-click
      .select(()->
        return false
      )
  
  _initModuleButton:(name,mode)->
    buttonPane = $(@element[0]).dialog("widget").find '.ui-dialog-titlebar-buttonpane'
    buttonPane.append('<a class="ui-dialog-titlebar-'+name+' ui-corner-all ui-state-default" href="#" title="'+name+'"><span class="ui-icon '+@options.icons[name]+'">'+name+'</span></a>')
      .find(".ui-dialog-titlebar-"+name)
        .attr("role", "button")
        .mouseover(()-> $(@).addClass("ui-state-hover"))
        .mouseout(()-> $(@).removeClass("ui-state-hover"))
        .focus(()-> $(@).addClass("ui-state-focus"))
        .blur(()-> $(@).removeClass("ui-state-focus"))
      .end()
      .find(".ui-dialog-titlebar-"+name)
        .toggle(@options[mode.option])
        .click((e)=>
          e.preventDefault()
          @[name]()
        )
        .end()

  _initTitleBar:()->
    switch @options.titlebar
        when false then 0
        when "none"
          # create new draggable-handle as substitute of title bar
          if $(@element[0]).dialog("option", "draggable")
            handle = $("<div />").addClass("ui-dialog-draggable-handle").css("cursor", "move").height(5)
            $(@element[0]).dialog("widget").prepend(handle).draggable("option", "handle", handle);
          # remove title bar and keep it draggable
          $(@element[0])
            .dialog("widget")
            .find(".ui-dialog-titlebar")
              # clear title text
              .find(".ui-dialog-title").html("&nbsp;").end()
              # keep buttons at upper-right-hand corner
              .css(
                "background-color" : "transparent"
                "background-image" : "none"
                "border" : 0
                "position" : "absolute"
                "right" : 0
                "top" : 0
                "z-index" : 9999
              )
            .end();
        when "transparent"
          # remove title style
          $(@element[0])
            .dialog("widget")
            .find(".ui-dialog-titlebar")
            .css(
              "background-color" : "transparent"
              "background-image" : "none"
              "border" : 0
            )
        else
          $.error( "jQuery.dialogExtend Error : Invalid <titlebar> value '" + @options.titlebar + "'" );

  state:()->
    return @_state

  restore:()->
    # trigger custom event
    @_trigger "beforeRestore"
    @_restore()
    # modify dialog buttons according to new state
    @_toggleButtons()
    # trigger custom event
    @_trigger "restore"
  
  _restore:()->
    unless @_state is "normal"
      @["_restore_"+@_state]()
      # mark new state
      @_setState "normal"
      # return focus to window
      $(@element[0]).dialog("widget").focus()
  
  _saveSnapshot:()->
    if @_state is "normal" 
      @original_config_resizable = $(@element[0]).dialog("option", "resizable")
      @original_config_draggable = $(@element[0]).dialog("option", "draggable")
      @original_size_height = $(@element[0]).dialog("widget").outerHeight()
      @original_size_width = $(@element[0]).dialog("option", "width")
      @original_size_maxHeight = $(@element[0]).dialog("option", "maxHeight")
      @original_position_mode = $(@element[0]).dialog("widget").css("position")
      @original_position_left = $(@element[0]).dialog("widget").offset().left-$('body').scrollLeft()
      @original_position_top = $(@element[0]).dialog("widget").offset().top-$('body').scrollTop()
      @original_titlebar_wrap = $(@element[0]).dialog("widget").find(".ui-dialog-titlebar").css("white-space")

  _loadSnapshot:()->
    {
      "config" :
        "resizable" : @original_config_resizable
        "draggable" : @original_config_draggable
      "size" :
        "height" : @original_size_height
        "width"  : @original_size_width
        "maxHeight" : @original_size_maxHeight
      "position" :
        "mode" : @original_position_mode
        "left" : @original_position_left
        "top"  : @original_position_top
      "titlebar" :
        "wrap" : @original_titlebar_wrap
    }
  
  _toggleButtons:(newstate)->
    state = newstate or @_state
    $(@element[0]).dialog("widget")
      .find(".ui-dialog-titlebar-restore")
        .toggle( state != "normal" )
        .css({ "right" : "1.4em" })
      .end()
    for name,mode of @modes
      $(@element[0]).dialog("widget")
      .find(".ui-dialog-titlebar-"+name)
      .toggle( state != mode.state && @options[mode.option] )
    # place restore button after current state button
    for name,mode of @modes
      if mode.state is state
        $(@element[0]).dialog("widget")
          .find(".ui-dialog-titlebar-restore")
            .insertAfter(
              $(@element[0]).dialog("widget")
              .find(".ui-dialog-titlebar-"+name)
            )
          .end()
