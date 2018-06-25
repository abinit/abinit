jquery-dialogextend 2.0.4 [![project status](http://stillmaintained.com/ROMB/jquery-dialogextend.png)](http://stillmaintained.com/ROMB/jquery-dialogextend) [![Build Status](https://travis-ci.org/ROMB/jquery-dialogextend.png?branch=master)](https://travis-ci.org/ROMB/jquery-dialogextend)
===
Download
===
[development build](https://raw.github.com/ROMB/jquery-dialogextend/master/build/jquery.dialogextend.js)

[minified library](https://raw.github.com/ROMB/jquery-dialogextend/master/build/jquery.dialogextend.min.js)

Compatible
===
- jQuery 1.11.1
- jQueryUI 1.11.0

Overview
===
- Neat, simple, and ABSOLUTELY unobtrusive
- Extending (instead of replacing) original jQuery UI dialog
- Maximize and minimize buttons
- Show/Hide close button
- Double-clickable title bar
- Enhanced title bar options
- Configurable icons
- Custom events

Demo
===
- Test Tool : [http://romb.github.io/jquery-dialogextend/example.html](http://romb.github.io/jquery-dialogextend/example.html)

Tested Browsers
===
- Chrome 35
- Firefox 14
- IE 8

Please support this project
===

Donate Bitcoins: 1G8T7Xh2AN5ceduHmHT5TpPFUeddsnQHLQ

Options
===

#### closable ####
Type: *Boolean*

Usage: enable/disable close button

Default: *true*


#### maximizable ####
Type: *Boolean*

Usage: enable/disable maximize button

Default: *false*

#### minimizable ####

Type: *Boolean*

Usage: enable/disable minimize button

Default: *false*

#### collapsable ####
Type: *Boolean*

Usage: enable/disable collapse button

Default: *false*

#### minimizeLocation ####

Type: *String*

Usage: sets alignment of minimized dialogues

Default: *'left'*

Valid: *'left'*, *'right'*

#### dblclick ####

Type: *Boolean*, *String*

Usage: set action on double click

Default: *false*

Valid: *false*, *'maximize'*, *'minimize'*, *'collapse'*


#### titlebar ####

Type: *Boolean*, *String*

Default: *false*

Valid: *false*, *'none'*, *'transparent'*


#### icons ####

Type: *Object*

Default:

    {
      "close" : "ui-icon-circle-closethick", // new in v1.0.1
      "maximize" : "ui-icon-extlink",
      "minimize" : "ui-icon-minus",
      "restore" : "ui-icon-newwin"
    }

Valid: *&lt;jQuery UI icon class&gt;*

Events
===

#### load ####

Type: *load*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "load" : function(evt, dlg) { ... }
	});
	//Bind to event by type
	//NOTE : You must bind() the <dialogextendload> event before dialog-extend is created
	$("#my-dialog")
	  .bind("dialogextendload", function(evt) { ... })
	  .dialogExtend();


#### beforeCollapse ####

Type: *beforeCollapse*

Example:

    //Specify callback as init option
    $("#my-dialog").dialogExtend({
      "beforeCollapse" : function(evt) { ... }
    });
    //Bind to event by type
    $("#my-dialog").bind("dialogextendbeforeCollapse", function(evt) { ... });

#### beforeMaximize ####

Type: *beforeMaximize*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "beforeMaximize" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendbeforeMaximize", function(evt) { ... });

#### beforeMinimize ####

Type: *beforeMinimize*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "beforeMinimize" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendbeforeMinimize", function(evt) { ... });

#### beforeRestore ####

Type: *beforeRestore*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "beforeRestore" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendbeforeRestore", function(evt) { ... });

#### collapse ####

Type: *collapse*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "collapse" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendcollapse", function(evt) { ... });

#### maximize ####

Type: *maximize*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "maximize" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendmaximize", function(evt) { ... });

#### minimize ####

Type: *minimize*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "minimize" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendminimize", function(evt) { ... });

#### restore ####

Type: *restore*

Example:

	//Specify callback as init option
	$("#my-dialog").dialogExtend({
	  "restore" : function(evt) { ... }
	});
	//Bind to event by type
	$("#my-dialog").bind("dialogextendrestore", function(evt) { ... });

Methods
===
#### collapse ####

Usage: Collapse the dialog without double-clicking the title bar

Trigger: *dialogextendbeforeCollapse*, *dialogextendcollapse*

Example:

	$("#my-dialog").dialogExtend("collapse");
#### maximize ####

Usage: Maximize the dialog without clicking the button

Trigger: *dialogextendbeforeMaximize*, *dialogextendmaximize*

Example:

	$("#my-dialog").dialogExtend("maximize");

#### minimize ####

Usage: Minimize the dialog without clicking the button

Trigger: *dialogextendbeforeMinimize*, *dialogextendminimize*

Example:

	$("#my-dialog").dialogExtend("minimize");

#### restore ####

Usage: Restore the dialog from maximized/minimized/collapsed state without clicking the button

Trigger: *dialogextendbeforeRestore*, *dialogextendrestore*

Example:

	$("#my-dialog").dialogExtend("restore");

#### state ####

Usage: Get current state of dialog

Return: *String*

Value: *'normal'*, *'maximized'*, *'minimized'*, *'collapsed'*

Example:

	switch ( $("#my-dialog").dialogExtend("state") ) {
	  case "maximized":
	    alert("The dialog is maximized");
	    break;
	  case "minimized":
	    alert("The dialog is minimized");
	    break;
	  case "collapsed":
	    alert("The dialog is collapsed");
	    break;
	  default:
	    alert("The dialog is normal");
	}

Theming
===
The dialog will have class according to its current state.

	<div class="ui-dialog">
	  <div class="ui-dialog-titlebar">...</div>
	  <div class="ui-dialog-content ui-dialog-{normal|maximized|minimized|collapsed}">...</div>
	</div>
The buttons are wrapped by title bar of jQuery UI Dialog.

*Note : After using dialogExtend, close button will not be a direct child of title bar anymore. It will be wrapped by a button pane element*

	<div class="ui-dialog-titlebar ui-widget-header ui-corner-all ui-helper-clearfix">
	  ...
	  <div class="ui-dialog-titlebar-buttonpane">
	    <a class="ui-dialog-titlebar-close ui-corner-all" href="#">...</a>
	    <a class="ui-dialog-titlebar-maximize ui-corner-all" href="#"><span class="ui-icon {icons.maximize}">maximize</span></a>
	    <a class="ui-dialog-titlebar-minimize ui-corner-all" href="#"><span class="ui-icon {icons.minimize}">minimize</span></a>
	    <a class="ui-dialog-titlebar-restore ui-corner-all" href="#"><span class="ui-icon {icons.restore}">restore</span></a>
	  </div>
	  ...
	</div>

Example - Basic Config
===
	$(function(){
	  $("#my-button").click(function(){
	    $("<div>This is content</div>")
	      .dialog({ "title" : "My Dialog" })
	      .dialogExtend({
	        "maximizable" : true,
	        "dblclick" : "maximize",
	        "icons" : { "maximize" : "ui-icon-arrow-4-diag" }
	      });
	  });
	});
Example - Full Config
===
	$(function(){
	  $("#my-button").click(function(){
	    $("<div>This is  content</div>")
	      .dialog({
	        "title" : "This is dialog title",
	        "buttons" : { "OK" : function(){ $(this).dialog("close"); } }
	       })
	      .dialogExtend({
	        "closable" : true,
	        "maximizable" : true,
	        "minimizable" : true,
	        "collapsable" : true,
	        "dblclick" : "collapse",
	        "titlebar" : "transparent",
	        "minimizeLocation" : "right",
	        "icons" : {
	          "close" : "ui-icon-circle-close",
	          "maximize" : "ui-icon-circle-plus",
	          "minimize" : "ui-icon-circle-minus",
	          "collapse" : "ui-icon-triangle-1-s",
	          "restore" : "ui-icon-bullet"
	        },
	        "load" : function(evt, dlg){ alert(evt.type); },
	        "beforeCollapse" : function(evt, dlg){ alert(evt.type); },
	        "beforeMaximize" : function(evt, dlg){ alert(evt.type); },
	        "beforeMinimize" : function(evt, dlg){ alert(evt.type); },
	        "beforeRestore" : function(evt, dlg){ alert(evt.type); },
	        "collapse" : function(evt, dlg){ alert(evt.type); },
	        "maximize" : function(evt, dlg){ alert(evt.type); },
	        "minimize" : function(evt, dlg){ alert(evt.type); },
	        "restore" : function(evt, dlg){ alert(evt.type); }
	      });
	  });
	});
