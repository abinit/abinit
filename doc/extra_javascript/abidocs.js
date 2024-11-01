
// Code executed on page ready
$(function() {

    // Activate popovers with tippyjs: https://atomiks.github.io/tippyjs/
    tippy('[data-toggle="popover"]', {
        position:'right',
    	animation:'scale',
        allowHTML: true,
        duration: [200, 100],
    	arrow:'true'
    });

    // This for the floating action button: https://bootsnipp.com/snippets/featured/inbox-by-gmail
    $('.fab').hover(function () {
        $(this).toggleClass('active');
    });

    // When arrow is clicked scroll to top of body
    $('#return-to-top').click(function() {      
        $('body,html').animate({scrollTop: 0}, 500);
    });

    // Add pip/unpin button to jquery dialog.
    // Based on http://appdevonsharepoint.com/how-to-pin-a-jquery-ui-dialog-in-place/
    $('button.PinDialog').click(function () {
        var CurrentDialogPosition = $(this).closest('.ui-dialog').offset();
        var DialogLeft = CurrentDialogPosition.left - $(window).scrollLeft();
        var DialogTop = CurrentDialogPosition.top - $(window).scrollTop();
        $(this).toggleClass('PinDialog DialogPinned').toggleClass('ui-state-highlight ui-state-default').children().toggleClass('ui-icon-pin-w ui-icon-pin-s').closest('.ui-dialog').css({ 'position': 'fixed', 'top': DialogTop, 'left': DialogLeft });
    });

    $('button.DialogPinned').click(function () {
        $(this).toggleClass('PinDialog DialogPinned').toggleClass('ui-state-highlight ui-state-default').children().toggleClass('ui-icon-pin-s ui-icon-pin-w').closest('.ui-dialog').css({ 'position': 'absolute', 'top': $(this).closest('.ui-dialog').offset().top, 'left': $(this).closest('.ui-dialog').offset().left });
    });

    // That selector matches all spans that have an id attribute and it starts with foo (e.g. fooblah
    $(".editor").each(function(index, element) {
        element.removeAttribute("hidden")
        var editor = ace.edit(element.id);
        //editor.setTheme("ace/theme/monokai");
        editor.setTheme("ace/theme/github");
        //editor.getSession().setMode("ace/mode/javascript");
        editor.setHighlightActiveLine(true);
        editor.setAutoScrollEditorIntoView(true);
        editor.setReadOnly(true);  // false to make it editable
    });

    // Allow users to edit markdown files on github. Based on
    // https://github.com/mkdocs/mkdocs/wiki/MkDocs-Recipes#associate-github-page-with-current-mkdoc-page
    // with minor modifications to select github social link.
    var git = 'https://github.com/abinit/abinit/edit/master/docs';
    var t1 = window.location.pathname;
    var url = null;
    if (t1 == '/') {
        url = git + '/index.md';
    } else {
        url = git + t1.substr(0, t1.length-1) + '.md';
    }   
    $("a.md-footer-social__link.fa-github-alt").attr('href', url).attr('target', '_blank');
});


// This for the table of variables with search bar implemented by Jordan

function TabLetterLinkDeactive() {
  // Get all elements with class="TabLetterLink" and remove the class "active"
  TabLetterLink = document.getElementsByClassName("TabLetterLink");
  for (var i = 0; i < TabLetterLink.length; i++) {
    TabLetterLink[i].className = TabLetterLink[i].className.replace(" active", "");
  }
}


function openLetter(evt, letter) {
  // Declare all variables
  var i, TabContentLetter, TabLetterLink;

  // Get all elements with class="TabContentLetter" and hide them
  TabContentLetter = document.getElementsByClassName("TabContentLetter");
  for (i = 0; i < TabContentLetter.length; i++) {
    TabContentLetter[i].style.display = "none";
  }

  TabLetterLinkDeactive();

  // Show the current tab, and add an "active" class to the button that opened the tab
  var myLetter = document.getElementById(letter)
  myLetter.style.display = "block";
  myLetter.getElementsByClassName('HeaderLetter')[0].style.display = "none";
  evt.currentTarget.className += " active";
} 


function searchInput() {
  // Declare variables
  var input, filter, ul, li, a, i, allLetters, letter;
  var TabContentLetter;

  TabLetterLinkDeactive();

  TabContentLetter = document.getElementsByClassName("TabContentLetter");
  for (i = 0; i < TabContentLetter.length; i++) {
    TabContentLetter[i].style.display = "none";
    TabContentLetter[i].getElementsByClassName("HeaderLetter")[0].style.display = "none";
  }

  input = document.getElementById('InputSearch');
  filter = input.value.toUpperCase();
  var ulLetters = document.getElementById('Letters');
  allLetters = ulLetters.getElementsByTagName('li');
  for (letter = 0; letter < allLetters.length; letter++) {
    liul = allLetters[letter];
    li = liul.getElementsByTagName('li');

    // Loop through all list items, and hide those who don't match the search query
    // Start at 1 to avoid finding the letters themselve
    for (i = 1; i < li.length; i++) {
      a = li[i].getElementsByTagName("a")[0];
      if (a.innerHTML.toUpperCase().indexOf(filter) > -1) {
        liul.getElementsByTagName('ul')[0].style.display = "block";
        li[0].style.display = "block";
        li[i].style.display = "block";
      } else {
        li[i].style.display = "none";
      }
    }
  }
}


function defaultClick(first) {
  if ( !first) {
    document.getElementById('InputSearch').value = '';
    searchInput();
  }
  if ( window.location.hash && ( /^#[a-zA-Z]$/.test(window.location.hash) ) ) {
    document.getElementById('click'+window.location.hash[1]).click();
  }
  else {
    document.getElementById("clickA").click();
  }
}

// Build customized jquery dialog.
// See http://api.jqueryui.com/dialog/ and https://github.com/ROMB/jquery-dialogextend
// The code to pin the dialog is based on http://appdevonsharepoint.com/how-to-pin-a-jquery-ui-dialog-in-place/

function abidocs_jqueryui_dialog(dialog_id, btn_id) { 

    //return;
    var e = $(dialog_id);
    e.dialog({
        width: 500,
        height: 500,
        //width: "auto",
        //height: "auto",
        autoOpen: false,
        show: {effect: "blind", duration: 200},
        //hide: {{effect: "explode", duration: 1000}},
        position: {my: 'center', at: 'center', of: window},
        buttons: [
            {text: "Close",
              icon: "ui-icon-close",
              classes: {"ui-button": "ui-corner-all"},
              click: function(){ $(this).dialog("close"); }
              //showText: false
            }
        ],
       create: function () {
         var titlebar = e.parent().children('.ui-dialog-titlebar');
         titlebar.prepend('<button class="ui-button ui-widget ui-state-default ui-corner-all ui-button-icon-only PinDialog" role="button" aria-disabled="false" title="pin down"><span class="ui-button-icon-primary ui-icon ui-icon-pin-w"></span></button>');

      }
   });

   // This does not work anylonger as it raises `dialogExtend is not a function`
   //e.dialogExtend({
   //    "maximizable": false, "minimizable": true, "collapsable": true, "minimizeLocation": "left",
   //    "dblclick": "collapse",
   //    "icons": {
   //         "close": "ui-icon-close",
   //         //"maximize": "ui-icon-extlink",
   //         "minimize": "ui-icon-minus",
   //         "restore": "ui-icon-newwin",
   //         "collapse": "ui-icon-triangle-1-s"
   //    }
   //});

   $(btn_id).click(function() { e.removeAttr('hidden').dialog('open'); });
}
