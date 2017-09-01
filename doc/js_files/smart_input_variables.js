function TabLetterLinkDeactive() {
  // Get all elements with class="TabLetterLink" and remove the class "active"
  TabLetterLink = document.getElementsByClassName("TabLetterLink");
  for (i = 0; i < TabLetterLink.length; i++) {
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
