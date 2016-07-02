/* Formats Mathematica output code after it has been run through the
   sed script to straighten out the continuation lines for free-format
   fortran 90 and "compactify" the code.  Note that 39 continuation
   lines are all that are allowed in the Fortran 95 standard, so more
   are flagged for hand-coding attention.  Lines up to 132 characters
   are allowed, but there is no sense in making the code even less
   readable than this already does. */

#include <stdio.h>

#define MAX_CHAR  65
#define MAX_CONTINUE  39

main()
{
  int c;
  int n_char,n_paren,n_continue;

  n_continue=0;

  printf("\n    ");
  n_char=4;

  while((c = getchar()) != EOF){
    if(c == '#'){
      printf("\n    ");
      n_char = 4;
      n_continue = 0;
      n_paren = 0;
    }
    else if(c == '\\'){
      if((c = getchar()) == '\n'){
        if((c = getchar()) == '#'){
          printf("\n    ");
          n_char = 4;
          n_continue = 0;
          n_paren = 0;
        }
      }
    }

    if(n_char >= MAX_CHAR && c != ')' && n_paren > 0){
      printf("&\n");
      n_continue += 1;
      if(n_continue >= MAX_CONTINUE){
        printf("!***** WARNING: continue lines > MAX_CONTINUE *****\n");
        n_continue = 0;
      }
      printf("&    ");
      n_char = 5;
      n_paren = 0;
    }

    if(c != ' ' && c != '\n' && c != '#'){
      putchar(c);
      n_char += 1;
    }

    if(n_char >= MAX_CHAR && c == ')') n_paren += 1;

  }
}
