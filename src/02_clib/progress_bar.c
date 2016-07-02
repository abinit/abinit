/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: varia.c 2259 2006-06-29 18:29:48Z acastro $
*/

#include "abi_clib.h"
#include "xmalloc.h"

#include <sys/time.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <termios.h>
#include <unistd.h>  /* tcsetpgrp */

/* returns true if process is in the foreground, copied from openssh scp source */
static int foreground_proc(void)
{
	int ctty_pgrp;
	static pid_t pgrp = -1;
	
	if (pgrp == -1) pgrp = getpgrp();
	
#if defined HAVE_TCGETPGRP
	return ((ctty_pgrp = tcgetpgrp(STDOUT_FILENO)) != -1 && ctty_pgrp == pgrp);
#else 
	return ((ioctl(STDOUT_FILENO, TIOCGPGRP, &ctty_pgrp) != -1 && ctty_pgrp == pgrp));
#endif
}

/* Returns the width of the terminal (set to 80 if both winsize and $COLUMNS fail */
int get_tty_width(void)
{
  char* env;
  int cols=0;
	struct winsize winsize;
	
	if (ioctl(fileno(stdout), TIOCGWINSZ, &winsize) != -1) cols = winsize.ws_col;

  if (!cols) { /* Try $COLUMNS */                   
    env = getenv("COLUMNS");   
    if (env) cols = atoi(env);        
  }                              

  return ( cols>0 ? cols : 80);
}

/* Displays progress bar with a percentage */
void 
FC_FUNC_(clib_progress_bar,CLIB_PROGRESS_BAR)
  (const int *actual, const int *max)
{
	static struct timeval start;
	struct timeval now;
#define BUFFER_SIZE  512
#define FMT_SIZE     64
	char buf[BUFFER_SIZE], fmt[FMT_SIZE];
	int i, ratio, barlength, remaining;
	double elapsed;

  int my_actual = *actual;
  
	if (my_actual < 0) { /* initialize the timer */
		(void) gettimeofday(&start, (struct timezone *) 0);
		my_actual = 0;
    /* return; */
  }

#if 0
	if (!foreground_proc()) return;
#endif
	
	if (*max > 0) {
		ratio = 100 * my_actual / *max;
		if (ratio < 0  ) ratio = 0  ;
		if (ratio > 100) ratio = 100;
  }
  else
		ratio = 100;

	sprintf(buf, "%d", *max);
	i = strlen(buf);
	if (i<3) i=3;
	sprintf(fmt, "\r[%%%dd/%%%dd]", i, i);

  sprintf(buf, fmt, my_actual, *max);
  sprintf(buf + strlen(buf), " %3d%%" , ratio);

	barlength = get_tty_width() - strlen(buf) - 15;
	if (barlength > 0) {
		i = barlength * ratio / 100;
    sprintf(buf + strlen(buf),
						 "|%.*s%*s|", i,
						 "*******************************************************"
						 "*******************************************************"
						 "*******************************************************"
						 "*******************************************************"
						 "*******************************************************"
						 "*******************************************************"
						 "*******************************************************",
						 barlength - i, "");
	}

	/* time information now */
	(void) gettimeofday(&now, (struct timezone *) 0);
	elapsed = now.tv_sec - start.tv_sec;
	
	if (elapsed <= 0.0 || my_actual <= 0) {
    sprintf(buf + strlen(buf),"     --:-- ETA");
  }
  else {
		remaining = (int) (*max / (my_actual / elapsed) - elapsed);
		if (remaining < 0) remaining = 0;
		i = remaining / 3600;
		if (i)
      sprintf(buf + strlen(buf),"%4d:", i);
		else
      sprintf(buf + strlen(buf),"     ");

		i = remaining % 3600;
    sprintf(buf + strlen(buf),"%02d:%02d%s", i / 60, i % 60," ETA");
  }

	printf("%s", buf);
	fflush(stdout);
}

char *append(char *dest, size_t *size_dest, const char *src);

/* Append SRC to DEST. Return a pointer to DEST  */
char *append(char *dest, size_t *size_dest, const char *src) 
{
  size_t space_left = *size_dest - strlen(dest) -1;

  if ( strlen(src) > space_left ) { /* Have to realloc DEST */
    *size_dest = strlen(dest) + strlen(src) + 1;
    xrealloc(dest, *size_dest);
  }

  return strcat(dest, src);
}
