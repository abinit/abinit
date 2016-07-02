/*{\src2tex{textfont=tt}}
 ****p* ABINIT/Timeout/timeout
 * NAME
 *  timeout
 *
 * FUNCTION
 *  Executes a command and imposes an elapsed time limit. The command
 *  is run in a separate POSIX process group so that the right thing
 *  happens with commands that spawn child processes.
 *
 * COPYRIGHT
 *  This program is part of SATAN.
 *  Copyright 1995 by Dan Farmer and Wietse Venema. All rights reserved.
 *
 *  Redistribution and use in source and binary forms are permitted
 *  provided that this entire copyright notice is duplicated in all
 *  such copies. No charge, other than an "at-cost" distribution fee,
 *  may be charged for copies, derivations, or distributions of this
 *  material without the express written consent of the copyright
 *  holders.
 *
 *  THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 *  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 *  WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR ANY PARTICULAR
 *  PURPOSE.
 *
 *  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, LOSS OF USE, DATA, OR
 *  PROFITS OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 *  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 *  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 *  DAMAGE.
 *
 * NOTES
 *  Command line: timeout [-<signal>] time command
 *
 *    -<signal>
 *      Specify an optional signal to send to the controlled process.
 *      By default, timeout sends SIGKILL, which cannot be caught or
 *      ignored.
 *
 *    time
 *      The elapsed time limit after which the command is terminated.
 *
 *    command
 *      The command to be executed.
 *
 *    The exit status is the one of the command (status 1 in case of a
 *    usage error).
 *
 * SOURCE
 */

/* System libraries. */

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>

extern int optind;

/* Application-specific. */

#define perrorexit(s) { perror(s); exit(1); }

static int kill_signal = SIGKILL;
static char *progname;
static char *commandname;
static int  time_to_run;

static void usage()
{
    fprintf(stderr, "usage: %s [-signal] time command...\n", progname);
    exit(1);
}

static void terminate(sig)
int     sig;
{
    signal(kill_signal, SIG_DFL);
    fprintf(stderr, "timeout.c: aborting command ``%s'' after %d seconds with signal %d\n",
	    commandname, time_to_run, kill_signal);
    kill(0, kill_signal);
}

int     main(argc, argv)
int     argc;
char  **argv;
{
    /*int     time_to_run;*/
    pid_t   pid;
    pid_t   child_pid;
    int     status;

    progname = argv[0];

    /*
     * Parse JCL.
     */
    while (--argc && *++argv && **argv == '-')
	if ((kill_signal = atoi(*argv + 1)) <= 0)
	    usage();

    if (argc < 2 || (time_to_run = atoi(argv[0])) <= 0)
	usage();

    commandname = argv[1];

    /*
     * Run the command and its watchdog in a separate process group so that
     * both can be killed off with one signal.
     */
    setsid();
    switch (child_pid = fork()) {
    case -1:					/* error */
	perrorexit("timeout: fork");
    case 0:					/* run controlled command */
	execvp(argv[1], argv + 1);
	perrorexit(argv[1]);
    default:					/* become watchdog */
	(void) signal(SIGHUP, terminate);
	(void) signal(SIGINT, terminate);
	(void) signal(SIGQUIT, terminate);
	(void) signal(SIGTERM, terminate);
	(void) signal(SIGALRM, terminate);
	alarm(time_to_run);
	while ((pid = wait(&status)) != -1 && pid != child_pid)
	     /* void */ ;
	return (pid == child_pid ? WEXITSTATUS(status) : -1);
    }
}
/*****/
