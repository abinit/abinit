---
authors: MG, XG
---

# HowTo use git and the gitlab server

This page is intended as a quick reference to *git* and its integration
with the Abinit project.
If you are not familiar with *git*, we would strongly advise to watch this tutorial:

<iframe width="1384" height="629" src="https://www.youtube.com/embed/HVsySz-h9r4" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

For further information about *git*, please consult 
the [official documentation](https://git-scm.com/).

!!! tip

    To access the online help for COMMAND, use:

        git COMMAND --help

    See also this page with the [most commonly used git tips and tricks](https://github.com/git-tips/tips)

In the next sections, we explain how to configure *git* to interoperate with the 
our [ABINIT gitlab server](https://gitlab.abinit.org/).
It is assumed you you already have an account on our **internal** gitlab server 
Note that having an account on gitlab.com or github.com is not enough since we run our own instance 
of the gilab server. You need to contact Jean-Michel to have an account created for your.
    

## Initial configuration

If this is the very first time you use *git*, please set the following global parameters before doing anything else:

```sh
mkdir -p $HOME/.config/git
git config --global user.name "Firstname Lastname"
git config --global user.email "someone@someserver.somedomain"
git config --global core.editor "my_preferred_editor"
git config --global core.excludesFile "$HOME/.config/git/ignore"
git config --global color.ui "auto"
git config --global merge.conflictstyle "diff3"
```

Replace *Firstname*, *Lastname*, *someone@someserver.somedomain*, and *my_preferred_editor*,
by your respective first name, last name, email address, and preferred editor.

To be able to **push your contributions** to the Abinit Forge, you need to add the
following section to your *~/.ssh/config* file

```
Host abinit-forge
HostName gitlab.abinit.org
    User git
    ServerAliveInterval 52
    Compression yes
```

so that one can use the **abinit-forge** hostname to **clone**, **pull**, and **push** to your repository.

For further info, please consult the [official documention](https://gitlab.abinit.org/help/ssh/README.md)

To clone your repository, execute:

    git clone abinit-forge:DEVELOPER/abinit

where DEVELOPER must be replaced by your Abinit Forge login.

In order to avoid typing your password every time you issue a command that accesses gitlab,
you have to introduce your public keys in your profile.
See <https://gitlab.abinit.org/profile/keys>.
On your local machine, generate a RSA ssh key **WITHOUT** passphrase:

    ssh-keygen -t rsa

and call it *id_rsa_gitlab*. Then add a section in the *~/.ssh/config file*:

```
host gitlab
  Hostname gitlab.abinit.org
  User git
  KeepAlive yes
  IdentityFile ~/.ssh/id_rsa_gitlab
```

Finally, copy the public key *id_rsa_gitlab.pub* on gitlab.
Now, you can use (on your local machine) the following syntax:

    git clone gitlab:user_id/abinit.git

instead of:

    git clone git@gitlab.abinit.org:user_id/abinit.git

To be sure the key is proposed each time git calls ssh, you can use ssh-agent:

    ssh-agent # this starts the agent, and provides the process id
    # execute the 3 lines of commands that ssh-agent proposes, e.g.
    SSH_AUTH_SOCK=/tmp/ssh-ngsERHER3K1HS/agent.15589; export SSH_AUTH_SOCK;
    SSH_AGENT_PID=15590; export SSH_AGENT_PID;
    echo Agent pid 15590;
    ssh add ~/.ssh/id_rsa_gitlab # add the corresponding ssh key for gitlab


Every developer has his/her specific gitlab **user_id** e.g. *gonze*.
An additional *virtual* developer, called **trunk**, is also defined.
<!--
===== Standard names (projects, ID, branches) =====
On the ABINIT gitlab server, an ABINIT project can be defined, specifically for you,
with address <color blue>git@gitlab.abinit.org</color>:<color red>user_id</color><color blue>/abinit.git</color>.
In order to do this, log on https://gitlab.abinit.org, go to "Explore projects", find "Trunk / Abinit" (+perhaps select "All"),
then click on the "Fork" button. You need to contact Jean-Michel Beuken to have access to the ABINIT gitlab.
In order to start to work locally (not on the ABINIT gitlab server, but on your own machine),
you should setup the SSH environment for gitlab, as described
In particular on some machines you need to have an ssh agent running with your rsa key already available,
so that git finds it when it runs ssh.
-->
You have by default a **master** branch and a **develop** branch in your repo,
but it is also possible to create and work in other branches.

To clone from the gitlab repository to your local repository, the recommended command is:

    git clone gitlab:user_id/abinit.git -b develop

where **user_id** is to be replaced by the adequate name.
Additional explanation: the normal *clone* command from the gitlab repository
to his/her gitlab repository on his/her local machine is:

    git clone git@gitlab.abinit.org:user_id/abinit.git

However, with the proper SSH environment, the following *clone* command can be used:

    git clone gitlab:user_id/abinit.git

This will clone the *master*, while the usual working branch is the
*develop* branch, that can be obtained directly with the above-mentioned command.
After some modifications, using the git commands (e.g. git add, git commit ...),
you can push on the gitlab server thanks to:

    git push  # Add --tags to push tags.

In order for the modifications to be merged in the trunk, a merge request (MR) has to be issued, as described later.
The same might be done for other branches, created from the local *develop* one, using:

    git branch another_branch
    git checkout another_branch

For pushing, the first time, use:

    git push -u origin another_branch

You are able to create your own additional branches, either locally or on gitlab.
The name is your own choice.


!!! important

    Note that, according to the [gitflow branching model](http://nvie.com/posts/a-successful-git-branching-model),
    the **master** branch will be of little use for the *physical* developers (including you)
    while the **develop** and **release-*** branches will be quite important for normal developments.
    In a nutshell, new developments are always merged in **trunk/develop** while **trunk/master** contains
    the stable version. Finally, **release-** branches are mainly used for bug fixes and minor changes
    that should be included in the next release.

    ![](https://miro.medium.com/max/1400/1*9yJY7fyscWFUVRqnx0BM6A.png)

### Git branches

In the ABINIT + git workflow:

  * In order to be merged in the trunk, the branch has to be **on-track**, as explained later
  * You have to issue a "merge request" on gitlab (usually to **trunk/develop**).
  * By default, all branches pushed on gitlab might be tested on-demand by the user.
  * Only branches that are **on-track** and that have succeeded (or passed) all the automatic tests
    are considered for merge (as was the case previously).

The *develop* branch corresponds to the development stage.
The *release-* branches corresponds to release candidates or hotfix branches.
It is slightly simpler than the gitflow naming scheme.

!!! tip

    In the git philosophy, branches are often created, then merged, then destroyed.
    Developers are encouraged to take advantage of this flexibility.

### How to trigger the action of buildbot? (On-demand)

The developer has to use the on-demand form on the ABINIT [buildbot portal](https://bbportal.abinit.org)
The authentification is made with the gitlab login/password.
One important feature has to be taken into account: with git, each commit is characterized
by a SHA-1 code (40 hexadecimal characters).
The four first characters actually define the commit with 1 combination out of 65536 4-character combinations,
while the five first characters define the commit with 1 combination out of 1048576 5-character combinations.
Gitlab knows the repository of the different users, the existing branches, as well as the SHA-1 code of each commit.
The three following fields must be defined by the developer for buildbot to know which is the revision of ABINIT to be run:

  - “User” : choose one of the existing User_ID from the menu. Default value: *user_id*.
  - “Branch” : choose one of the existing User_ID from the menu. Default value: *develop*.
  - “Commit” : the developer can enter a free text with the first hexadecimal characters
     of the SHA-1 code of the commit (in case different commits of the user on the branch
     have the same first characters, the last commit will be taken by buildbot).
     Default value: the last commit of User/Branch.

### What is an on-track branch?

For a branch to be **on-track**, a specific commit must be contained in the history of the branch.
This specific commit will be created in *trunk/develop* and corresponds to the release
of a new version (not necessarly a public version).
For a *release-* branch to be **on-track**, not only a specific commit must
be contained in the past history of the branch, but also the commit that starts
the next development version **cannot** be present.
Only one among the *release-* branches is **on-track** at any time.

!!! note

As an example, suppose we were on ABINIT v8.9.x, and want to start preparing a release 8.10.0 (for production)
and a new v8.11.0 (for development):

  * a branch entitled *release-8.10.0* will be forked from the *develop* branch
  * after this branching, the first commit in the *release-8.10.0* branch will be tagged start-8.10.0,
    while the first commit in the *develop* branch will be tagged start-8.11.0,
  * for a *develop* branch to be considered **on-track** (and eligible for a merge request), the commit
    tagged start-8.11.0 will have to be present;
  * for a *release-8.10.0* branch to be considered **on-track** (and eligible for a merge request),
    the commit tagged start-8.10.0 will have to be present, while it will be forbidden to contain
    the commit tagged start-8.11.0 .

In complement to the start of a X.Y.Z version being tagged as "start-X.Y.Z",
the commit that ends some X.Y.Z version of ABINIT will be tagged "X.Y.Z".

### What is shown in the Buildbot Status table?

The [Buildbot Status table](https://bbportal.abinit.org/#/status)
shows a summary of the recent results obtained by the test farm for alll active branches.
If you want to see all the results, select the **filters**: all + all, then click on the **update button**.
The default selection of bots, that is the *night* set, should be adequate for all merge requests.
The selection capabilities of the Buildbot Status table are rather extended, you have to play a bit wit them.

By the way, the behaviour of this Buildbot Status table when changing the selections has still some problems -as of May 2018-.
Do not hesitate to click on the update button, and the "Grouping" button -the latter back and forth-, sorry for the inconvenience.

### How and when will the merge in the master branch be done?

In order for your development to be incorporated, there should be **merge request**
from your specific branch to the *trunk*.
For most branches, the merge request will be to the *trunk/develop* branch.
The list of merge requests is available [here](https://bbportal.abinit.org/#/mr).

However, when a *release-* branch is ready to be merged,
the merge request should target the corresponding *trunk/release-*</color> branch.
The master branch is only used by the trunk.
So, never issue a merge request to trunk/master.
You are supposed to merge, inside your branches, specific tagged branches from trunk,
in order to avoid divergences. These specific tagged branches will be advertised.
As mentioned earlier,
the presence of such tagged commits allows one to identify whether the branch is **on-track**,
before an automatic testing to be done.

### How to synchronize with the trunk?

In order to keep your branches up to date with those
of the trunk, you should first register *git@gitlab.abinit.org:trunk/abinit.git* as a remote with:

    git remote add trunk gitlab:trunk/abinit.git

At this point, one can fetch the branches in trunk with:

    git fetch trunk

then, if the develop branch is to be updated, supposing it is checked out,
To merge *trunk/develop* in your develop branch:

    git checkout develop
    git merge trunk develop

You can combine the last two commands in one as:

    git pull trunk develop

If, on the contrary, a new branch (e.g. a release branch, let's says 8.8 to fix the ideas) has to be created:

    git branch release-8.8 start-8.8.1    # this creates the branch release-8.8 from the start-8.8.1 tag
    git checkout release-8.8
    git merge remotes/trunk/release-8.8
    git push -u origin release-8.8

That's it! You can now make modifications in your release-8.8, then issue a merge request to the trunk/release-8.8.

<!--
### Additional info: Setup of the SSH environment
## How to clone your repository with git and track `trunk`

To clone your repository on your `localhost`:

    git clone

To track `trunk`:

    git remote add trunk

To show the list of remote branches:

    git remote -v

To merge the `develop` branch of `trunk` in your branch:

    git checkout develop
    git pull trunk develop

To push to the gilab server:

    git push origin develop


!!! tip

    To access the online help use: `git COMMAND --help`


!!! important

    gitflow: You should always send pull requests to trunk/develop
-->
