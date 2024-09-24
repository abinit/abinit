# git(lab) : ABINIT specificities

## Standard names (projects, ID, branches)

Let us suppose that you are one of the ABINIT developers ...

Every "physical" developer of ABINIT using gitlab has his/her specific <span style="color:red">user_id</span> (e.g. <span style="color:red">gonze</span>).
Moreover, an additional "virtual" developer, called <span style="color:red">trunk</span>, is also defined.

On the ABINIT gitlab server, an ABINIT project can be defined, specifically for you, with address <span style="color:blue">git@gitlab.abinit.org</span>:<span style="color:red">user_id</span><span style="color:blue">/abinit.git</span>. In order to do this, log on [https://gitlab.abinit.org](https://gitlab.abinit.org), go to "Explore projects", find "Trunk / Abinit" (+perhaps select "All"), then click on the "Fork" button. You need to contact Jean-Michel Beuken to have access to the ABINIT gitlab.

In order to start to work locally (not on the ABINIT gitlab server, but on your own machine), you should setup the SSH environment for gitlab, as described  [in the last section of this document](#Additional_infosetup_of_the_ssh_environment). In particular on some machines you need to have an ssh agent running with your rsa key already available, so that git finds it when it runs ssh.

You have by default a <span style="color:magenta">master</span> branch and a <span style="color:magenta">develop</span> branch in your repository, but can create or work on other branches, following the philosophy of git. In particular, according to  [gitflow](http://nvie.com/posts/a-successful-git-branching-model), this <span style="color:magenta">master</span> branch will be of little use for the "physical" developers (including you), as we will see later, while the <span style="color:magenta">develop</span> and <span style="color:magenta">release-*</span> branches will be quite important for you.

The visibility level of the ABINIT projects of the <span style="color:red">user_id</span> developers and of the <span style="color:red">trunk</span> developer is "internal". The <span style="color:red">trunk</span><span style="color:blue">/abinit.git</span> (well, only the <span style="color:magenta">master</span> branch) is however mirrored on github.

## (Quick start) Cloning from gitlab

In order to clone from the gitlab repository to your gitlab repository on your local machine (either for the first time, or in order to restart from a fresh copy of the gitlab repository), the recommended command is : \\
<span style="color:blue">git clone gitlab:</span><span style="color:red">user_id</span><span style="color:blue">/abinit.git</span> -b <span style="color:magenta">develop</span>\\
where <span style="color:red">user_id</span> is to be replaced by the adequate name.\\
\\
[Additional explanation : the normal "clone" command from the gitlab repository to his/her gitlab repository on his/her local machine is : \\
<span style="color:blue">git clone git@gitlab.abinit.org:</span><span style="color:red">user_id</span><span style="color:blue">/abinit.git</span> \\
However, with the proper SSH environment, the following "clone" command can be used : \\
<span style="color:blue">git clone gitlab:</span><span style="color:red">user_id</span><span style="color:blue">/abinit.git</span> \\
This will clone the <span style="color:magenta">master</span>, while the usual working branch is the
<span style="color:magenta">develop</span> branch, that can be obtained directly with the above-mentioned command.]\\
\\
After some modifications, using the git commands (e.g. git add, git commit ...), you can push on the gitlab server thanks to : \\
<span style="color:blue">git push --tags  </span>     <span style="color:red">( double-dash before tags )</span>\\
\\
In order for the modifications to be merged in the trunk, a merge request has to be issued, as described later.\\
\\
The same might be done for other branches, created from the local <span style="color:magenta">develop</span> one, using the commands : \\ 
<span style="color:blue">git branch </span><span style="color:magenta">another_branch</span> \\
<span style="color:blue">git checkout </span><span style="color:magenta">another_branch</span> \\
For pushing, the first time, use : \\
 <span style="color:blue">git push -u origin</span> <span style="color:magenta">another_branch</span>.\\
You are able to create your own additional branches, either locally or on gitlab. The name is your own choice.

## Git branches

<span style="color:brown">In the ABINIT+git workflow</span>: 
  * By default, all branches pushed on gitlab might be tested on-demand by the user.
  * In order to be merged in the trunk, the branch has to be **on-track**, [as explained later](#what_is_an_on-track_branch). You have to issue a "merge request" on gitlab (to the <span style="color:red">trunk</span> virtual developer). 
  * Only branches that are **on-track** and that have succeeded (or passed) all the automatic tests are considered for merge (as was the case previously).

The <span style="color:magenta">develop</span> branch in the ABINIT naming scheme corresponds to the development stage. The <span style="color:magenta">release-*</span> branches in the ABINIT naming scheme  corresponds to release candidates or hotfix branches. It is slightly simpler than the gitflow naming scheme.

In the git philosophy, branches are often created, then merged, then destroyed. The ABINIT developer is encouraged to take advantage of this flexibility. However, many developers work only on one main project, and will already be satisfied by simply using their <span style="color:magenta">develop</span> branch.

## How to trigger the action of buildbot ? (On-demand)

The developer has to use the //git on-demand form on the ABINIT buildbot portal ([https://bbportal.abinit.org](https://bbportal.abinit.org))//. The authentification is made with login/password of gitlab.

One important feature has to be taken into account : with git, each commit is characterized by a SHA-1 code (40 hexadecimal characters). The four first characters actually define the commit with 1 combination out of 65536 4-character combinations, while the five first characters define the commit with 1 combination out of 1048576 5-character combinations.

Gitlab knows the repository of the different users, the existing branches, as well as the SHA-1 code of each commit. The three following fields must be defined by the developer for buildbot to know which is the revision of ABINIT to be run:
  - “User” : choose one of the existing User_ID from the menu. Default value : <span style="color:red">user_id</span> .
  - “Branch” : choose one of the existing User_ID from the menu. Default value : <span style="color:magenta">develop</span>.
  - “Commit” : the developer can enter a free text with the first hexadecimal characters of the SHA-1 code of the commit (in case different commits of the user on the branch have the same first characters, the last commit will be taken by buildbot). Default value : the last commit of User/Branch.

## What is an "on-track" branch ?

For a <span style="color:magenta">develop</span> branch to be **on-track**, a specific commit must be contained in the past history of the branch. This specific commit will be tagged in the <span style="color:red">trunk</span>/<span style="color:magenta">develop</span> branch (actually, corresponding to the definition of the start of the version).

For a <span style="color:magenta">release-*</span> branch to be **on-track**, not only a specific commit must be contained in the past history of the branch, but also the commit that starts the next development version **cannot** be present. Only one among the <span style="color:magenta">release-*</span> branches is **on-track** at any time.

As an example, suppose we were on ABINIT v9.9.3, and want to start preparing a release 9.10.1 (for release) and a new v9.11.0 (for development)   :
  * a branch entitled <span style="color:magenta">release-9.10</span> will be forked from the <span style="color:magenta">develop</span> branch ;
  * after this branching, suppose the first commit in the <span style="color:magenta">release-9.10</span> branch is 1d5efe42, while the first commit in the <span style="color:magenta">develop</span> is the commit 34dbe7f2d
  * for a <span style="color:magenta">develop</span> branch to be considered **on-track** (and eligible for a merge request), the commit 34dbe7f2d must be present ;
  * for a <span style="color:magenta">release-9.10</span> branch to be considered **on-track** (and eligible for a merge request), the commit 1d5efe42 must be present, while it will be forbidden to contain the commit 34dbe7f2d .

In complement to the start of a X.Y.Z version for release being tagged as "start-X.Y.Z",
the commit that ends some X.Y.Z version for release will be tagged "X.Y.Z".


## What is shown in the Buildbot Status table ?

(See the Buildbot Status at [https://bbportal.abinit.org/#/status](https://bbportal.abinit.org/#/status))

The Buildbot Status table show a summary of the recent results obtained by the test farm for active branches.
If you want to see all the results, select the "filters" : all + all , then click on the update button.
The default selection of bots, that is the "night" set, should be adequate for all merge requests.
The selection capabilities of the Buildbot Status table are rather extended, you have to play a bit wit them.

## How and when will the merge in the master branch be done ?

In order for your development to be incorporated, there should be a gitlab “merge request” from your specific <span style="color:red">user_id</span> to the <span style="color:red">trunk</span>. 
For most branches of <span style="color:red">user_id</span>, the merge request will be to the <span style="color:red">trunk</span>/<span style="color:magenta">develop</span> branch. 
(See the list of merge requests at [https://bbportal.abinit.org/#/mr](https://bbportal.abinit.org/#/mr))

However, when a <span style="color:magenta">release-*</span> branch is ready to be merged, the "merge request" should target the corresponding <span style="color:red">trunk</span>/<span style="color:magenta">release-*</span> branch. 

The master branch is only used by the trunk. So, never issue a merge request to trunk/master.

You are supposed to merge, inside your branches, specific tagged branches from the <span style="color:red">trunk</span> user, in order to avoid divergences. These specific tagged branches will be advertised. As [mentioned earlier](#what_is_an_on-track_branch), the presence of such tagged commits allows one to identify whether the branch is **on-track**, before an automatic testing to be done.

## How to synchronize with the "trunk" virtual user ?

In order to keep your branches <span style="color:red">user_id</span>/<span style="color:magenta">develop</span> or 
<span style="color:red">user_id</span>/<span style="color:magenta">release-*</span> up to date with those of the <span style="color:red">trunk</span> virtual user, you should first define git@gitlab.abinit.org:trunk/abinit.git as a remote :
<code>
   git remote add trunk gitlab:trunk/abinit.git
</code>
Then, in order to synchronize, in a first step, issue :
<code>
   git fetch --tags trunk
</code>
then, if the develop branch is to be updated, supposing it is checked out, merge the trunk/develop in your develop : 
<code>
   git merge remotes/trunk/develop
</code>

## How to set up a "hotfix" branch to be merged in trunk/release-x.y ?

If you have a hotfix, a new branch has to be created. To fix the ideas, let's suppose you want to communicate a bug fix to ABINITv9.8, you have to issue :
<code>
   git branch release-9.8 start-9.8.0    (this creates the branch release-9.8 from the start-9.8.0 tag)
   git checkout release-9.8
   git merge remotes/trunk/release-9.8
   git push -u origin release-9.8 --tags
</code>
That's it ! You can now make modifications in your release-9.8, then issue a merge request to the trunk/release-9.8 .

## Additional info: how to not mess with your branches ? ... show-branch ...

To see which branches are present in your current copy of abinit type:
git show-branch -a 
(documented at [https://git-scm.com/docs/git-show-branch](https://git-scm.com/docs/git-show-branch)).
This should also list the remotes (if not add -r). An example output:
<code>
* [develop] Merge branch 'develop' of gitlab:trunk/abinit into develop
 ! [master] Minor modif, to initialize v8.7.3 modified:   KNOWN_PROBLEMS
  ! [origin/HEAD] Minor modif, to initialize v8.7.3 modified:   KNOWN_PROBLEMS
   ! [origin/develop] Merge branch 'develop' of gitlab:trunk/abinit into develop
    ! [origin/master] Minor modif, to initialize v8.7.3 modified:   KNOWN_PROBLEMS
     ! [trunk/develop] Import inside gitlab develop, the modifications from v8.6.2 to v8.6.3 
      ! [trunk/master] Minor modif, to initialize v8.7.3 modified:   KNOWN_PROBLEMS
</code>

The syntax is compact but a bit barbaric: the first section gives you the branches which are actually present, with a star for the one currently checked out. The other branches are offset by 1 space.

<code>
-  -    [develop] Merge branch 'develop' of gitlab:trunk/abinit into develop
*  + +  [trunk/develop] Import inside gitlab develop, the modifications from v8.6.2 to v8.6.3 
</code>


The second section shows commits, and can tell you the relation between the branches. Commits are marked with a + in each column for each branch they are included in, and merges with a -. Note that the master branches are never used, the top commit is a merge of trunk into the current develop, and the second commit is contained in the currently active, the origin/develop (of which it is just a clone), and the trunk/develop branches.

## Additional info: Setup of the SSH environment

In order to avoid typing your password every time you issue a command that accesses gitlab,
you have to introduce your public keys in your profile ( [https://gitlab.abinit.org/profile/keys](https://gitlab.abinit.org/profile/keys) ).

On your local machine, generate a ssh key of <span style="color:red">RSA</span> type only <span style="color:red">WITHOUT passphrase</span> :
<code>
  ssh-keygen -t rsa
</code>  
and call it //id_rsa_gitlab// .

Then add a section in the  ~/.ssh/config file :

<code>
  host gitlab
   Hostname gitlab.abinit.org
   User git
   KeepAlive yes
   IdentityFile ~/.ssh/id_rsa_gitlab
</code>   

and then, copy the public key //id_rsa_gitlab.pub// on gitlab.

Now, you can use (on your local machine) the following syntax :\\
<span style="color:blue">git clone gitlab:</span><span style="color:red">user_id</span><span style="color:blue">/abinit.git</span> \\   
instead of the above-mentioned\\
<span style="color:blue">git clone git@gitlab.abinit.org:</span><span style="color:red">user_id</span><span style="color:blue">/abinit.git</span>

To be sure the key is proposed each time git calls ssh, you can use ssh-agent:
<code>
  ssh-agent # this starts the agent, and provides the process id
  execute the 3 lines of commands that ssh-agent proposes, e.g.
  SSH_AUTH_SOCK=/tmp/ssh-ngsERHER3K1HS/agent.15589; export SSH_AUTH_SOCK;
  SSH_AGENT_PID=15590; export SSH_AGENT_PID;
  echo Agent pid 15590;
  ssh add ~/.ssh/id_rsa_gitlab # add the corresponding ssh key for gitlab
</code>
For further details, please consult the official documentation [https://gitlab.abinit.org/help/ssh/README.md](https://gitlab.abinit.org/help/ssh/README.md)
