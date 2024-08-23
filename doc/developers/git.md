ABINIT developer path to Git

# Introduction

For a long time, ABINIT has relied on the [Bazaar](https://en.wikipedia.org/wiki/GNU_Bazaar)  Version Control System 
(VCS). However, this software project is apparently not maintained anymore.
By contrast, [git](https://git-scm.com/) shows advantages with respect to other VCS, in 
particular, an impressive speed, thanks to its commit storage model. See 
[Wikipedia](https://en.wikipedia.org/wiki/Git_%28software%29): "The Eclipse 
Foundation reported <...> that as of May 2014, Git is now the most widely used source code management tool, with 42.9% of professional software developers reporting that they use Git as their primary source control system”.

So, after consulting the most active developers of the ABINIT community, the decision has been taken in 2015 to switch to git.

# What should I do to jump in ?

It will take a few hours to learn git and the corresponding ABINIT development model. 
- First step : become familiarized with the basics of git. This means simply: read chapters 1 to 3 of the "Reference Manual" of git, that is available at [https://git-scm.com/doc](https://git-scm.com/doc). Install a version of git on your local development computer, and play a bit with it. This will likely take two hours.
- Second step : become familiarized with the so-called [gitflow](http://nvie.com/posts/a-successful-git-branching-model). Indeed, the working mode that 
  is adopted for ABINIT is quite close to this flow (see later).
- Third step : spend some time to get familiarized with gitlab. in particular,
  please, log in on page [https://gitlab.abinit.org](https://gitlab.abinit).
  org). In order to do this, you have first to arrange an access with 
  Jean-Michel Beuken (jean-michel.beuken@uclouvain.be): send him an e-mail 
  and he will create for you an SSH access.
- Fourth step, read [the specificities of the use of git for the ABINIT 
  project](developers_git_specificities_git_abinit_standard_names_projects_branches.md).
- Fifth step : you might have a look at the [blueprint from Yann, entitled 
  "Gitlab specifications"](maintainers_blueprints_gitlab_specs.md).
- Sixth step : practice !


