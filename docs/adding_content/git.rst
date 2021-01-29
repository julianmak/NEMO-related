.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Git commands
============

`Git <git-scm.com>`_ is s version control software. Similar software exist (e.g.
`Mercurial <www.mercurial-scm.org >`_, `Subversion <subversion.apache.org>`_),
but I almost exclusively use Git now for my own things (NEMO uses subversion but
the only thing I ever do is ``svn checkout LOCATION -r VERSION``, so it's not
really using it...) 

I personally use Git for backup mostly, occasionally reverting files, as well as
hosting websites (e.g. `here <https://julianmak.github.io/>`_ and `here
<https://jmak-omfg.github.io/>`_). For my kind of files (mostly text files, as
source LaTeX files, html or bits of code) I find it much more convenient and
safer than saying using Dropbox (manual version control is too error prone for
me). I haven't personally used Git that much in terms of collabroative work at
the moment, so the commands below are going to be skimpy on those related
commands.

Repositories
------------

`Github <https://github.com/>`_ is my go to for making repositories, partly
because it can render Jupyter notebooks I use a lot. Github used to only have
public repositories, but now they have private ones too so I migrated from
`Bitbucket <https://bitbucket.org/product/>`_. Make an account and create a
repository so there is a target to push and pull files from.

Keep the files small! Github doesn't accept anything larger than 100 Mbs I
think. e.g. commit LaTeX source files but not necessarily the compiled version.

Basic commands
--------------

The regular commands I use are:

* ``git add`` registers files that Git should track and note changes
* ``git commit -m "SOME DEEP MESSAGE"`` actually registers the changes made since last commit
* ``git push [origin master]`` pushes the commits up to the repository
* ``git pull`` pulles the commits from the repository to the local computer

Occasionally I screw something up so I need to do:

* ``git mv`` to move the files around by telling Git to still track them
* ``git rm [--cached] FILES`` to make git to stop tracking the files (the ``--cached`` is so that the physical files are not removed; leave it out if you actually want to get rid of it)
* ``git checkout HEAD FILES`` if I screw up the ``rm``, ``mv`` or ``git rm`` commands to recover the removed physical files. ``HEAD`` can be replaced by revision number
* ``git log`` to check the log of commits and revision numbers

[TO ADD] some branching and merging commands

Access tokens
-------------

Git is phasing out password logins on terminal access, so you either
have to do two factor authorisation (2FA), use a SSH key, or a token (there are
others presumably). The following bits of scrap code documents how to use a
token.

1. log into Git in the browswer, click your profile picture, then Settings -> Developer settings -> Personal access tokens
2. give your token a descriptive name and give permissions for the access token
3. when you click OK you should get to a screen with that tells you to copy the token (don't close this page yet!)
4. open terminal and do

.. code-block:: bash

  git config --global crediential.helper 'cache --timeout=31104000'
  
where you can change the ``timeout`` entry to something that works for you (something large if you want to keep the token active for longer, units are in seconds). I tried using ``store`` on Ubuntu but it doesn't seem to do anything (``store`` saves an extra file with the credentials in) 

5. go back to webpage, copy the access token, git as normal, but when it asks you for a password, paste the access token in instead
6. if it worked properly then now you get to bypass the username and password typing until the timeout period

At any point you can revoke the access token on the Git webpage.


