.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

git (or bitbucket) commands
===========================

[TO TIDY] Git is phasing out password logins on terminal access, so you either
have to do two factor authorisation (2FA), use a SSH key, or a token (there are
others presumably). The following bits of scrap code documents how to use a
token.

1. log into Git in the browswer, click your profile picture, then Settings -> Developer settings -> Personal access tokens
2. give your token a descriptive name and give permissions for the access token
3. when you click OK you should get to a screen with that tells you to copy the token (don't close this page yet!)
4. open terminal and do

.. code-block:: bash

  git config --global credential.helper cache
  git config --global crediential.helper 'cache --timeout=31104000'
  
where you can change the ``timeout`` entry to something that works for you (something large if you want to keep the token active, units are in seconds)
5. go back to webpage, copy the access token, git as normal, but when it asks you for a password, paste the access token in instead
6. if it worked properly then now you get to bypass the username and password typing until the timeout period

At any point you can revoke the access token on the Git webpage.


