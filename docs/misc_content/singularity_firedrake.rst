.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:singularity:

Singularity
===========

`Singularity <https://sylabs.io/>`_ is as far as I can tell similar to `Docker <https://www.docker.com/>`_ and the idea is sometimes you want to have software working on a machine (e.g. a computer cluster) but it is difficult to install it natively on the cluster for whatever reason (e.g. no admin rights, clashing compilers/OS/libraries etc.) The idea here is you build the software package on a machine where you have the administrator rights as some image/container, then send and load that image/container to some production machine so that the image/container can piggyback on the available resources there.

Of course it doesn't have to be just used this way, e.g. you can make a singularity container (with a Linux OS) within a Mac and then use it on the same Mac, which bypasses some potentially annoying problems Macs have with compilers.

See the `Singularity <https://sylabs.io/>`_ website for the manual.

Key singularity commands used
-----------------------------

* ``[sudo] singularity build [--sandbox]``

build a sandbox, build an image from a sandbox, or build an image from a definition file

* ``[sudo] singularity shell [--writable] [--writable-tmpfs] [--overlay IMAGE]``

shell into the contained (image and/or sandbox), options respectively are

1. persistent write in sandbox/image itself
2. non-persistent write in sandbox/image itself
3. persistent write but in an overlay image (created separately)

* ``[sudo] singularity exec [--writable] [--writable-tmpfs] [--overlay IMAGE]``

run the executable build in container, with write privileges accordingly if required

Installing singularity
----------------------

Things needed really are a computer with admin/root privileges, and some familiarity with shell commands. On my laptop I run Linux Ubuntu 20.04 LTS, and doing ``uname -r`` tells me that my kernel version is ``5.4.0-40-generic``, so I'm going to do the following with Singularity ver. 3 (ver. 2 is old), specifically with version 3.5.2. 

.. warning ::

  The kernel version is important because it seems Singularity versions prior to and including 3.5.0 do not play with the 5.4 kernels because of some deflation problem (e.g. see `here <https://github.com/hpcng/singularity/issues/4801>`_). When trying to do something like ``singularity build`` (see later) the error you get is something like:
  
  .. code :: bash
  
    kernel reported a bad superblock for squashfs image partition
    possible causes are that your kernel doesn't support the compression algorithm or the image is corrupted
    
  Just use a newer version to fix this.

Installing singularity can be done through ``apt`` or ``yum`` or equivalent, though you might be limited by what versions are available. To build you own follow instructions on `here <https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps>`_ (you will need to also install the language ``Go``). In my case I effectively copy and pasted the code `here <https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps>`_ which gets me Singularity ver. 3.5.2.

Making the image
----------------

Once you have singularity the easiest (but by no means the most reproducible) way to do things is use the ``--sandbox`` option with the ``build`` command:

.. code :: bash
  
  sudo singularity build --sandbox sandbox_name docker://ubuntu:20.04
  
which downloads in this case ``ubuntu 20.04`` from ``docker`` to be used as the container operating system, which will be within a folder called ``sandbox_name`` in whatever directory you triggered the command from (e.g. ``./sandbox_name/usr/bin`` would be where the container binaries sit). Other options available may be found `here <https://sylabs.io/guides/3.5/user-guide/build_a_container.html>`_.

Then you can go into the container by

.. code :: bash
  
  sudo singularity shell --writable sandbox_name
  
The ``--writable`` command allows you to modify the things in the container as if you are whatever OS you decided to have (in this case Ubuntu 20.04 because I did ``docker://ubuntu:20.04`` in the ``build`` command). ***It is highly recommended you do this as ``sudo`` or ``root``*** (skip below warning if you don't care why).

.. warning ::

  First reason to do it as ``sudo`` is that to update the libraries and packages in the container presumably you would want to use package managers like ``apt`` (e.g. ``apt update && apt install python2.7``), which can't be done without ``sudo`` because the follow up ``dpkg`` command needs root privileges.
  
  Second is more subtle and should be fixable but I haven't bothered to figure out how to fix it. If you don't use ``sudo`` then when you shell in, you end up exactly where you triggered the command from (e.g. ``pwd`` will tell you something like ``/home/jclmak/singularity_images/sandbox_name``, compared to doing it with ``sudo``, where ``pwd`` should then be in ``/root``). If you do the former and start installing things then ``$PATH`` and other things you install will most likely have the wrong ``PATH`` when you build it into an image later to send somewhere else. For example:
  
  .. code :: bash
  
    ### without sudo ###
    (ubuntu)      singularity shell --writable sandbox_name/
    (singularity) pwd
                  >> /home/jclmak/singularity_images/sandbox_name
    
    # install Firedrake into what should be the container /home/firedrake
    (singularity) cd home/ && python3 firedrake-install
                  >> (whole load of stuff) DONE
    
    (singularity) ls /home/jclmak/singularity_images/sandbox_name/home  # because we can see the host folders
                  >> firedrake
    
    # build image
    (singularity) exit
    (ubuntu)      sudo singularity build testing.sif sandbox_name/
                  >> (builds sandbox into an image) DONE
    
    # send image
    (ubuntu)      rsync -arv testing.sif jclmak@hpc3.ust.hk:~/
                  >> DONE
    (ubuntu)      ssh jclmak@hpc3.ust.hk   # log onto external machine
    (hpc3)        module load singularity  # load the singularity module
    
    # look for where firedrake should have been installed in the container
    (hpc3)        singularity exec testing.sif ls /home
                  >> NOTHING 
    # because it actually got installed in /home/jclmak/singularity_images/sandbox_name/home which doesn't exist on the remote machine
    
    ### with sudo ###
    (ubuntu)      sudo singularity shell --writable sandbox_name/
    (singularity) pwd
                  >> /root
    
    # install Firedrake into what should be the container /home/firedrake
    (singularity) cd ../home/ && python3 firedrake-install
                  >> (whole load of stuff) DONE
    
    (singularity) ls /home/jclmak/singularity_images/sandbox_name/home
                  >> NOTHING (because folder does not exist in the container)
    (singularity) ls /home/
                  >> firedrake
    
    # build image
    (singularity) exit
    (ubuntu)      sudo singularity build testing.sif sandbox_name/
                  >> (builds sandbox into an image) DONE
    
    # send image
    (ubuntu)      rsync -arv testing.sif jclmak@hpc3.ust.hk:~/
                  >> DONE
    (ubuntu)      ssh jclmak@hpc3.ust.hk   # log onto external machine
    (hpc3)        module load singularity  # load the singularity module
    
    # look for where firedrake should have been installed in the container
    (hpc3)        singularity exec testing.sif ls /home/
                  >> firedrake
  
  Doing it without sudo seems to let it work within the machine where the container was created (because the folders are there as part of the sandbox on the creation machine), but doesn't work when you send the container elsewhere because it can no longer find the files in the sandbox.
  
Once you are in the shell (as ``sudo`` probably) you can do whatever you might normally do. For example doing ``apt update && apt install python3.8`` installs ``python3`` (a link of ``python3.8``) into ``/usr/bin``.

Using the image
---------------

Then you can trigger what you installed by (as a sandbox here)

.. code :: bash

  (ubuntu)      sudo singularity shell --writable sandbox_name/
  (singularity) apt update && apt install python3.8
                >> (install some packages) DONE
  (singularity) exit
  (ubuntu)      singularity exec sandbox_name/ /usr/bin/python3 --version
                >> 3.8.?
                
``/usr/bin`` should be in ``$PATH`` anyway so you could do

.. code :: bash

  (ubuntu)      singularity exec sandbox_name/ which python3
                >> /usr/bin/python3
  (ubuntu)      singularity exec sandbox_name/ python3 --version
                >> 3.8.?
                
To use it elsewhere, ``build`` it as an image via

.. code :: bash

  (ubuntu)      sudo singularity build sandbox_image.sif sandbox_name/
                 >> DONE
  (ubuntu)      singularity exec sandbox_image.sif python3 --version
                 >> 3.8.?

The image can be sent to other machines, e.g.

.. code :: bash

  (ubuntu)      rsync -arv sandbox_image.sif jclmak@hpc3.ust.hk:~/
  (ubuntu)      ssh jclmak@hpc3.ust.hk
  (hpc3)        module load singularity
  (hpc3)        singularity exec sandbox_image.sif python3 --version
                >> python 3.8.?
  (hpc3)        python3 --version # calling the HPC3 python3
                >> python 3.6.8
                
If doing something like the above where you are just running an executable, the outputs will be in the host machine and not in the container (because the container should now be an immutable image, and the image has not been specified as ``--writable``). For the Firedrake example below the program wants to write some cache but directly into the image itself, so I have a slight work around there.

Working example 1: Julia
------------------------

`Julia <https://julialang.org/>`_ is a new computing language that appears to combine the user benefits of Python and computationally optimised languages such as Fortran.

There is already an existing guide for `installing Julia as a container <https://github.com/sylabs/examples/tree/master/lang/julia>`_. The thing to be aware of is that it may be worthwhile creating an overlay with more space (see next working example) or keeping the sandbox folder lying around to make Julia images with on the creation computer in case more packages need to be downloaded, since the singularity image shouldn't be written to.


Working example 2: Firedrake
----------------------------

`Firedrake <https://www.firedrakeproject.org/>`_ is a finite element based computational framework with automatic code generation capabilities. High/abstract level code in Python is converted into machine code at the production stage for performance reasons. Parallelised code is automatically generated when more cores are provided at the run/compile stage, without a need to change any of the high/abstract level code.

The following will get Firedrake working in the container (given in ``--sandbox`` mode here for the time being):

.. code :: bash

  (ubuntu)      sudo singularity build --sandbox firedrake docker://ubuntu:20.04
  (ubuntu)      sudo singularity shell --writable firedrake
  (singularity) apt update && apt install curl python2.7 python3.8
  
  # I am going to install it in /home rather than /root because 
  # /root can't be seen when the container will be run in due course
  (singularity) cd /home   
  (singularity) curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
  (singularity) python3 firedrake-install --disable-ssh
                >> YES TO INSTALL ALL THE SUGGESTED DEPENDENCIES
                >> TAKES BLOODY AGES (>= 1hr) TO INSTALL (PETSc)
  (singularity) ls -lh /home
                >> drwxr-xr-x root root firedrake/
                
So here the owner and the group of ``firedrake/`` is both ``root`` (because I installed it as ``sudo``) and we see that we have 755 access (no write for non-owner).

This will be fine ***if*** not for the fact that firedrake needs to write to ``./firedrake/.cache``, so no one but owner of the folder (``root``) can actually run with the firedrake modules once it's wrapped into a container because only ``root`` can create the ``.cache`` folder. The way I got round this to do the lazy thing:

.. code :: bash

  (singularity) chmod -R 777 /home/firedrake && ls -lh /home
                >> drwxrwxrwx root root firedrake/
                
So while the firedrake folder can now in principle be written to, in reality when it is packaged in singularity image the stuff within the image itself cannot be touched without ``sudo`` access (because it was created with ``sudo`` and is supposed to be immutable) and the ``--writable`` flag.

This creates another dilemma because then we have nowhere to write ``.cache`` since the image is supposed to be untouchable. There are two ways around this, and I recommend the second one for firedrake:

1. when ``singularity exec`` is called, give it the ``--writable-tmpfs`` flag, which will provide a non-persistent write option (so it writes to some temporary space but everything is discarded when the ``singularity exec`` command is finished)
2. create an overlay so the writing is done to some allocated space, but this overlay can be re-used, which means the cache that is created is persistent across runs. To do this,

.. code :: bash

  (singularity) ls -lh /home
                >> drwxrwxrwx root root firedrake/
  (singularity) exit
  (ubuntu)      sudo singularity build firedrake.sif firedrake/
                >> DONE
                
  # create an 500MB overlay with dd and mkfs-ext3 (no need for sudo access)
  (ubuntu)      dd if=/dev/zero of=firedrake_cache bs=1M count=500 && mkfs.ext3 firedrake_cache
  (ubuntu)      singularity exec --overlay firedrake_cache firedrake.sif /home/firedrake/bin/python linear_SWE2d.py
  
where ``linear_SWE2d.py`` is a file I have that solves something so requiring a write to ``firedrake/.cache``. If you shell into the image with the overlay option you should be able to see ``firedrake/.cache`` there (and it should not exist in the over shelled into without the overlay).

Installation can be done through a definition file too (see below), triggered by ``sudo singularity build firedrake.sif firedrake.def``. It takes upwards of an hour on my laptop largely because PETSc takes ages to build. The resulting ``firedrake.sif`` file is about 2GBs (can probably be smaller, I haven't figured out which of the files in src/PETSc are safe to delete...)

.. code :: bash

  BootStrap: docker
  From: ubuntu:20.04

  %post
      export DEBIAN_FRONTEND=noninteractive
      apt -y update
      apt -y upgrade
      apt -y install git curl nano python2.7 python3.8 build-essential
      apt -y install dialog apt-utils autoconf automake bison flex cmake gfortran 
      apt -y install libblas-dev liblapack-dev libtool
      apt -y install python3-dev python3-pip python3-tk python3-venv
      apt -y install zlib1g-dev libboost-dev
      cd /home/
      curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
      python3 firedrake-install --disable-ssh
      
      echo "changing permissions to enable temporary writes if need be (particularly in .cache)"
      chmod -R 777 /home/firedrake

  %environment
      export LC_ALL=C
      export PATH="/home/firedrake/bin:$PATH"
      
  %runscript
      echo "in firedrake.sif"
      echo "attempting to call python --version"
      python --version

  %labels
      Author Julian Mak

