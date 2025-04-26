.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:jupyter:

Jupyter notes
=============

Jupyter notebook (and lab)
--------------------------

Remote instance and login
-------------------------

Sometimes it is useful to open a Jupyter instance on the computer cluster but open it on the local machine, which saves having to download the data onto the machine and eating up space (and there is possibility to do parallelised data processing through `DASK <https://www.dask.org/>`__ for example). The way below is some overall instructions; details will differ depending on firewall configurations and machine details.

0. With anaconda, miniconda or similar, make a virtual environment and install jupyter and other packages accordingly.
1. SSH or similar onto the remote machine, and submit a job that opens that virtual environment, trigger jupyter with some specified port (normally ``8888``, but you can choose), and keep it open. See below for a sample submit script (in my case the remote machine is called ``hpc3`` and the environment is called ``py311``).
2. The job should be running on some node(s). Still on the remote machine, query the IP address of the master node via ``nslookup hhnode-ib-21`` or ``dig +short hhnode-ib-21.local`` (replace ``hhnode-ib-21`` with whatever the node name is). Also query the security token that gets generated for use later; it should be in the output (in my case it's in ``stdouterr_${job_number}``).
3. Go back to the local machine, and open a SSH tunnel to that IP address with 

.. code-block:: bash

  ssh -N -f -L ${LOCAL_PORT}:${IP_ADDRESS}:${REMOTE_PORT} ${USER}@${CLUSTER}

So as an example, ``ssh -N -f -L 4167:10.1.2.126:4167 jclmak@hpc3.ust.hk`` means I bind the remote port ``4167`` at the IP address ``10.1.2.126`` to my local machine's port ``4167`` using my appropriate credentials on the cluster.

.. note::

  I use port ``4167`` because if I have a local instance of jupyter that usually by default opens at ``8888`` and lead to a clash. I could have done ``ssh -N -f -L 4167:10.1.2.126:8888`` or similar to avoid the clash I suppose.

4. Use ``lsof -i :{$LOCAL_PORT}`` to see if the port is open. If there are already ones open, you might want to kill those to avoid clashing.
5. Open a browser and enter ``http://localhost:4167`` (or whatever you decided to substitute ``${LOCAL_PORT}`` for. If all goes well then a jupyter instance will open but ask you for a token; enter the token from above in and you should have a jupyter instance working on the external machine but controlled on your the local machine.

.. note::

  The above works ok for me whether my laptop firewall is on or not. However, I seem to need to enter the token manually, while for my post-doc she could do something like ``http://hhnode-ib-21:8888/?token=WHATEVER`` and get on directly. Not sure what the deal is.
  
Sample submit script for a SLURM system:

.. code-block:: bash

  #SBATCH -o stdouterr_%j # output and error file name
  #SBATCH -n 1            # total number of mpi tasks requested
  #SBATCH -N 1            # total number of nodes requested
  #SBATCH -p cpu          # queue (partition) -- standard, development, etc.
  #SBATCH -t 12:00:00     # maximum runtime
  ##SBATCH -x hhnode-ib-[10,24,52,53,48,36]    # avoid some nodes

  # Setup runtime environment if necessary
  # module load anaconda3  # commented here because I use a separate miniconda manager

  # just make sure the environment really is off otherwise it seems to fail to load the environment here
  python --version
  source /home/jclmak/miniconda3/bin/activate
  source /home/jclmak/miniconda3/bin/activate py311
  python --version

  # Set Tunneling information
  node=$(hostname)
  user=$(whoami)
  cluster="hpc3.ust.hk"
  port=4167  # provide a specific port to avoid possible clashing

  # Print tunneling instructions
  echo -e "

  # Command to create SSH tunnel:
  ssh -N -f -L ${port}:IP_TO_FILL:${port} ${user}@${cluster}

  # Use a browser on your local machine to go to:
  http://localhost:${port}/

  use output of dig +short ${node}.local to replace IP_TO_FILL

  "

  jupyter-notebook --no-browser --ip=${node} --port=${port}

  # keep job alive so it can be tunneled in
  sleep 36000
