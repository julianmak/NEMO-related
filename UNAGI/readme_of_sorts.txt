Rough step-by-step guide in building UNAGI (for NEMO ver >3.6)

* overall plan of action
  -- make use of the NEMO "TOOLS/DOMAINcfg" to generate a "domcfg.nc" that
     specifies all the model parameters etc.
     -- 'cause ain't nobody got time to write the model in from scratch (in
        "USR/") or hack up a "domcfg.nc" from scratch"...
  -- read that "domcfg.nc" to create forcing fields etc.
  -- add in some modifications in "MY_SRC" as required
  
* NOTE: I am not completely precise about specifying paths where the commands
        should be trigger below, so if you blindly follow this it will most
        likely not work (you shouldn't ever be blindly doing things generally
        anyway...)
  
1) create some directory and dump "gen_NEMO_UNAGI_fields.ipynb" there, go into
   the file (or just copy some code from there) to specify the domain grid etc.
   and bathymetry, and run the bit of code that gives you a "bathy_meter.nc"
   
                       OUTCOME: a "some_dir/bathy_meter.nc"
   
2) go into the relevant "TOOLS/DOMAINcfg" (or "tools/DOMAINcfg" in the newer
   NEMOs) and compile it as e.g.
   
                    ./maketools -n DOMAINcfg -m HKUST_HPC3_xios1
   
   -- according to the README you need to use XIOS1.0, so you may need to
      rebuild this yourself
     
                       OUTCOME: a "DOMAINcfg/make_domain_cfg.exe"
                       
3) copy the "some_dir/bathy_meter.nc" you made in 1) into the "DOMAINcfg"
   folder, and in "DOMAINcfg" probably do

                    cp -p namelist_cfg namelist_cfg.orig
                    
   just to be safe, and modify the "namelist_cfg" so that the parameters in
   there are consistent with the numbers given in
   "some_dir/gen_NEMO_UNAGI_fields.ipynb"
   
   -- there is a sample of "namelist_cfg" in this repository for reference, do a
      visual compare; the one that isn't clear is "jphgr_msh = 3" (beta-plane
      with regular spacing according to "src/dom_oce.f90")
      
                       OUTCOME: a modified "DOMAINcfg/namelist_cfg" file
                       
4) do (either on your laptop or on some other machine)

                    mpirun -np 1 ./make_domain_cfg.exe
                    
   to get the "domcfg.nc"
   
   -- if it fails, look into "ocean.output" and trace errors as usual
      
                       OUTCOME: a "DOMAINcfg/domcfg.nc" file
                       
5) now copy the "DOMAINcfg/domcfg.nc" back into the directory in "some_dir", and
   run the rest of the script accordingly to generate the initial states and
   forcing files
   
   -- for this experiment we only have zonal wind
      
                       OUTCOME: a "some_dir/domcfg.nc" file
                                a "some_dir/forcing.nc" file
                                a "some_dir/state.nc" file
                                
   -- I prefixed these files with experiment name and suffixed with resolution,
      something like e.g. "domcfg_UNAGI_R100.nc", for book keeping purposes and
      partly because later on when you run the experiment the "namelist_cfg"
      will ask you for the file name and you specify it as e.g.
      
                    cn_domcfg = "domcfg_UNAGI"
      
6) we now have the files to run the experiment but we need to compile the
   experiment first

   -- in the NEMO ver x (3.6 < x < 4.0; if folders have capitalised names this
      is probably you) you can just do something like
      
                    ./makenemo -r GYRE_PISCES -n UNAGI -m HKUST_HPC3 -j0
                    
      which will create a folder but not compile  
      
   -- if you are in NEMO 4.0 onwards then I tend to do the following:
   
                * go into ref_cfgs.txt and add in an entry for UNAGI
                     (just OCES is fine, UNAGI not currently configured for TOP)
                * go out and do
                     rsync -arv GYRE_PISCES UNAGI
                * go into UNAGI, change the cpp file name so it is suffixed with
                  UNAGI
                * for now within "cpp_UNAGI.fcm", get rid of "key_top" and
                  probably put in "key_nosignedzero" (compiler dependent)
                  
      then do
      
                    ./makenemo -r UNAGI -m HKUST_HPC3 -j0
                    
      which will create some folders within UNAGI but not compile
      
   -- copy some content from "respository/MY_SRC/*" to "UNAGI/MY_SRC/", and
      compile/debug as usual (e.g. ./makenemo -r UNAGI -m HKUST_HPC3 -j4)
                    
                       OUTCOME: a "UNAGI/EXP00/{opa,nemo}" executable
                       
7) I normally copy the contents in "EXP00" into some other folder (e.g.
   "EXP_R100"), so assuming you do that, then copy the created
   "{state,forcing,domcfg}.nc" files into (say) "EXP_R100" so NEMO can read it,
   and modify "namelist_cfg" as usual for setting up experiment
   
   -- there is a "EXP00/namelist_cfg" in repository, check differences for
      seeing what I changed (the one I have is for NEMO 3.7 I think so there are
      subtle differences if you are using NEMO 4.0 because some namelist entries
      have changed)
      -- if the run fails to intialise, look in "ocean.output" to have an idea
         where it failed
      -- you need to tell NEMO to read stuff by
      
                 &namcfg  
                    ln_read_cfg = .true.
                       cn_domcfg = "domcfg_UNAGI"          (that's what I called it) 
                       
                 &namtsd
                    sn_temp = "state_UNAGI" etc. "toce"    (my variable name was toce)
                    sn_sal  = "state_UNAGI" etc.
                    
                 &namsbc
                    ln_usr  = .false.
                    ln_flx  = .true.  (because I am specifying wind stress here)
                    ln_blk  = .false.
                    ln_ssr  = .true.  (for surface restoring)
                 
                 &namesbc_flx
                    sn_utau = "forcing_UNAGI" etc. "utau"  (my variable name was utau)
                    sn_vtau = etc.
                    
                 &namsbc_ssr
                    sn_sst  = "state_UNAGI" etc. "sst"     (I have sst there too)
                    sn_sst
                      
      -- I hacked in an enhanced vertical diffusivity over a sponge region for
         this experiment (cf. Munday et al 2015, GRL), with magnitude specified
         in
         
                 &namzdf
                    rn_avt_amp = 250.0   (i.e. 250x amplification over a sponge)
                    
         if you don't like that etc you can always get rid of the stuff in
         "MY_SRC" accordingly
         
                       OUTCOME: an experiment that we solve all the world's problems
