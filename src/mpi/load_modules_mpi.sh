if [ ! "$MODULES_LOADED" ] ; then
   module load gnu
   module load mpi
   module load cuda
   export MODULES_LOADED=yes
fi
