if [ ! "$MODULES_LOADED" ] ; then
   module load gnu
   module load mpich
   module load cuda
   export MODULES_LOADED=yes
fi
