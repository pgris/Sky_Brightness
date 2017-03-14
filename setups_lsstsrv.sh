#!/bin/bash

source ../../lsst_stack/loadLSST.bash

#does this needs to be done each time or once for all ?
#eups declare mysql system -m none -r none -c

export PYTHONPATH=../palpy-1.7.0_installed/lib/python2.7/site-packages:$PYTHONPATH
#setup sims_operations phg_sims_operations_v11
export PYTHONPATH=../pyephem-3.7.6.0_installed/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=/data/pgris/sims_operation/healpy-1.9.0_installed/lib/python2.7/site-packages:$PYTHONPATH
#export PYTHONPATH=../pymssql-2.0.1_installed/lib/python2.7/site-packages/pymssql-2.1.2.dev0-py2.7-linux-x86_64.egg
export PYTHONPATH=/data/pgris/sims_operation/pymssql_installed/lib/python2.7/site-packages:$PYTHONPATH
setup  -k -r ../freetds
#setup pymssql phg_pymssql_v11
setup -k -r ../pykg_config
setup -k -r ../healpy
setup -k -r ../sims_utils
setup -k -r ../sims_catalogs_generation
setup -k -r ../sims_sed_library
setup -k -r ../sims_dustmaps
setup -k -r ../throughputs
setup -k -r ../sims_maps
setup -k -r ../sims_photUtils
setup -k -r ../sims_coordUtils
setup -k -r ../sims_maf
setup -k -r ../sims_operations
export PYTHONPATH=../sims_operations/python/lsst/sims/operations:$PYTHONPATH
#eups declare pyephem system -m none -r none -c