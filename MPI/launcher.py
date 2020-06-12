import os

#----------------------------------------------------------------
# Deployment function: launch a MPI4PY pgm on a set of cluster nodes
# + autocompute the number of allocated nodes
#----------------------------------------------------------------
def deployOS():
  # - "socket" deployment strategy on Kyle
  # nbProcesses = 2*nbNodes  # Max
  # nbProcesses = 4
  # os.system("mpirun -np " + str(nbProcesses) +
  #           " -map-by ppr:4:socket -rank-by socket -bind-to socket" +
  #           " python3 parallel.py")
  # - "core" deployment strategy on Kyle
  # nbProcesses = 16*nbNodes  # Max
  nbProcesses = 4
  os.system("mpirun -np " + str(nbProcesses) +
            " -map-by ppr:1:core -rank-by core -bind-to core" +
            " python3 parallel.py")


#-----------------------------------------------------------------
# Main code
#-----------------------------------------------------------------

print("Deployment using OS module")
deployOS()




