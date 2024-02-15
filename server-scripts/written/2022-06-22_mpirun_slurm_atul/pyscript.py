from mpi4py import MPI

import socket

import subprocess

comm = MPI.COMM_WORLD

# The variable comm.rank contains a different number for each node. This will be output in out.txt file:
print( "Hello! I'm rank " + str(comm.rank) +' at hostname ' + str(socket.gethostname()) )

# This is a dummy worker. I'm simply calling "python worker.py" which runs 1000's of for loops to keep the nodes busy.
# Replace the stuff inside subprocess.call() with the real work that needs to be done.
# I'm passing a command line argument to worker.py here which is simply comm.rank. In a real implementation, this
# behavior can be used to do more interesting things like choosing with folder to process.
# We can also have an if statement that changes based on the value of comm.rank and does different things.
cmd = subprocess.call(['python', 'worker.py', str(comm.rank)])

# This is not strictly necessary, but it just shows how to create different
# output files if needed, one for each node:
with open("log.out"+str(comm.rank), 'w') as f:
    f.write(str(comm.rank) + str(socket.gethostname()))
