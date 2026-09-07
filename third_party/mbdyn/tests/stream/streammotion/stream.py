from __future__ import print_function
from subprocess import call
from time import sleep, clock
from tempfile import TemporaryFile
from socket import socket, AF_INET, SOCK_STREAM
from struct import unpack

f1 = TemporaryFile()
command = 'mbdyn -s -f stream &'
elapsed = 0.
inc = 0.001  # same as mbdyn timestep
call(command, shell=True, stdout=f1)
start = clock()
sock = socket(AF_INET, SOCK_STREAM)
sock.bind(('', 5501))				
sock.listen(5)
newSocket, address = sock.accept()
while True:
	while clock() < start + elapsed:
		pass
	elapsed += inc
	receivedData = newSocket.recv(24)
	if not receivedData:
		break
	else:
#		print("{}".format(len(receivedData)))
		print("Data={} elapsed={}".format(unpack('d'*3, receivedData), elapsed))
newSocket.close()
sock.close()
if f1:
	f1.seek(0)
	print("{}".format(f1.read()))
f1.close()

