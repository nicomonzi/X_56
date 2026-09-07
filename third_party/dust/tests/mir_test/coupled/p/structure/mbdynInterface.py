
import sys
from mbc_py_interface import mbcNodal

import precice
from precice import *

import numpy as np

class MBDynInterface:
  """ MBDyn interface collecting kinematic and forces """
  #> Constructor ------------------------------------------------
  def __init__(self):
    self.initialized = False
    self.data = { "Position":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "Velocity":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "Rotation":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "AngularVelocity":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "Force":{ \
                    "type":"vector", "io":"read" , "data":[] }, \
                  "Moment":{ \
                    "type":"vector", "io":"read" , "data":[] } \
                }

  #> Inner class ------------------------------------------------
  class Socket:
    """ Class containing the socket parameters """
    """ for comm between MBDyn and mbc_py      """
    def __init__(self, \
                 path="", host="", port=0, \
                 timeout=-1, verbose=0, data_and_next=1, \
                 refnode=0, nnodes=1, labels = 0, rot=0x100, \
                 accels=0 ):
      self.path          = path
      self.host          = host          
      self.port          = port          
      self.timeout       = timeout       
      self.verbose       = verbose       
      self.data_and_next = data_and_next 
      self.refnode       = refnode       
      self.nnodes        = nnodes        
      self.labels        = labels        
      self.rot           = rot           
      self.accels        = accels

      print(' in Socket __init__. self.nnodes: ', self.nnodes)

  #> Methods ----------------------------------------------------
  #> initialize -------------------------------------------------
  def initialize(self,path='.mbdyn/mbdyn.sock', \
                      verbose=0, \
                      nnodes=1, \
                      accels=0, \
                      forceType='Nodal', \
                      dumpAuxFile=False, \
                      dumpAuxFilen='./../nnodes.dat'):

    #> Initialize communication
    self.socket = self.Socket( path=path, verbose=0, \
                               nnodes=nnodes, accels=accels )

    if ( forceType == 'Nodal' ):
      self.nodal = mbcNodal(self.socket.path, \
                            self.socket.host, \
                            self.socket.port, \
                            self.socket.timeout, \
                            self.socket.verbose, \
                            self.socket.data_and_next, \
                            self.socket.refnode, \
                            self.socket.nnodes, \
                            self.socket.labels, \
                            self.socket.rot, \
                            self.socket.accels )
      print(' self.socket.nnodes: ', self.socket.nnodes)
      #> Negotiate communication and receive kinematics from MBDyn
      self.nodal.negotiate()
      self.nodal.recv()

      if ( dumpAuxFile ):
        fid = open(dumpAuxFilen, "w")
        fid.write('%d' % ( self.socket.nnodes ) )
        fid.close()

      #> Save coordinates of the MBDyn nodes
      n = self.socket.nnodes; nd = 3
      self.data["Position"       ]["data"] = np.reshape( self.nodal.n_x    , (n, nd) )
      self.data["Velocity"       ]["data"] = np.reshape( self.nodal.n_xp   , (n, nd) )
      self.data["Rotation"       ]["data"] = np.reshape( self.nodal.n_theta, (n, nd) )
      self.data["AngularVelocity"]["data"] = np.reshape( self.nodal.n_omega, (n, nd) )

      #> Initialize read fields
      self.data["Force" ]["data"] = np.zeros((n, nd))
      self.data["Moment"]["data"] = np.zeros((n, nd))

      #> Initialize nodal.n_f, n_m
      self.nodal.n_f[:] = 0.;  self.nodal.n_m[:] = 0.
      # self.nodal.send(True)

    else:
      sys.exit(" So far, only forceType='Nodal' is implemented")

  #> finalize ---------------------------------------------------
  def finalize(self):
    self.nodal.destroy()
    
  #> refConfigNodes ---------------------------------------------
  def refConfigNodes(self, filen='./refConfigNodes.in'):
    rr = np.loadtxt( fname = filen )

    # loadtxt import a (1,n) array as a vector (with len(rr.shape)=1),
    # but an array is needed. If only one node is defined, save rr in
    # a (1,3) array
    if ( len(rr.shape) == 1 ):
      rr = rr.reshape( 1,3 )

    return rr
