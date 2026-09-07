
import sys
import numpy as np
from mbc_py_interface import mbcNodal

import precice
from precice import *


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
      # rot = 0x100 -> orientation vector
      # rot = 0x200 -> orientation matrix
      # rot = 0x400 -> Euler 123

  #> Methods ----------------------------------------------------
  #> initialize -------------------------------------------------
  def initialize(self, path='.mbdyn/mbdyn.sock', \
                       verbose=1, \
                       nnodes=1, \
                       accels=0, \
                       forceType='Nodal', \
                       dumpAuxFile=False, \
                       dumpAuxFilen='./../nnodes.dat'):

    #> Initialize communication
    self.socket = self.Socket( path=path, verbose=verbose, \
                               nnodes=nnodes, accels=accels )
    #print(' self.socket.accels: ', self.socket.accels)
    if ( forceType == 'Nodal' ):
      self.nodal = mbcNodal( self.socket.path, \
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
      # print(' self.data["Position"]["data"]: ', self.data["Position"]["data"] )
      # print(' self.data["Velocity"]["data"]: ', self.data["Velocity"]["data"] )
      # print(' self.data["Rotation"]["data"]: ', self.data["Rotation"]["data"] )
      # print(' self.data["AngularVelocity"]["data"]: ', \
      #         self.data["AngularVelocity"]["data"] )

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

  #> sendMBDynMesh ----------------------------------------------
  # *** old *** old function, aiming at using the initial configuration of a
  # MBDyn model as the reference configuration for the "Lagrangian" coupling
  # trhough PreCICE
  def sendMBDynMesh(self, chord, twist):
 
    n = self.socket.nnodes
    #> ~~~ PreCICE-aux communication: start ~~~~~~~~~~~~~~~~~~~~~
    #> Precice params: communication
    comm_rank = 0; comm_size = 1
    config_file_name = "./../precice-config-aux.xml"
    solver_name   = "MBDyn"
    mesh_name_mbd = "AuxilMBDynMesh"
    # mesh_name_flu = "AuxilFluidMesh"
    pos_name      = "Position"; rot_name    = "Rotation"
    cho_name      = "Chord"   ; twi_name    = "Twist"
    dum_name      = "Dummy"
    
    print(" Configure preCICE-aux ...", end='')
    interface = precice.Interface( solver_name, config_file_name, \
                                   comm_rank, comm_size )
    dimensions  = interface.get_dimensions(); nd = dimensions
    mesh_id_mbd = interface.get_mesh_id( mesh_name_mbd )
    pos_id    = interface.get_data_id( pos_name, mesh_id_mbd )
    rot_id    = interface.get_data_id( rot_name, mesh_id_mbd )
    cho_id    = interface.get_data_id( cho_name, mesh_id_mbd )
    twi_id    = interface.get_data_id( twi_name, mesh_id_mbd )
    dum_id    = interface.get_data_id( dum_name, mesh_id_mbd )
    print(" done. \n")
   
    #> Define reference configuration in the PreCICE global reference
    # frame
    rr_mbd = np.zeros((n, nd))
    for i in np.arange(n):    # todo with reshape
      rr_mbd[i,:] = np.array([ i, 0., 0. ]) # nodal.n_x[i*nd:(i+1)*nd]
    node_id_mbd = interface.set_mesh_vertices( mesh_id_mbd, rr_mbd )
    
    #> Fields
    pos = np.ones((n, nd)); cho = np.ones((n))
    rot = np.ones((n, nd)); twi = np.ones((n))
    
    for i in np.arange(n):
      pos[i,:] = self.nodal.n_x[i*nd:(i+1)*nd]
      rot[i,:] = self.nodal.n_theta[i*nd:(i+1)*nd]  # check !
      cho[i  ] = chord
      twi[i  ] = twist
    
    # print(' rr_mbd:', rr_mbd)
    # print(' pos   :', pos)
    # print(' rot   :', rot)
    
    dum = np.ones((n)) # dummy field read from fluidSolver
    
    dt_precice = interface.initialize()
    is_ongoing = interface.is_coupling_ongoing()
    print(" is_ongoing: ", is_ongoing)
    
    # cowic = interface.
    
    while ( is_ongoing ):
    
      interface.write_block_vector_data( pos_id, node_id_mbd, pos ); # print(' pos: ', pos)
      interface.write_block_vector_data( rot_id, node_id_mbd, rot ); # print(' rot: ', rot)
      interface.write_block_scalar_data( cho_id, node_id_mbd, cho ); # print(' cho: ', cho)
      interface.write_block_scalar_data( twi_id, node_id_mbd, twi ); # print(' twi: ', twi)
    
      dum = interface.read_block_scalar_data( dum_id, node_id_mbd ); # print(' dum: ', dum)
    
      dt_precice = interface.advance(dt_precice)
      is_ongoing = interface.is_coupling_ongoing()


    interface.finalize()
    print(" Finalize auxiliary coupling ")
  
    return pos, rot, interface
      
