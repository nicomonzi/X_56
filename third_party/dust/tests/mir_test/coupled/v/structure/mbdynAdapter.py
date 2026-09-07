
import numpy as np
import os
import precice
import pdb
import time
class MBDynAdapter:
  class Participant:
    def __init__(self, name='MBDyn'):
      """ Do nothing """
      self.name  = name
      self.field = [] 

    class Mesh:
      def __init__(self, name='MBDynNodes', id=0):
        self.name = 'MBDynNodes'
        self.id   = id
    class Field:
      def __init__(self, name='', id=0, data=[]):
        self.name = name 
        self.id   = id   
        self.data = data 


  def __init__(self, mbdInterface, \
               config_file_name = './../precice-config.xml' ):

    self.debug = True

    #> Initialize PreCICE participant dictionary, p for participant
    self.p = { \
               'name':'MBDyn', \
               'mesh':{ \
                  'name':'MBDynNodes', 'id':[], 'node_ids':[], \
                  'nodes':[], 'nnodes':[], 'dim':[] }, \
               'fields':{}  \
             }
    field_dict = { 'name':'field_name', 'id':[], 'data':[],
             'type':'scalar/vector', 'io':'read/write' }
  
    #> Use MBDyn interface, containing info about MBDyn-mbc_py

    self.mbd = mbdInterface
        
    # All the data that can be exchanged with MBDyn
    # TODO#1: move in MBDyn interface
    fieldlist = [ 'Position', 'Velocity', 'Rotation', \
                  'AngularVelocity', 'Force', 'Moment' ]
    typelist  = [ 'vector', 'vector', 'vector', \
                  'vector', 'vector', 'vector'  ]
    iolist    = [ 'write', 'write', 'write', 'write', 'read', \
                  'read' ]
     
    #> === PreCICE interface ====================================
    comm_rank = 0; comm_size = 1   
    
    self.interface = precice.Interface( self.p['name'], \
                                        config_file_name, \
                                        comm_rank, comm_size )

    self.dim = self.interface.get_dimensions()

    self.p['mesh']['id'] = \
              self.interface.get_mesh_id( self.p['mesh']['name'] )

    #> === Fields ===============================================
    # get_data_id
    # TODO#2: improve it, after having redefined the fields that can
    # be used for communications (see TODO#1)
    for i in np.arange( len(fieldlist) ):
      fieldn = fieldlist[i]
      if ( self.interface.has_data( fieldn, self.p['mesh']['id'] ) ):
        field = field_dict.copy()
        field['name'] = fieldn
        field['id']   = self.interface.get_data_id( field['name'], \
                                            self.p['mesh']['id'] )
        field['type'] = self.mbd.data[fieldn]['type']
        field['io']   = self.mbd.data[fieldn]['io']
        self.p['fields'][field['name']] = field

    #> === Mesh =================================================
    #> Nodes: set_mesh_vertices
    # new: read from file the reference configuration needed for
    # PreCICE communication
    self.p['mesh']['nodes'] = self.mbd.refConfigNodes()
    # old: initial configuration as the reference configuration
    self.p['mesh']['nnodes'] = np.size( self.mbd.data['Position']['data'], 0 )
    self.p['mesh']['dim']    = np.size( self.mbd.data['Position']['data'], 1 )
    self.p['mesh']['node_id'] = self.interface.set_mesh_vertices( \
                    self.p['mesh']['id'], self.p['mesh']['nodes'] )

    if ( self.debug ): 
      for i in np.arange(self.mbd.socket.nnodes):
        print('  ', i,': ', self.p['mesh']['node_id'][i],' , ', \
                            self.p['mesh']['nodes'][i,:] )

    #> === Initialize ===========================================
    self.dt_precice = self.interface.initialize()
    self.is_ongoing = self.interface.is_coupling_ongoing()
    if ( self.debug ): 
      print(' interface.is_coupling_ongoing: \033[0;32m ▶ %s' % self.is_ongoing )

    #> === Initialize data ======================================
    cowid = precice.action_write_initial_data()
    if ( self.interface.is_action_required( cowid ) ):
      if ( self.p['fields'][fie]['io'] == 'write' ):
        if ( self.p['fields'][fie]['type'] == 'scalar' ):
          self.interface.write_block_scalar_data( \
                                  self.p['fields'][fie]['id'], \
                                  self.p['mesh']['node_id'],   \
                                  self.mbd.data[fie]['data'] )
        if ( self.p['fields'][fie]['type'] == 'vector' ):
          self.interface.write_block_vector_data( \
                                  self.p['fields'][fie]['id'], \
                                  self.p['mesh']['node_id'],   \
                                  self.mbd.data[fie]['data'] )
      
      self.interface.mark_action_fulfilled( cowic )

    self.interface.initialize_data()


  def runPreCICE(self, dt_set):

    n = self.mbd.socket.nnodes; nd = 3

    cowic = precice.action_write_iteration_checkpoint()
    coric = precice.action_read_iteration_checkpoint()

    dt_precice = self.dt_precice

    force = np.zeros((n, nd))

    t = 0.
    niter = 0.
    is_ongoing = self.interface.is_coupling_ongoing()
    while ( is_ongoing ):
      niter = niter + 1
      if (niter == 1):
        start = time.time()
      print('\033[0m   Iteration ▶ ', niter)
      
      if ( self.interface.is_action_required( cowic ) ):
        pos_t = self.mbd.data['Position']['data'][:,:] 
        vel_t = self.mbd.data['Velocity']['data'][:,:] 
        rot_t = self.mbd.data['Rotation']['data'][:,:] 
        ome_t = self.mbd.data['AngularVelocity']['data'][:,:]
        self.interface.mark_action_fulfilled( cowic )
   
      #> Set MBDyn nodal values of forces and moments
      for i in np.arange(n):
        self.mbd.nodal.n_f[i*nd:(i+1)*nd] = self.mbd.data['Force' ]['data'][i,:] 
        self.mbd.nodal.n_m[i*nd:(i+1)*nd] = self.mbd.data['Moment']['data'][i,:] 
      
      dt = min( dt_set, dt_precice )
    
      #> === Communication with MBDyn ===========================
      if ( self.mbd.nodal.send(False) ):
        break
     
      # Receive data from MBDyn
      if ( self.mbd.nodal.recv() ):
        print('**** break, after nodal.recv() ****'); break
    
      #> Read position and velocity from MBDyn
      self.mbd.data['Position'       ]['data'] = np.reshape( self.mbd.nodal.n_x    , (n, nd)) 
      self.mbd.data['Velocity'       ]['data'] = np.reshape( self.mbd.nodal.n_xp   , (n, nd))
      self.mbd.data['Rotation'       ]['data'] = np.reshape( self.mbd.nodal.n_theta, (n, nd))
      self.mbd.data['AngularVelocity']['data'] = np.reshape( self.mbd.nodal.n_omega, (n, nd))

      #> Write to AeroSolver
      for fie in self.p['fields']:
        if ( self.p['fields'][fie]['io'] == 'write' ):
          if ( self.p['fields'][fie]['type'] == 'scalar' ):
            self.interface.write_block_scalar_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id'],   \
                                    self.mbd.data[fie]['data'] )
          if ( self.p['fields'][fie]['type'] == 'vector' ):
            self.interface.write_block_vector_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id'],   \
                                    self.mbd.data[fie]['data'] )
    
      is_ongoing = self.interface.is_coupling_ongoing()
      dt_precice = self.interface.advance(dt)
    
      #> Receive data from AeroSolver and set nodal.n_f field
      for fie in self.p['fields']:
        if ( self.p['fields'][fie]['io'] == 'read' ):
          if ( self.p['fields'][fie]['type'] == 'scalar' ):
            self.mbd.data[fie]['data'] = \
              self.interface.read_block_scalar_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id']    )
          if ( self.p['fields'][fie]['type'] == 'vector' ):
            self.mbd.data[fie]['data'] = \
              self.interface.read_block_vector_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id']    )

      #> Check convergence: iterate or finalize the timestep
      if ( self.interface.is_action_required( coric ) ): # dt not converged
        self.mbd.data['Position'       ]['data'] = pos_t
        self.mbd.data['Velocity'       ]['data'] = vel_t
        self.mbd.data['Rotation'       ]['data'] = rot_t
        self.mbd.data['AngularVelocity']['data'] = ome_t
        self.interface.mark_action_fulfilled( coric )
      else: # dt converged
        if ( self.mbd.nodal.send(True) ):
          break
        # Receive data from MBDyn
        if ( self.mbd.nodal.recv() ):
          print('**** break, after nodal.recv() ****'); break
        t = t + dt
        niter = 0
        end = time.time()
        print('\033[0m   Elapsed Time:', round((end - start),8), 'sec')
        print('\033[0;32m ------------------------------- ')
        print('\033[0;32m ◀ Simulation Time ▶ ', round(t, 6))
        print('\033[0;32m ------------------------------- \033[0m ')
            
    print(' before finalize() ')
    self.interface.finalize()
    print(' Finalize auxiliary coupling ')
    
    #> ~~~ PreCICE coupling between mbdyn and aero solver ~~~~~~~~~~~   
    self.mbd.nodal.destroy()

        


