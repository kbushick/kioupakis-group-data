#this is one of main scripts for setting up indirect auger runs. Generate the q list, and nscf points and unique mapping
#the mapping is necessary because given symmetry there are repeated kpoints, and we only keep unique ones for the klist
#note that the order that the points are iterated through (108-119) is important and reflected in the auger source code as well

import numpy as np
import sys
import os

num_arg = len(sys.argv)
if num_arg != 5:
  sys.exit('write_grid.py is expecting command line arguments in the format: eeh/hhe full/irr auger_input cb_vb_pts_dir')

calc_type = str(sys.argv[1]) #eeh/hhe
grid_type = str(sys.argv[2]) #full/irr
auger_fname = str(sys.argv[3]) #usually auger.in
bep_dir = str(sys.argv[4]) #usually $SCRIPTS_DIR/input

#extract material params from the auger input file
recip_lat = np.zeros((3,3))
with open(auger_fname,'r') as fin:
  for row in fin:
    if str(row.split('=')[0].strip()) == 'nq(:)':
      ks = [int(x.strip()) for x in row.split('=')[1].split(',')[:3]]
      print('read q grid as {}x{}x{} from input'.format(ks[0],ks[1],ks[2]))
    elif str(row.split('=')[0].strip()) == 'alat':
      alat = float(row.split('=')[1].split(',')[0])
      print('read alat as {} from input'.format(alat))
    elif str(row.split('=')[0].strip()) == 'b(:,1)':
      recip_lat[0,:] = [float(x.strip()) for x in row.split('=')[1].split(',')[:3]]
    elif str(row.split('=')[0].strip()) == 'b(:,2)':
      recip_lat[1,:] = [float(x.strip()) for x in row.split('=')[1].split(',')[:3]]
    elif str(row.split('=')[0].strip()) == 'b(:,3)':
      recip_lat[2,:] = [float(x.strip()) for x in row.split('=')[1].split(',')[:3]]

tpa = 2*np.pi/alat #we don't actually want to multiply by this because phonon input is in units of 2pi/a
b = np.transpose(recip_lat) #can multiply this by tpa for full coordinates
print('read b matrix as (transposed):')
print(b)

#band extrema points
bep = np.loadtxt(bep_dir+'/cb_vb_points.dat')
num_cb = int(bep[0][0])
cbs = bep[1:num_cb+1]
vbs = bep[num_cb+1:]

if calc_type == 'eeh':
  k1ps = cbs
  k2ps = cbs
  k3ps = vbs
elif calc_type == 'hhe':
  k1ps = vbs
  k2ps = vbs
  k3ps = cbs
else:
  sys.exit('Calc type is not recognized')
  
if grid_type == 'full':
  qlist = np.zeros((int(ks[0]*ks[1]*ks[2]),3))
  qlist_cart = np.zeros((int(ks[0]*ks[1]*ks[2]),3))
  counter = 0
  #generate a full q grid from [0,1) 
  for i1 in range(0,int(ks[0])):
    for i2 in range(0,int(ks[1])):
      for i3 in range(0,int(ks[2])):
        q = np.array([i1/ks[0],i2/ks[1],i3/ks[2]])
        q_cart = np.matmul(b,q)
        qlist[counter,:] = q
        qlist_cart[counter,:] = q_cart
        counter += 1
  num_q = ks[0]*ks[1]*ks[2]
elif grid_type == 'irr':
  try: 
    #take q points from the previously run kirr output
    qw = np.loadtxt('kirr_coord_weights_full.dat',skiprows=1)
    with open('kirr_coord_weights_full.dat', 'r') as fin:
      num_q = int(fin.readline().strip())
    qlist_cart = np.zeros((num_q,3))
    qlist = qw[:,0:-1]
    #weights are important to preserve for the irr grid implementation
    qweights = np.reshape(qw[:,-1],(-1,1))
    for i,q in enumerate(qlist): 
      qlist_cart[i,:] = np.matmul(b,q)
  except: 
    sys.exit('Irr grid requires kirr_coord_weights_full.dat to be present in the root directory')
else: 
  sys.exit('Grid type is not recognized')

if os.path.isdir('phonon_2'):
  skip_dir_make = True
else:
  skip_dir_make = False

max_nscf = 0
for counter in range(1,num_q+1):
  if counter % 50 == 0:
    print('Now working on preparing input for phonon_{i}'.format(i=counter))

  q = qlist[counter-1,:]
  q_cart = qlist_cart[counter-1,:]
  if not skip_dir_make:
    os.mkdir('phonon_{i}'.format(i=counter))
    np.savetxt('phonon_{i}/qpt.dat'.format(i=counter),np.reshape(q,(1,-1)),fmt='%.18f',delimiter=' ')
    np.savetxt('phonon_{i}/qpt_cart.dat'.format(i=counter),np.hstack((np.reshape(q_cart,(1,-1)),np.array([[1]]))),fmt='%.18f',header='1',delimiter=' ',comments='')

  nscf_points = np.zeros((12*k1ps.shape[0]*k2ps.shape[0]*k3ps.shape[0],3))
  i = 0
  for k1p in k1ps:
    for k2p in k2ps:
      for k3p in k3ps:
        #the order of point generation must be consistent throughout the calculation workflow
        nscf_points[12*i+0,:] = k1p #k1 
        nscf_points[12*i+1,:] = k2p #k2
        nscf_points[12*i+2,:] = k3p #k3
        nscf_points[12*i+3,:] = k1p+k2p+q-k3p #k4_abs
        nscf_points[12*i+4,:] = k1p+q #m15_abs
        nscf_points[12*i+5,:] = k2p+q #m26_abs
        nscf_points[12*i+6,:] = k3p-q #m37_abs
        nscf_points[12*i+7,:] = k1p+k2p-k3p #m48
        nscf_points[12*i+8,:] = k1p+k2p-q-k3p #k4_emit
        nscf_points[12*i+9,:] = k1p-q #m15_emit
        nscf_points[12*i+10,:] = k2p-q #m26_emit
        nscf_points[12*i+11,:] = k3p+q #m37_emit
        i += 1

  #keep only the unique points, and keep a mapping back to the full list (inv)
  unique_points, inv = np.unique(nscf_points,axis=0,return_inverse=True)
  umklapp_points = np.zeros(unique_points.shape)
  #for the unique points, find the point which is in the first BZ with the smallest vector
  for i,k in enumerate(unique_points):
    kmin = k
    k_cart = np.matmul(b,k)
    q_min = np.sqrt(np.dot(k_cart,k_cart))
    for u1 in range(-3,4):
      for u2 in range(-3,4):
        for u3 in range(-3,4):
          ktemp = [u1,u2,u3]+k
          ktemp_cart = np.matmul(b,ktemp)
          qtemp = np.sqrt(np.dot(ktemp_cart,ktemp_cart))
          if qtemp < q_min: #or (np.abs(qtemp-q_min) < 1e-8 and np.sum(ktemp) > np.sum(kmin)):
            q_min = qtemp
            kmin = ktemp
    umklapp_points[i,:] = kmin

  #make sure the Umklapp process has not further reduced the number of points
  final, f_inv = np.unique(umklapp_points,axis=0,return_inverse=True)
  #chain the mappings and increment by one for fortran readability of all points (432 for eeh, 72 for hhe for Si)
  final_map = [f_inv[x]+1 for x in inv] #transition to index starting at 1 for fortran
  w = np.ones((final.shape[0],1)) #add weights to the end of the kpoints
  np.savetxt('phonon_{i}/{ct}_nscf_points.dat'.format(i=counter,ct=calc_type),np.hstack((final,w)),fmt='%.18f',delimiter=' ',header='{i}'.format(i=final.shape[0]),comments='')
  np.savetxt('phonon_{i}/{ct}_index_map.dat'.format(i=counter,ct=calc_type),final_map,fmt='%i',delimiter=' ',comments='')
  if final.shape[0] > max_nscf: 
    max_nscf = final.shape[0]

try: 
  with open('qlist.dat','r') as fin:
    qpts, nscf_max = [int(x) for x in fin.readline().strip().split()]
  if nscf_max >= max_nscf:
    print('Previous run had a larger max nscf, no rewrite of qlist necessary')
  else: 
    print('This run has a larger max nscf, rewriting qlist')
    np.savetxt('qlist.dat',qlist,fmt='%.18f',delimiter=' ',header='{i} {j}'.format(i=num_q,j=max_nscf),comments='')
    np.savetxt('qlist_cart.dat',qlist_cart,fmt='%.18f',delimiter=' ',header='{i} {j}'.format(i=num_q,j=max_nscf),comments='')
except:
  print('No existing qlist file found, writing qlists')
  np.savetxt('qlist.dat',qlist,fmt='%.18f',delimiter=' ',header='{i} {j}'.format(i=num_q,j=max_nscf),comments='')
  np.savetxt('qlist_cart.dat',qlist_cart,fmt='%.18f',delimiter=' ',header='{i} {j}'.format(i=num_q,j=max_nscf),comments='')

if grid_type == 'irr':
  np.savetxt('qweights.dat',qweights,fmt='%.18f',delimiter=' ')

