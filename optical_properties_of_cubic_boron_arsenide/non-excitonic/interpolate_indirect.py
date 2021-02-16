import numpy as np
for col in range(8,9):
	with open('indirect_eps2_col{num}.dat'.format(num=col),'r') as fin:
		prev_row = [0.0,0.0]
		dx = 0.01	
		initial_point_total = 100
		interp_point_number = 5
		last_energy = 5.0
		last_target_energy = 40
		decay_factor = 0.75 #smaller is faster
		with open('indirect_eps2_col{num}_interpolated.dat'.format(num=col),'w') as fout:
			fout.write('{x} {y}\n'.format(x=prev_row[0],y=prev_row[1]))
			for row in fin:
				eng,eps2 = [float(num) for num in row.split()]
				xs = np.linspace(prev_row[0]+dx,eng,interp_point_number)
				ys = np.linspace(prev_row[1],eps2,interp_point_number+1)[1:]
				for x,y in zip(xs,ys):
					fout.write('{x} {y}\n'.format(x=x,y=y))
				prev_row = [eng,eps2]	
			for x in np.linspace(last_energy+dx,last_target_energy,last_target_energy/dx-interp_point_number*initial_point_total):
				thisy = prev_row[1]*np.exp(-(x-(last_energy+dx))/decay_factor)
				fout.write('{x} {y}\n'.format(x=x,y=thisy))
