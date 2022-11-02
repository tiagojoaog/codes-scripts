from ase.io import *
from math import cos, sin, pi
from ase.build import surface
from ase.atoms import *
from ase.constraints import FixAtoms
import numpy as np
import os,sys,ase
import ase.data
import copy,numbers

class myatoms(Atoms): #this is just inheriting Atoms
    def import_from_Atoms(self,atoms): #same def
        new_atoms = self.__class__(cell=atoms.cell, pbc=atoms.pbc, info=atoms.info,
                               celldisp=atoms._celldisp.copy())

        new_atoms.arrays = {}
        for name, a in atoms.arrays.items():
            new_atoms.arrays[name] = a.copy()
        new_atoms.constraints = copy.deepcopy(atoms.constraints)
        return new_atoms

    def rotate(self, a, v, center=(0, 0, 0), rotate_cell=False): #same def
        if not isinstance(a, numbers.Real):
            a, v = v, a

        norm = np.linalg.norm
        v = string2vector(v)

        normv = norm(v)

        if normv == 0.0:
            raise ZeroDivisionError('Cannot rotate: norm(v) == 0')

        if isinstance(a, numbers.Real):
            a *= pi / 180
            v /= normv
            c = cos(a)
            s = sin(a)
        else:
            v2 = string2vector(a)
            v /= normv
            normv2 = np.linalg.norm(v2)
            if normv2 == 0:
                raise ZeroDivisionError('Cannot rotate: norm(a) == 0')
            v2 /= norm(v2)
            c = np.dot(v, v2)
            v = np.cross(v, v2)
            s = norm(v)
            eps = 1e-7
            if s < eps:
                v = np.cross((0, 0, 1), v2)
                if norm(v) < eps:
                    v = np.cross((1, 0, 0), v2)
                assert norm(v) >= eps
            elif s > 0:
                v /= s

        center = self._centering_as_array(center)

        p = self.arrays['positions'] - center
        self.arrays['positions'][:] = (c * p -
                                       np.cross(p, s * v) +
                                       np.outer(np.dot(p, v), (1.0 - c) * v) +
                                       center)
        if rotate_cell:
            rotcell = self.get_cell()
            rotcell[:] = (c * rotcell -
                          np.cross(rotcell, s * v) +
                          np.outer(np.dot(rotcell, v), (1.0 - c) * v))    #TODO: This class is Inherited from Atoms to improve this part here later.
            self.set_cell(rotcell)

class slab(): # My Slab definition to generate surfaces with different terminations
    def __init__(self,atoms,cell_threshold=0.01,bond_threshold=0.3,cubic_check=True,repeat=(1,1,1),miller=(1,0,0),vacuum=18,layers=4,frozen_layers=2):
        self.bulk           = atoms
        self.cell_threshold = cell_threshold
        self.bond_threshold = bond_threshold
        self.miller         = miller
        self.repeat         = repeat
        self.vacuum         = vacuum
        self.layers         = layers
        self.frozen_layers  = frozen_layers
        self.cell_matrix    = self.bulk.get_cell()
        if cubic_check:
            self.check_conventional()       
        slab_init=ase.build.surface(self.bulk, self.miller,1,periodic=True)
        make_slab=ase.build.surface(self.bulk, self.miller,self.layers,vacuum=self.vacuum/2,periodic=True)
        make_extra_layer=ase.build.surface(self.bulk, self.miller,self.layers+1,vacuum=self.vacuum/2,periodic=True)
        N=len(slab_init)
        extra_indices=[i+len(make_slab) for i in range(N)]
        bottom_indices=[i for i in range(N)]
        dummy_slab_dict=self.dict_pos(slab_init)
        connect_dict=self.connect_simple(extra_indices,bottom_indices)
        repeated_slab=make_slab.copy()
        repeated_slab=repeated_slab.repeat(self.repeat)
        repeated_slab.center(vacuum=self.vacuum/2,axis=2)
        self.freeze(repeated_slab,'lower_half')
        repeated_slab.write('slab_original.traj')
        rotated_slab=myatoms()
        rotated_slab=rotated_slab.import_from_Atoms(repeated_slab)
        rotated_slab.rotate(180,'x',center=(0,0,0),rotate_cell=True)
        self.freeze(rotated_slab,'lower_half')
        rotated_slab.write('slab_original_rotated.traj')
        save_iter=1
        print("{:.2f}".format(0.0),' % done.',end='\r')
        for i in range(N):
            index=self.lowest_point(dummy_slab_dict)
            dummy_slab_dict.pop(index)
            extra_layer=make_extra_layer.get_positions()[connect_dict[index]]
            pos=self.pos_list(make_slab,slab_init,index,extra_layer)
            make_slab.set_positions(pos)
            self.is_bonded(make_slab)
            if self.bonded:
                save_slab=make_slab.copy()
                save_slab.center(vacuum=self.vacuum/2,axis=2)
                rotated_slab=myatoms()
                rotated_slab=rotated_slab.import_from_Atoms(save_slab)
                rotated_slab.rotate(180,'x',center=(0,0,0),rotate_cell=True)
                rotated_slab=rotated_slab.repeat(self.repeat)
                self.freeze(rotated_slab,'lower_half')
                #rotated_slab.rotate(-180,'z',center=(0,0,0),rotate_cell=True)
                rotated_slab.write('slab_rotated_termination_'+str(save_iter)+'.traj')
                save_slab=save_slab.repeat(self.repeat)
                self.freeze(save_slab,'lower_half')
                save_slab.write('slab_termination_'+str(save_iter)+'.traj')
                save_iter+=1
            print("{:.2f}".format(100*i/N),' % done.',end='\r')
    def freeze(self,slab,string): #This definition is not ideal for well-ordered crystals such as Cu, where z does not change from atom to atom within one layer. (atoms frozen are not consistent with the free ones from layer to layer)
        properties= {}            #For such crystals, it is perhaps better to freeze by hand. TODO: Improve this.
        elements = []
        count=0
        for i in slab.get_chemical_symbols():
            properties[i+'_'+str(count)]=slab.get_positions()[count][2]
            if i not in elements:
                elements.append(i)
            count+=1
        full_frozen = []
        for i in elements:
            frozen = {}
            for j in properties:
                element=j.split('_')[0]
                index=j.split('_')[1]
                z=properties[j]
                if i==element:
                    frozen[int(index)]=z
            frozen_list=[k for k, v in sorted(frozen.items(), key=lambda item: item[1])]
            if string == 'lower_half':
                frozen_list=frozen_list[0:int(len(frozen)*self.frozen_layers/self.layers)]
            elif string == 'upper_half':
                frozen_list=frozen_list[int(len(frozen)*self.frozen_layers)/self.layers:]
            for k in frozen_list:
                full_frozen.append(k)
        c=FixAtoms(indices=full_frozen)
        slab.set_constraint(c)

    def check_conventional(self):
        a=self.cell_matrix[0]
        b=self.cell_matrix[1]
        c=self.cell_matrix[2]

        if abs(a[0]-np.linalg.norm(a)) > self.cell_threshold:
            self.error('''Cell is not cubic. Conventional cells are often cubic but not always!
                          If you are certain this is a conventional cell, 
                          please use cubic_check=False option.''')
        elif abs(b[1]-np.linalg.norm(b)) > self.cell_threshold:
            self.error('''Cell is not cubic. Conventional cells are often cubic but not always!
                          If you are certain this is a conventional cell, 
                          please use cubic_check=False option.''')
        elif abs(c[2]-np.linalg.norm(c)) > self.cell_threshold:
            self.error('''Cell is not cubic. Conventional cells are often cubic but not always!
                          If you are certain this is a conventional cell, 
                          please use cubic_check=False option.''')

    def error(self,message):
        print('ERROR:',message)
        sys.exit(2)

    def dict_pos(self,atoms):
        index=0
        dictionary={}
        for i in atoms.get_positions():
          dictionary[index]=i
          index+=1
        return dictionary
     
    def pos_list(self,atoms,slab_init,index,extra_layer):
        pos=[]
        iterate=0
        for i in atoms.get_positions():
            if iterate!=index:
                pos.append(i)
            else:
                pos.append(extra_layer)
            iterate+=1
        return pos

    def coordinate_transformation(self,indices,atoms):
        center=np.array((0.0,0.0,0.0))
        for i in indices:
            center+=np.array(atoms.get_positions()[i])
        center=center/len(indices)
        return center

    def connect(self,make_slab,make_extra_layer,extra_indices,bottom_indices):
        connect_indices={}
        center_slab=self.coordinate_transformation(bottom_indices,make_slab)
        center_extra=self.coordinate_transformation(extra_indices,make_extra_layer)
        for i in bottom_indices:
            center_pos1=np.array(make_slab.get_positions()[i])-center_slab
            for j in extra_indices:
                center_pos2=np.array(make_extra_layer.get_positions()[j])-center_extra
                if abs(center_pos2[0]-center_pos1[0])<0.2 and abs(center_pos2[1]-center_pos1[1])<0.2 and abs(center_pos2[2]-center_pos1[2])<0.2:
                    connect_indices[i]=j
        return connect_indices

    def connect_simple(self,extra_indices,bottom_indices):
        connect_indices={}
        for i,j in zip(bottom_indices,extra_indices):
            connect_indices[i]=j
        return connect_indices

    def lowest_point(self,dictionary):
        z=None
        index=None
        for i in dictionary:
            if z==None or dictionary[i][2]<z:
                z=dictionary[i][2]
                index=i
        return index

    def is_bonded(self,atoms):
        self.bonded=True
        bondlist=[]
        limit=len(atoms)
        indices=range(4*limit,5*limit)
        repeated_atoms=atoms.repeat((3, 3, 1))
        for ind1 in indices:
            i=repeated_atoms.get_positions()[ind1]
            atomic_number_i=repeated_atoms.get_atomic_numbers()[ind1]
            bond='no'
            covalent_i = ase.data.covalent_radii[atomic_number_i]
            for ind2 in range(len(repeated_atoms)):
                j=repeated_atoms.get_positions()[ind2]
                atomic_number_j=repeated_atoms.get_atomic_numbers()[ind2]
                covalent_j = ase.data.covalent_radii[atomic_number_j]
                bond_threshold=covalent_j+covalent_i+self.bond_threshold
                bond_length=np.linalg.norm(j-i)
                if bond_length!=0 and bond_length<=bond_threshold:
                    bond='yes'
            bondlist.append(bond)
        if 'no' in bondlist:
            self.bonded=False
                               
