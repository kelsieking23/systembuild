import os
import subprocess
import re

class SystemBuild:

    # Initializes paramaters of system
    def __init__(self, nmer, isfragments, iscap, nterm, cterm, gromacs, structure, water, ff, ter, ndx_pdb2gmx, ignh, ndx_solv_ions, bt, box, angles, 
                center, pbc, ions_mdp, pname, nname, neutral, np, nn, conc, minim_mdp, 
                nvt_mdp, npt_mdp, md_mdp):
        
        # GROMACS location 
        self.gromacs = gromacs

        # pdb2gmx
        self.nmer = nmer
        self.structure = structure
        self.water = water
        self.iscap = iscap
        self.ff = ff
        self.ndx_pdb2gmx = ndx_pdb2gmx
        self.ignh = ignh
        self.isfragments = isfragments
        self.nterm = nterm
        self.cterm = cterm
        self.ter = ter
        self.known_aas = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'PHE', 'GLY', 'GLU', 'GLN', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'PRO', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

        # solvation & ions
        self.ndx_solv_ions = ndx_solv_ions
        self.bt = bt # box type
        self.box = box # vector
        self.angles = angles # angles between box vectors
        self.center = center
        self.pbc = pbc
        self.ions_mdp = ions_mdp
        self.pname = pname
        self.nname = nname
        self.neutral = neutral
        self.np = np
        self.nn = nn
        self.conc = conc

        # EM, NVT, NPT, MDrun
        self.minim_mdp = minim_mdp
        self.nvt_mdp = nvt_mdp
        self.npt_mdp = npt_mdp
        self.md_mdp = md_mdp


    # Find GROMACS directory relative to current working directory
    def get_gromacs_location(self):
        gromacs = self.gromacs
        gromacs_parts = gromacs.split('/')
        cwd_parts = os.getcwd().split('/')
        gromacs_version = gromacs_parts[-1]
        i = 0
        same = 0
        for item in cwd_parts:
            if item == gromacs_parts[i]:
                same +=1
            if i < (len(gromacs_parts) - 1):
                i +=1
            else:
                break
        difference = len(cwd_parts) - same
        cd = '../'
        cd = cd * difference 
        gromacs = cd + gromacs_version
        return gromacs

    # isolates fragments
    def make_fragments(self):

        # only execute if one wishes to simulate a fragment of a peptide
        isfragments = self.isfragments
        if isfragments == True:

            # initialize variables
            nterm = self.nterm
            cterm = self.cterm
            #ncap = nterm - 1
            #ccap = cterm + 1
            structure = self.structure

            # atoms is a dictionary that will hold all the atoms in the fragment.
            # its key will be the atom number, the value will be a list all the info
            # from the PDB for the atom. 
            atoms = {}

            # open the pdb
            f = open(structure, 'r')

            # n is an incrementer for later
            n = 0

            #initialize a list to hold pdb file contents
            pdb = []

            # put the contents of the PDB file into a list
            for line in f:
                pdb.append(line)
            f.close()

            # parse pdb contents and isolate fragment
            for line in pdb:
                line_parts = line.split()

                # TER indicates the end of a chain. only want chain A, so break if chain A ends
                if 'TER' in line_parts[0]:
                    break
                # we only want atoms, and not any other line
                if line_parts[0] != 'ATOM':
                    n += 1
                    continue
                # if the line is an atom, and atom is a part of the residue that will become part of the n-terminal cap,
                # collect info about that atom, and add it to the dictionary
                if 'ATOM' in line_parts[0]:
                    if int(line_parts[5]) == nterm + 1:
                        atom_info = []
                        current_atom = line_parts[5]
                        while int(line_parts[5]) <= cterm + 1:
                            temp_line = pdb[n]
                            line_parts = temp_line.split()
                            current_atom = line_parts[5]
                            atom_info.append(line_parts)
                            n += 1
                            next_line = pdb[n]
                            next_line_parts = next_line.split()
                            next_atom = next_line_parts[5]
                            if current_atom != next_atom:
                                if int(current_atom) == cterm + 1:
                                    atoms[int(current_atom)] = atom_info
                                    break
                                else:
                                    atoms[int(current_atom)] = atom_info
                                    atom_info = []
                    else:
                        n +=1
                        continue
            keys = []
            for key in atoms:
                if key > cterm:
                    keys.append(key)
            for key in keys:
                del atoms[key]
            # return a dictionary containing all atom info 
            return atoms

    # returns header info on pdb file
    def get_structure_header(self, filename):

        header = []
        f = open(filename, 'r')
        for line in f:
            line_parts = line.split()
            if 'ATOM' in line_parts[0]:
                continue
            elif 'HETATM' in line_parts[0]:
                continue
            if 'MODEL' in line_parts[0]:
                break
            else:
                header.append(line)
        f.close()
        return header

    # caps fragments
    def cap(self):

        # get the fragment dictionary
        atoms = self.make_fragments()

        # only execute if atoms dict exists
        if atoms:
            #initialize variables
            nterm = self.nterm
            cterm = self.cterm
            ncap = nterm - 1
            ccap = cterm + 1

            # specify ACE
            ncap_atoms = atoms[ncap]
            ace_cap = []
            for atom in ncap_atoms:
                #if atom[2] == 'CA':
                    #atom[3] = 'ACE'
                    #ace_cap.append(atom)
                if atom[2] == 'C':
                    atom[3] = 'ACE'
                    ace_cap.append(atom)
                if atom[2] == 'O':
                    atom[3] = 'ACE'
                    ace_cap.append(atom)
            atoms[ncap] = ace_cap

            # specify NH2
            ccap_atoms = atoms[ccap]
            nh2_cap = []
            for atom in ccap_atoms:
                if atom[2] == 'N':
                    atom[3] = 'NME'
                    nh2_cap.append(atom)
            atoms[ccap] = nh2_cap
            # return new capped fragment
            return atoms

    # write a pdb file with fragment
    def write_fragment_pdb(self, filename):

        # initialize variables
        if self.iscap == True:
            atoms = self.cap()
        else:
            atoms = self.make_fragments()
        nterm = self.nterm
        cterm = self.cterm
        #ncap = nterm - 1
        #ccap = cterm + 1
        header = self.get_structure_header(filename)
        structure = self.structure
        structure = structure[0:-4]
        new_pdb_file = os.getcwd() + "/" + structure + "_fragment.pdb"
        f = open(new_pdb_file, 'w')
        for line in header:
            f.write(line)
            f.write('\n')
        f.write('MODEL 1')
        f.write('\n')
        for i in range(nterm, cterm+1):
            residue = atoms[i]
            for item in residue:
                string = 'ATOM'
                # place serial number
                serial = item[1]
                serial_slots = 7 - len(serial)
                spaces = ' ' * serial_slots
                string += spaces + str(serial) + ' '
                # place atom name
                atom_name = item[2]
                atom_name_slots = 4 - len(atom_name)
                spaces = ' ' * atom_name_slots
                string += str(atom_name) + spaces + ' '
                # place res name
                res_name = item[3]
                res_name_slots = 3 - len(res_name)
                spaces = ' ' * res_name_slots
                string += spaces + str(res_name) + ' '
                # place chain id
                chain_id = item[4]
                string += chain_id
                # place res sequence #
                res_seq = item[5]
                res_seq_slots = 4 - len(res_seq)
                spaces = ' ' * res_seq_slots
                string += spaces + res_seq + '    '
                # x val
                xval = item[6]
                xval_slots = 8 - len(xval)
                spaces = ' ' * xval_slots
                string += spaces + xval 
                # y val 
                yval = item[7]
                yval_slots = 8 - len(yval)
                spaces = ' ' * yval_slots
                string += spaces + yval
                # z val
                zval = item[8]
                zval_slots = 8-len(zval)
                spaces = ' ' * zval_slots
                string += spaces + zval
                # place occupancy 
                occupancy = item[9]
                occupancy_slots = 6 - len(occupancy)
                spaces = ' ' * occupancy_slots
                string +=spaces + occupancy
                # place t factor
                tfactor = item[10]
                tfactor_slots = 6 - len(tfactor)
                spaces = ' ' * tfactor_slots
                string += spaces + tfactor + '          '
                # place symbol
                symbol = item[11]
                symbol_slots = 2 - len(symbol)
                spaces = ' ' * symbol_slots
                string += spaces + symbol
                f.write(string)
                f.write('\n')
        f.write('END')
        f.close()
        return new_pdb_file

    def get_ff_index(self):
        ff = self.ff
        gromacs = self.get_gromacs_location()
        ffs = gromacs + '/share/top'
        directory_files = os.listdir(ffs)
        forcefields = []
        for item in directory_files:
            if item[-2:] == 'ff':
                forcefields.append(item)
        forcefields.sort()
        n = 1
        for item in forcefields:
            if ff in item:
                print(n)
                return n
            n +=1 

    def get_water_index(self):
        water = self.water
        ff = self.ff
        gromacs = self.get_gromacs_location()
        ffs = gromacs + '/share/top'
        directory_files = os.listdir(ffs)
        for item in directory_files:
            if ff in item:
                forcefield = item
        ffs = ffs + '/' +forcefield
        ff_files = os.listdir(ffs)
        for item in ff_files:
            if '.dat' in item:
                water_models = item
        water_models = ffs + '/' + water_models
        f = open(water_models, 'r')
        available_water_models = []
        for line in f:
            if line == '\n':
                continue
            else:
                line_parts = line.split()
                available_water_models.append(line_parts[0])
        f.close()
        n = 1
        for item in available_water_models:
            if water in item:
                return n
            n +=1
    
    def separate_chains(self, nmer):
        structure = self.structure
        nterm = self.nterm
        cterm = self.cterm
        ncap = nterm - 1
        ccap = cterm + 1
        len_peptide = ccap - ncap
        if nmer == 6:
            filename = structure[0:-4] + '_hexamer.gro'
        f = open(filename, 'r')
        file_contents = []
        #reses = []
        for line in f:
            file_contents.append(line)
        f.close()
        renumbers = []
        n = 0
        reses = {}
        for item in file_contents:
            if item == file_contents[0]:
                n +=1
                continue
            elif item == file_contents[1]:
                n +=1
                continue
            elif item == file_contents[-1]:
                n +=1
                continue
            else:
                item_parts = item.split()
                resnum = item_parts[0]
                resnum = resnum.strip()
                resnum = resnum[0:-3]
                if int(resnum) > ccap:
                    newncap = resnum
                    newccap = int(resnum) + int(len_peptide)
                    break
                elif resnum not in reses:
                    reses[resnum] = 0
                else:
                    renumbers.append(resnum)
                    n +=1
                    continue
        for item in file_contents:
            if item == file_contents[0]:
                n +=1
                continue
            elif item == file_contents[1]:
                n +=1
                continue
            elif item == file_contents[-1]:
                n +=1
                continue
            else:
                item_parts = item.split()
                resnum = item_parts[0]
                resnum = resnum.strip()
                resnum = resnum[0:-3]
                if int(resnum) > ccap:
                    break
                else:
                    reses[resnum] += 1
        filename = filename[0:-4] + '_renumber.gro'
        f = open(filename, 'w')
        i = 0
        k = 0
        f.write(file_contents[k])

        f.write(file_contents[k+1])

        keys = []
        for item in reses.keys():
            keys.append(item)
        while i < nmer:
            for val in reses:
                print(val)
                n = 1
                l = 0
                while n < reses[val]:
                    line = file_contents[k+2]
                    line_parts = line.split()
                    resinfo = line_parts[0]
                    resinfo = resinfo[-3:]
                    new_resinfo = str(keys[l]) + resinfo
                    newline = new_resinfo + line[8:]
                    f.write(newline)
                    n +=1
                    k +=1
                l +=1
            i +=1
        f.close()

    def check_mindist(self, cutoff):

        seeds_directory = os.getcwd() + '/../structure/seeds/'
        all_mindist_directories = []
        seeds = os.listdir(seeds_directory)
        for item in seeds:
            if '.DS_Store' not in item:
                mindist = seeds_directory + item + '/mindist/'
                all_mindist_directories.append(mindist)
            else:
                continue
        for item in all_mindist_directories:
            files = os.listdir(item)
            print('Checking ' + item + '...')
            rejected = 0
            accepted = 0
            for name in files:
                filename = item + name
                f = open(filename, 'r')
                for line in f:
                    line_parts = line.split()
                    if '@' in line_parts[0]:
                        continue
                    elif '#' in line_parts[0]:
                        continue
                    else:
                        if float(line_parts[1]) < 2:
                            print('Rejected: ' + line_parts[1])
                            rejected +=1
                        else:
                            print('Accepted: ' + line_parts[1])
                            accepted +=1
                f.close()
            print('Total accepted: ' + str(accepted))
            print('Total rejected: ' + str(rejected))

    def update_top_insert_molecules(self):
        nmer = self.nmer
        seeds_directory = os.getcwd() + '/../structure/seeds/'
        all_top_directories = []
        seeds = os.listdir(seeds_directory)
        for item in seeds:
            if '.DS_Store' in item:
                continue
            else:
                top = seeds_directory + item + '/topol.top'
                all_top_directories.append(top)
        for item in all_top_directories:
            contents = []
            n = 0
            if 'topol.top' in item:
                print(item)
                f = open(item, 'r')
            else:
                continue
            for line in f:
                contents.append(line)
            f.close()
            for line in contents:
                if '[ molecules ]' not in line:
                    n +=1
                    continue
                else:
                    mols_data = contents[n+2]
                    mols_data_split = mols_data.split()
                    nmols = mols_data_split[1]
                    nmols_stripped = nmols.strip()
                    if int(nmols_stripped) != nmer:
                        fixed_line = mols_data_split[0] + '     ' + str(nmer)
                        contents[n+2] = fixed_line
                        filename = item[0:-4] + '_' + str(nmer) + '.top'
                        f = open(filename,'w')
                        for line in contents:
                            f.write(line)
                        f.close()
                        cmd = 'mv ' + item + ' ' + item[0:-4] + '_old.top'
                        os.system(cmd)
                        cmd = 'mv ' + filename + ' ' + item
                        os.system(cmd)
                    else:
                        continue
    def get_atoms_pdb(self, filename):
        atoms = {}
        f = open(filename, 'r')
        for line in f:
            line_parts = line.split()
            if 'ATOM' not in line_parts[0]:
                continue
            else:
                info = []
                for part in line_parts:
                    info.append(part)
                atom_name = line_parts[1]
                atoms[str(atom_name)] = info
        f.close()
        print(atoms)
        return atoms
    
    def get_atoms_gro(self, filename):
        
        atoms = {}
        f = open(filename, 'r')
        n = 0
        for line in f:
            if n <2:
                n +=1
                continue
            line_parts = line.split()
            if len(line_parts) < 4:
                continue
            if 'TIP3' in line:
                break
            print(line)
            info = []
            line_parts = line.split()
            for part in line_parts:
                info.append(part)
            atom_name = line_parts[0][-3:]
            atoms[str(atom_name)] = info
            
        f.close()
        #nmer = self.nmer
        #peplength = int(atom_name)
        #print(n/nmer)
        return atoms

    
    def convert_pdb_to_gro(self, filename, newfilename, box_vectors):
        atoms = self.get_atoms(filename)
        keys = []
        for key in atoms.keys():
            keys.append(key)
        values = []
        for value in atoms.values():
            values.append(value)
        f = open(newfilename, 'w')
        f.write('GENERATED BY CONVERT_PDB_TO_GRO FROM SYSTEM BUILD\n')
        f.write(str(len(keys)) + '\n')
        n = 0
        for value in values:
            space = ' '
            res_num = value[5]
            diff_res_num = 5 - len(res_num)
            spaces_res_num = space * diff_res_num

            res_name = value[3]
            diff_res_name = 5 - len(res_name)
            spaces_res_name = space * diff_res_name
            
            atom_name = value[2]
            diff_atom_name = 5 - len(atom_name)
            spaces_atom_name = space * diff_atom_name 
            
            atom_num = value[1]
            diff_atom_num = 5 - len(atom_num)
            spaces_atom_num = space * diff_atom_num

            x_position_float = float(value[6]) / 10
            x_position_rounded = str(round(x_position_float, 3))
            x_pos_parts = x_position_rounded.split('.')
            if len(x_pos_parts[1]) == 2:
                x_position_rounded = x_position_rounded + '0'
            elif len(x_pos_parts[1]) == 1:
                x_position_rounded = x_position_rounded + '00'
            diff_x_position = 8 - len(x_position_rounded)
            spaces_x_pos = space * diff_x_position
            
            y_position_float = float(value[7]) / 10
            y_position_rounded = str(round(y_position_float, 3))
            y_pos_parts = y_position_rounded.split('.')
            if len(y_pos_parts[1]) == 2:
                y_position_rounded = y_position_rounded + '0'
            elif len(y_pos_parts[1]) == 1:
                y_position_rounded = y_position_rounded + '00'
            diff_y_position = 8 - len(y_position_rounded)
            spaces_y_pos = space * diff_y_position
            
            z_position_float = float(value[8]) / 10
            z_position_rounded = str(round(z_position_float, 3))
            z_pos_parts = z_position_rounded.split('.')
            if len(z_pos_parts[1]) == 2:
                z_position_rounded = z_position_rounded + '0'
            elif len(z_pos_parts[1]) == 1:
                z_position_rounded = z_position_rounded + '00'
            diff_z_position = 8 - len(z_position_rounded)
            spaces_z_pos = space * diff_z_position

            string = spaces_res_num + res_num + res_name + spaces_res_name + spaces_atom_name + atom_name + spaces_atom_num + atom_num + spaces_x_pos + x_position_rounded + spaces_y_pos + y_position_rounded + spaces_z_pos + z_position_rounded
            f.write(string)
            f.write('\n')
        x_vector = box_vectors[0]
        y_vector = box_vectors[1]
        z_vector = box_vectors[2]
        f.write(str(x_vector) + ' ' + str(y_vector) + ' ' + str(z_vector))
            

        f.close()
        

    def renumber_chains(self, filename, newfilename):
        atoms = self.get_atoms(filename)
        values = atoms.values()
        header = self.get_structure_header(filename)
        f = open(newfilename, 'w')
        for line in header:
            f.write(line)
        f.write('MODEL 1')
        f.write('\n')

        for value in values:
            if 'REMARK' in value:
                continue
            elif 'END' in value:
                continue
            elif 'TER' in value:
                continue
            else:
                if 'PROA' in value:
                    value[4] = 'A'
                if 'PROB' in value:
                    value[4] = 'B'
                if 'PROC' in value:
                    value[4] = 'C'
                if 'PROD' in value:
                    value[4] = 'D'
                if 'PROF' in value:
                    value[4] = 'F'
                data = value
                string = self.write_strings(data)
                f.write(string)
                f.write('\n')
        f.write('END')
        f.close()

        
        

    def write_strings(self, data):
        string = 'ATOM'
        # place serial number
        serial = data[1]
        serial_slots = 7 - len(serial)
        spaces = ' ' * serial_slots
        string += spaces + str(serial) + ' '
        # place atom name
        atom_name = data[2]
        atom_name_slots = 4 - len(atom_name)
        spaces = ' ' * atom_name_slots
        string += str(atom_name) + spaces + ' '
        # place res name
        res_name = data[3]
        res_name_slots = 3 - len(res_name)
        spaces = ' ' * res_name_slots
        string += spaces + str(res_name) + ' '
        # place chain id
        chain_id = data[4]
        string += chain_id
        # place res sequence #
        res_seq = data[5]
        res_seq_slots = 4 - len(res_seq)
        spaces = ' ' * res_seq_slots
        string += spaces + res_seq + '    '
        # x val
        xval = data[6]
        xval_slots = 8 - len(xval)
        spaces = ' ' * xval_slots
        string += spaces + xval 
        # y val 
        yval = data[7]
        yval_slots = 8 - len(yval)
        spaces = ' ' * yval_slots
        string += spaces + yval
        # z val
        zval = data[8]
        zval_slots = 8-len(zval)
        spaces = ' ' * zval_slots
        string += spaces + zval
        # place occupancy 
        occupancy = data[9]
        occupancy_slots = 6 - len(occupancy)
        spaces = ' ' * occupancy_slots
        string +=spaces + occupancy
        # place t factor
        tfactor = data[10]
        tfactor_slots = 6 - len(tfactor)
        spaces = ' ' * tfactor_slots
        string += spaces + tfactor + '      '
        # place symbol
        symbol = data[11]
        symbol_slots = 2 - len(symbol)
        spaces = ' ' * symbol_slots
        string += spaces + symbol
        print(string)
        return string


    def renumber_gro_separate_peptides(self, ndx, filename):
        i = 0
        index = []
        f = open(ndx, 'r')
        for line in f:
            index.append(line)
        f.close

        groups = []
        i = 0
        group = []
        for item in index:
            if '[' in item:
                if i > 1:
                    groups.append(group)
                    group = []
                else:
                    group = []
            else:
                item_parts = item.split()
                for part in item_parts:
                    group.append(part.strip())
                i +=1
        groups.append(group)

        f = open(filename, 'r')
        contents = []
        for line in f:
            contents.append(line)
        f.close()

        edited_contents = []
        i = 0
        it = 0
        groupnum = 0
        seq = []
        for line in contents:
            limit = len(contents) - 1
            if it < 2:
                it +=1
                continue
            elif line == contents[limit]:
                continue
            else:
                line_parts = line.split()
                res_info = line_parts[0].strip()
                atom_num = line_parts[2].strip()
                group1 = groups[0]
                if atom_num == group1[i]:
                    pair = (res_info, group1[i])
                    seq.append(pair)
                    i +=1
                limit = len(group1) -1
                if i > limit:
                    break 
        edited_contents = []
        groupnum +=1
        i = 0
        it = 0
        for line in contents:
            limit = len(contents) - 1
            if i == limit:
                edited_contents.append(line.strip('\n'))
                break
            if 'TIP3' in line:
                break
            if it < 2:
                edited_contents.append(line.strip('\n'))
                it +=1
                continue
            elif line == contents[limit]:
                edited_contents.append(line.strip('\n'))
                continue
            else:
                stripped = line.strip()
                constant_data = stripped[6:]
                newres = seq[i][0]
                check = seq[i][1]
                newres = '   ' + newres 
                newline = newres + constant_data
                edited_contents.append(newline)
                if check == '305':
                    i = 0
                else:
                    i +=1
        f = open('check.gro','w')
        for item in edited_contents:
            f.write(item)
            f.write('\n')
        f.close()

        solvent = []
        f = open('../structure/seeds/build2/hexamer_translated_solvate.gro', 'r')
        for line in f:
            if 'TIP3' in line:
                solvent.append(line)
            if 'SOD' in line:
                solvent.append(line)
            if 'CLA' in line:
                solvent.append(line)
        f.close()
        
        f = open('../structure/seeds/build2/check.gro','w')
        for item in edited_contents:
            f.write(item)
            f.write('\n')
        for item in solvent:
            f.write(item)
        f.close()
            




            
    


    
#system = SystemBuild(6, True, False, 16, 35, '/Users/kelsieking/Desktop/gromacs-2019.4', '../structure/2mxu.pdb', 'tip3p', 'charmm36', True, None, True, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None)
