'''
***
Writes a new .pdb file with chain information.
Returns a list containing edited file contents. 
Assumes the .pdb file to be edited does not currently have chain information.
***
* filename: the path to the .pdb file to be edited
* newfilename: the path to the resulting .pdb filename
* cterm: the c-terminal residue number 
***
'''
def add_chainid(filename, newfilename, cterm):
    
    # initialize variables
    contents = [] # hold contents of .pdb file
    new_contents = [] # hold contents of new .pdb file
    next_resnum = None # hold the next residue in the .pdb file
    chain_id = 'A' # hold current chain ID
    i = 0 # iterator
    
    f = open(filename, 'r') # open .pdb file
    for line in f:
        contents.append(line) # append contents of .pdb file to a list
    f.close()
    
    for line in contents:
        if 'ATOM' in line: # check for ATOM line
            line_parts = line.split()
            res_num = line_parts[4] # current residue number
            if int(res_num) % cterm == 0: # check if the residue number is the last in the chain
                parts = contents[i+1].split() # check next atom info
                if 'TER' not in parts: # check that next line is an atom
                    next_resnum = parts[4]
                else:
                    # format the new line with chain ID, append to newcontents
                    if len(line_parts[2]) < 4: 
                        line_format = '{:<4s}{:>7s}  {:<3} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                        formatted_list = line_format.format(*line_parts)
                        new_contents.append(formatted_list)
                    if len(line_parts[2]) == 4:
                        line_format = '{:<4s}{:>7s} {:<4} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                        formatted_list = line_format.format(*line_parts)
                        new_contents.append(formatted_list)
                    if next_resnum != res_num: # check if the next atom is on the next chain
                        chain_id = chr(ord(chain_id) + 1) # increment chain ID
            else:
                # format the new line with chain ID, append to new_contents
                if len(line_parts[2]) < 4:
                    print(line)
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
        else:
            # append remark/comment lines to new_contents
            new_contents.append(line)
        i += 1 # incrementor

    # write new .pdb file
    f = open(newfilename, 'w') 
    for item in new_contents:
        f.write(item)
        f.write('          \n') 
    f.close()

    # return list of edited contents
    return new_contents


'''
***
Renumbers several peptides to uniform residue numbering. 
Returns a list containing edited file contents. 
***
* filename: the path to the .pdb file to be edited
* newfilename: the path to the resulting .pdb filename
* cterm: the c-terminal residue number 
* nterm: the n-terminal residue number
***
'''     
def renumber_pdb(filename, newfilename, nterm, cterm):
    # initialize variables
    contents = [] # hold .pdb file contents
    new_contents = [] # hold edited contents
    num = 1 # hold current residue number
    i = 0 # iterator
    
    # open .pdb, append contents to contents list
    f = open(filename, 'r')
    for line in f:
        contents.append(line)
    f.close()


    for line in contents:
        if 'ATOM' in line: # check for atom line
            line_parts = line.split() # hold line contents
            res_num = line_parts[5] # hold current residue number
            if int(res_num) % cterm == 0: # check if current residue is the last in the chain
                parts = contents[i+1].split() # check next atom 
                if 'TER' not in parts: # check if next atom is atom
                    next_resnum = int(parts[5]) # hold next residue number
                # format the new line with appropriate residue number, append to newcontents
                if len(line_parts[2]) < 4:
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                # check if next atom is on a new peptide
                if int(next_resnum) != int(res_num):
                    num = nterm
            else:
                parts = contents[i+1].split() # check next atom
                if 'TER' not in parts: # check if next atom is atom
                    next_resnum = parts[5] # hold next residue number
                # format the new line with appropriate residue number, append to newcontents
                if len(line_parts[2]) < 4: 
                    print(line)
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if int(next_resnum) != int(res_num):
                    num += 1
        else:
            # append comment lines to new_contents
            new_contents.append(line)
        i += 1 # iterator

    # write new .pdb file
    f = open(newfilename, 'w') 
    for line in new_contents:
        f.write(line)
        f.write('      \n')
    f.close()

    # return list of edited contents
    return new_contents
        

