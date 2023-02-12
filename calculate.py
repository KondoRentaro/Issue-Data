import math

lis_in_file=[] 

lis_atomname=[]
lis_chainid=[]
lis_num=[]

#CA
lis_ca_x2=[]
lis_ca_y2=[]
lis_ca_z2=[]
lis_ca_chain=[]
lis_ca_num=[]  
lis_ca_atom=[]
lis_ca_dis=[]
lis_ca_p=[]

#CB
lis_cb_x2=[]
lis_cb_y2=[]
lis_cb_z2=[]
lis_cb_chain=[]
lis_cb_num=[]  
lis_cb_atom=[]
lis_cb_dis=[]
lis_cb_p=[]

def main():

    lis_chainid,lis_num = data_lis() 

    print(set(lis_chainid)) 

    chainid=input_chainid()
    flg = 0
    lis_num_new = []
    for chain in lis_chainid:
        if chain == chainid:
            lis_num_new.append(lis_num[flg])        
        flg += 1

    print(set(lis_num_new))
    
    lis_in_file = file_load()

    
    num=input_num()

    for data in lis_in_file:
        if data[0] == "ATOM":
            if data[5] == num:
                if data[4] == chainid:
                    
                    if data[2] == "CA":
                        
                        c_ca_x1=float(data[6])
                        c_ca_y1=float(data[7])
                        c_ca_z1=float(data[8])  
                        
                    if data[2] == "CB": 

                        c_cb_x1=float(data[6])
                        c_cb_y1=float(data[7])
                        c_cb_z1=float(data[8])  
                        
            if not data[4] == chainid:
                if data[2] == "CA":

                    e_ca_x2=float(data[6])
                    e_ca_y2=float(data[7])
                    e_ca_z2=float(data[8])  
                    e_ca_num=data[5]   
                    e_ca_atom=data[2] 
                    e_ca_p=data[3]   
                    e_ca_chain=data[4]   

                    lis_ca_x2.append(e_ca_x2)
                    lis_ca_y2.append(e_ca_y2)
                    lis_ca_z2.append(e_ca_z2)
                    lis_ca_num.append(e_ca_num)  
                    lis_ca_atom.append(e_ca_atom)   
                    lis_ca_p.append(e_ca_p)
                    lis_ca_chain.append(e_ca_chain)
                    
                if data[2] == "CB":    

                    e_cb_x2=float(data[6])
                    e_cb_y2=float(data[7])
                    e_cb_z2=float(data[8])  
                    e_cb_num=data[5]  
                    e_cb_atom=data[2]
                    e_cb_p=data[3]
                    e_cb_chain=data[4]

                    lis_cb_x2.append(e_cb_x2)
                    lis_cb_y2.append(e_cb_y2)
                    lis_cb_z2.append(e_cb_z2)
                    lis_cb_num.append(e_cb_num)   
                    lis_cb_atom.append(e_cb_atom)
                    lis_cb_p.append(e_cb_p)
                    lis_cb_chain.append(e_cb_chain)
                        
                elif data[2] == "CA":  #E鎖のCA原子のGLY
                    if data[3] == "GLY":
                                    
                        gly_ca_x=float(data[6])
                        gly_ca_y=float(data[7])
                        gly_ca_z=float(data[8])
                        gly_ca_num=data[5]  
                        gly_ca_atom=data[2]
                        gly_ca_p=data[3]
                        gly_ca_chain=data[4]

                        lis_cb_x2.append(gly_ca_x)
                        lis_cb_y2.append(gly_ca_y)
                        lis_cb_z2.append(gly_ca_z)
                        lis_cb_num.append(gly_ca_num)   
                        lis_cb_atom.append(gly_ca_atom)
                        lis_cb_p.append(gly_ca_p)    
                        lis_cb_chain.append(gly_ca_chain)

                    l1=len(lis_ca_num)  
                    l2=len(lis_cb_num)

                    ca_cnt=0 
                    cb_cnt=0 

    else:
        for e_ca_len in lis_ca_x2:  
            e_ca_dis=length(c_ca_x1,c_ca_y1,c_ca_z1,lis_ca_x2[ca_cnt],lis_ca_y2[ca_cnt],lis_ca_z2[ca_cnt])
            
            lis_ca_dis.append(e_ca_dis)  
            ca_cnt+=1

            if ca_cnt==(l1+1):   
                break

        for e_cb_len in lis_cb_x2:
            e_cb_dis=length(c_cb_x1,c_cb_y1,c_cb_z1,lis_cb_x2[cb_cnt],lis_cb_y2[cb_cnt],lis_cb_z2[cb_cnt])
            lis_cb_dis.append(e_cb_dis) 
            cb_cnt+=1
        
            if cb_cnt==(l2+1):  
                break

        atom_dis=input_atom_dis()
        atom=input_atom()
    for c_e_ca_cal in lis_ca_dis:
       if c_e_ca_cal <= float(atom_dis):   
           if atom == "CA":
                ca_decimal_point=2
                ca_cal_result=math.floor(c_e_ca_cal * 10 ** ca_decimal_point) / (10 ** ca_decimal_point)
                print(lis_cb_num[lis_ca_dis.index(c_e_ca_cal)],":",lis_ca_atom[lis_ca_dis.index(c_e_ca_cal)],":",lis_ca_p[lis_ca_dis.index(c_e_ca_cal)],":",lis_ca_chain[lis_ca_dis.index(c_e_ca_cal)],":",ca_cal_result)
    for c_e_cb_cal in lis_cb_dis:
        if c_e_cb_cal <= float(atom_dis):
            if atom == "CB":
                cb_decimal_point=2
                cb_cal_result=math.floor(c_e_cb_cal * 10 ** cb_decimal_point) / (10 ** cb_decimal_point)
                print(lis_cb_num[lis_cb_dis.index(c_e_cb_cal)],":",lis_cb_atom[lis_cb_dis.index(c_e_cb_cal)],":",lis_cb_p[lis_cb_dis.index(c_e_cb_cal)],":",lis_cb_chain[lis_cb_dis.index(c_e_cb_cal)],":",cb_cal_result)
    

def file_load(): 
    pdbfile=open('pdb_data/cdk5p35p27.pdb','r')
    for in_file in pdbfile:
        in_file=" ".join(in_file.split())  
        in_file=in_file.split()  
        lis_in_file.append(in_file)
    
    return lis_in_file

def data_lis():
    lis_in_file = file_load()
    for data in lis_in_file:
        if data[0] == "ATOM":
            lis_chainid.append(data[4])
            lis_num.append(data[5]) 

    return lis_chainid,lis_num

def length(x1,y1,z1,x2,y2,z2):
    distance=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    return distance  

def input_chainid():
    chainid=input("chainid: ") 
    return chainid

def input_num():
    num=str(input("residue number: ")) 
    return num

def input_atom_dis():
    atom_dis=input("distance between residues: ")  
    return atom_dis

def input_atom():
    atom=input("atom name(CA or CB): ")  
    return atom

main()