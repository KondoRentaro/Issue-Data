#  pdbfile のポリペプチド鎖 Chain の残基番号 num から 
# residue_dis Å 以内に位置する他のポリペプチド鎖の残基リストと残基間距離
# これが1番

import math

pdbfile=open('3bcc.pdb','r')

lis_model=[] #ファイルのpdbデータをリストに格納

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

chain_id=input("起点とするタンパク質: ") #chainid指定(ない場合あり)
num=str(input("残基番号: ")) #残基番号指定(ない場合あり)
residue_dis=input("距離: ")  #距離指定(ない場合あり)
atom=input("原子名(CAかCB): ")  #CαかCβか指定

def length(x1,y1,z1,x2,y2,z2):
    distance=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    return distance  #三点間の距離の公式 (残基間距離)

# def main():
#     cal()
#     result()

def cal():
    for model in pdbfile:
        model=" ".join(model.split())  #空白をつなげる
        model=model.split()  #要素で区切り1原子ごとにリスト化
        lis_model.append(model)
    
    for data in lis_model:
        if data[0] == "ATOM":
            if data[5] == num:
                if data[4] == chain_id:
                    if data[2] == "CA":
                    
                        ca_x1=float(data[6])
                        ca_y1=float(data[7])
                        ca_z1=float(data[8])  #指定されたデータの座標を格納
                        
                    if data[2] == "CB": 

                        cb_x1=float(data[6])
                        cb_y1=float(data[7])
                        cb_z1=float(data[8])  #C鎖抜き出す作業終わり
                
            if not data[4] == chain_id:
                if data[2] == "CA":

                    ca_x2=float(data[6])
                    ca_y2=float(data[7])
                    ca_z2=float(data[8])  #E(CA)の座標を格納
                    ca_num=data[5]  #Eの残基の番号を抽出  
                    ca_atom=data[2]  #CA
                    ca_p=data[3]   #アミノ酸名
                    ca_chain=data[4]   

                    lis_ca_x2.append(ca_x2)
                    lis_ca_y2.append(ca_y2)
                    lis_ca_z2.append(ca_z2)
                    lis_ca_num.append(ca_num)   #E(CA)の座標、残基の番号をリストに格納
                    lis_ca_atom.append(ca_atom)   
                    lis_ca_p.append(ca_p)
                    lis_ca_chain.append(ca_chain)
                    
                if data[2] == "CB":    #E鎖のCB原子を全て抜き出す

                    cb_x2=float(data[6])
                    cb_y2=float(data[7])
                    cb_z2=float(data[8])  
                    cb_num=data[5]  
                    cb_atom=data[2]
                    cb_p=data[3]
                    cb_chain=data[4]

                    lis_cb_x2.append(cb_x2)
                    lis_cb_y2.append(cb_y2)
                    lis_cb_z2.append(cb_z2)
                    lis_cb_num.append(cb_num)   
                    lis_cb_atom.append(cb_atom)
                    lis_cb_p.append(cb_p)
                    lis_cb_chain.append(cb_chain)
                        
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
                        lis_cb_p.append(gly_ca_p)    #GLYのCAをCBに代用　(これで混ぜるのがよくない??)
                        lis_cb_chain.append(gly_ca_chain)

                    l1=len(lis_ca_num)  #E鎖の長さをl1、l2に入れてる   iが要素以上になったら止まるようにするやつ(E鎖で測ってる)
                    l2=len(lis_cb_num)

                    ca_cnt=0 #Eの座標の1番目から順番にカウント
                    cb_cnt=0 

    else:
        for ca_len in lis_ca_x2:  
            ca_dis=length(ca_x1,ca_y1,ca_z1,lis_ca_x2[ca_cnt],lis_ca_y2[ca_cnt],lis_ca_z2[ca_cnt])
            
            lis_ca_dis.append(ca_dis)  #Cのnum番目とEの距離の結果をリストに格納   (やり方違うかも)
            ca_cnt+=1

            if ca_cnt==(l1+1):  #ca_cntがl1+1の長さと同じになったらやめる  
                break

        for cb_len in lis_cb_x2:
            cb_dis=length(cb_x1,cb_y1,cb_z1,lis_cb_x2[cb_cnt],lis_cb_y2[cb_cnt],lis_cb_z2[cb_cnt])
            lis_cb_dis.append(cb_dis) 
            cb_cnt+=1
        
            if cb_cnt==(l2+1):  
                break

def result():
    for ca_cal in lis_ca_dis:
        if ca_cal <= float(residue_dis):   
            if atom == "CA":
                ca_decimal_point=2
                ca_cal_result=math.floor(ca_cal * 10 ** ca_decimal_point) / (10 ** ca_decimal_point)
                print(lis_cb_num[lis_ca_dis.index(ca_cal)],":",lis_ca_atom[lis_ca_dis.index(ca_cal)],":",lis_ca_p[lis_ca_dis.index(ca_cal)],":",lis_ca_chain[lis_ca_dis.index(ca_cal)],":",ca_cal_result)

    for cb_cal in lis_cb_dis:
        if cb_cal <= float(residue_dis):
            if atom == "CB":
                cb_decimal_point=2
                cb_cal_result=math.floor(cb_cal * 10 ** cb_decimal_point) / (10 ** cb_decimal_point)
                print(lis_cb_num[lis_cb_dis.index(cb_cal)],":",lis_cb_atom[lis_cb_dis.index(cb_cal)],":",lis_cb_p[lis_cb_dis.index(cb_cal)],":",lis_cb_chain[lis_cb_dis.index(cb_cal)],":",cb_cal_result)
        
main()
