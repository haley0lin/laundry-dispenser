#!/usr/bin/env python
# coding: utf-8

# In[1]:


from scipy.optimize import linprog
import numpy as np
import math
import sys
from queue import Queue


class ILP():
    def __init__(self, c, A_ub, b_ub, A_eq, b_eq, bounds):
        # 全局参数
        self.LOWER_BOUND = -sys.maxsize
        self.UPPER_BOUND = sys.maxsize
        self.opt_val = None
        self.opt_x = None
        self.Q = Queue()

        # 这些参数在每轮计算中都不会改变
        self.c = -c
        self.A_eq = A_eq
        self.b_eq = b_eq
        self.bounds = bounds

        # 首先计算一下初始问题
        r = linprog(-c, A_ub, b_ub, A_eq, b_eq, bounds)

        # 若最初问题线性不可解
        if not r.success:
            raise ValueError('Not a feasible problem!')

        # 将解和约束参数放入队列
        self.Q.put((r, A_ub, b_ub))

    def solve(self):
        while not self.Q.empty():
            # 取出当前问题
            res, A_ub, b_ub = self.Q.get(block=False)

            # 当前最优值小于总下界，则排除此区域
            if -res.fun < self.LOWER_BOUND:
                continue

            # 若结果 x 中全为整数，则尝试更新全局下界、全局最优值和最优解
            if all(list(map(lambda f: f.is_integer(), res.x))):
                if self.LOWER_BOUND < -res.fun:
                    self.LOWER_BOUND = -res.fun

                if self.opt_val is None or self.opt_val < -res.fun:
                    self.opt_val = -res.fun
                    self.opt_x = res.x

                continue

            # 进行分枝
            else:
                # 寻找 x 中第一个不是整数的，取其下标 idx
                idx = 0
                for i, x in enumerate(res.x):
                    if not x.is_integer():
                        break
                    idx += 1

                # 构建新的约束条件（分割
                new_con1 = np.zeros(A_ub.shape[1])
                new_con1[idx] = -1
                new_con2 = np.zeros(A_ub.shape[1])
                new_con2[idx] = 1
                new_A_ub1 = np.insert(A_ub, A_ub.shape[0], new_con1, axis=0)
                new_A_ub2 = np.insert(A_ub, A_ub.shape[0], new_con2, axis=0)
                new_b_ub1 = np.insert(
                    b_ub, b_ub.shape[0], -math.ceil(res.x[idx]), axis=0)
                new_b_ub2 = np.insert(
                    b_ub, b_ub.shape[0], math.floor(res.x[idx]), axis=0)

                # 将新约束条件加入队列，先加最优值大的那一支
                r1 = linprog(self.c, new_A_ub1, new_b_ub1, self.A_eq,
                             self.b_eq, self.bounds)
                r2 = linprog(self.c, new_A_ub2, new_b_ub2, self.A_eq,
                             self.b_eq, self.bounds)
                if not r1.success and r2.success:
                    self.Q.put((r2, new_A_ub2, new_b_ub2))
                elif not r2.success and r1.success:
                    self.Q.put((r1, new_A_ub1, new_b_ub1))
                elif r1.success and r2.success:
                    if -r1.fun > -r2.fun:
                        self.Q.put((r1, new_A_ub1, new_b_ub1))
                        self.Q.put((r2, new_A_ub2, new_b_ub2))
                    else:
                        self.Q.put((r2, new_A_ub2, new_b_ub2))
                        self.Q.put((r1, new_A_ub1, new_b_ub1))


                        
                        
                        
                        
                        
def caltimes(mass,water):
    
    c = np.array([79, 31])
    A = np.array([[-13, -7], [99, 86]])
    b = np.array([-mass, water])#前面負的是最小要有mass 後面是少於水量（正）
    Aeq = None
    beq = None
    bounds = [(0, None), (0, None)]
    

    solver = ILP(-c, A, b, Aeq, beq, bounds)
    solver.solve()
    time=-solver.opt_val
    times=solver.opt_x
    #print("Result:", solver.opt_val, solver.opt_x)
    ma=times[0]
    mb=times[1]
    wa=88*ma+66*mb# 公升
    #p="您所需要的最小洗衣時間為"+ time +"。"+"\n\n"+"請使用標準模式清洗"+ ma +"次，使用快速模式清洗"+mb+ "次方可達成"
    p="您所需要的最小洗衣時間為" + str(time) +"分鐘。"+"\n\n用水量為"+str(wa)+"公升"+"\n\n"+"請使用標準模式清洗"+ str(ma) +"次，使用快速模式清洗"+str(mb)+ "次，方可達成"
    return p





#import pulp

import tkinter as tk


thsr = tk.Tk()
thsr.title("洗衣快速分配")
thsr.geometry('750x850')

#label

label00=tk.Label(thsr, text="【使用說明】\n\n 本程式可以幫助您計算洗衣之最短時間及使用程序，請輸入您所需要盥洗的衣物數量、重量以及最大用水量。\n\n 提醒您，本程式使用之計算標準為標準洗衣最大洗衣量12公斤；快速洗衣最大洗衣量6公斤。\n\n按照本程式計算之結果進行洗衣程序時請盡量接近最大洗衣量以保證您可以達成最快的洗衣效率。\n\n（若您輸入數字後程式沒有反應則您的最大用水量要求可能無法達成，請更改數據後再次嘗試）",font=("Arial",12)) 
label00.grid(row=0,column=0,columnspan=2,padx=20, pady=20)

label0=tk.Label(thsr, text="請輸入您有幾件極薄衣物，本程式將以0.2公斤計算  \n\n（包含短褲、裙子、睡衣等）",font=("Arial",12)) 
label0.grid(row=1,column=0,padx=20, pady=20)


label2=tk.Label(thsr, text="請輸入您有幾件輕薄衣物，本程式將以0.55公斤計算 \n\n（包含襯衫、輕薄外套、洋裝及女性長褲等）",font=("Arial",12)) 
label2.grid(row=2,column=0,padx=20, pady=20)


label4=tk.Label(thsr, text="請輸入您有幾件厚重衣物，本程式將以0.9公斤計算 \n\n（包含男性長褲、厚重外套、牛仔褲等，不包括西裝外套及皮草大衣等特殊衣物）",font=("Arial",12)) 
label4.grid(row=3,column=0,padx=20, pady=20)

label5=tk.Label(thsr, text="若您已知衣服重量也可直接輸入公斤數 \n\n（此數值會與上面衣物進行加總，若全部以重量計算請於上方衣物數量之輸入格填寫0）",font=("Arial",12)) 
label5.grid(row=4,column=0,padx=20, pady=20)



label6=tk.Label(thsr, text="請輸入最大用水量為幾公升？\n\n（一度等於1000公升，水費會因用水量而浮動，一度水為7~12元，請自行估算水費成本）",font=("Arial",12)) 
label6.grid(row=5,column=0,padx=20, pady=20)

#text
lightvar = tk.StringVar()#inputA=light=0.2
light= tk.Entry(thsr, textvariable=lightvar)
light.grid(row=1, column=1)

middlevar = tk.StringVar()#inputB=middle=0.55
middle = tk.Entry(thsr, textvariable=middlevar)
middle.grid(row=2, column=1)

heavyvar = tk.StringVar()#inputC=heavy=0.9
heavy = tk.Entry(thsr, textvariable=heavyvar)
heavy.grid(row=3, column=1)

weightvar = tk.StringVar()#inputC=heavy=0.9
weight = tk.Entry(thsr, textvariable=weightvar)
weight.grid(row=4, column=1)

watervar = tk.StringVar()#water
water = tk.Entry(thsr, textvariable=watervar)
water.grid(row=5, column=1)




def tp():
    #print("start1:",startvar.get())
    #print("end1:",endvar.get())
   
    light=lightvar.get()
    middle=middlevar.get()
    heavy=heavyvar.get() 
    weight=weightvar.get() 
    water=watervar.get() 

    
    light=light.strip()
    middle=middle.strip()
    heavy=heavy.strip()
    weight=weight.strip()
    water=water.strip()
    
    mass=int(light)*0.2+int(middle)*0.55+int(heavy)*0.9+int(weight)
    
    water=int(water)
    #print("start:",start)
    #print("end:",end)
    
    
    ans=caltimes(mass,water)
    
    ret.set(ans)
    

        
#botton

button=tk.Button(thsr, text="計算最短洗衣時間及次數 ",font=("Arial",12), command = tp)
button.grid(row=6,column=0,columnspan=2 ,padx=20, pady=20)#位置

ret = tk.StringVar() 

ans=tk.Label(thsr, textvariable=ret,font=("Arial",12)) 
ans.grid(row=7,column=0,columnspan=2,padx=20, pady=20)


thsr.mainloop()


# In[ ]:





# In[ ]:





# In[ ]:




