from tkinter import *  # 导入tkinter模块【必要步骤】
from PIL import Image, ImageTk
from  tkinter import filedialog
import tkinter as tk
from alignment import *
#ATCG
substitutionMatrixs = [
    [
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
    ],
    [
        [1,-5,-5,-1],
        [-5,1,-1,-5],
        [-5,-1,1,-5],
        [-1,-5,-5,1]
    ],
    [
        [5,-4,-4,-4],
        [-4,5,-4,-4],
        [-4,-4,5,-4],
        [-4,-4,-4,5]
    ]
]

#将碱基转换为集合下标
replace = {'A':0,'T':1,'C':2,'G':3}




root = Tk()  # 创建主窗口【必要步骤】
# 将该窗口赋值给root变量，方便后续使用
root.title("基因序列比对")
root.geometry('1500x900+200+0')#设置窗口大小及位置

alignment=StringVar()
useMatrix = StringVar()
s1 = StringVar()
s2 = StringVar()
gap = StringVar()
result = StringVar()
runtime =  StringVar()

label = Label(root,text="选择比对方式")
label.place(width=100,height=30,x=100,y=50)
alignmentContent = Spinbox(root,textvariable=alignment,values=("全局比对","局部比对"))
alignmentContent.place(width=80, height=30,x=250,y=50)

mlabel = Label(root,text="选择计分矩阵")
mlabel.place(width=100,height=30,x=100,y=100)
mlabelContent = Spinbox(root,textvariable=useMatrix,values=("DNA等价矩阵","转换-颠换矩阵","BLAST矩阵"))
mlabelContent.place(width=100, height=30,x=250,y=100)

gaplabel = Label(root,text="输入gap值")
gaplabel.place(width=100,height=30,x=100,y=140)
gaplabelContent = Entry(root,textvariable=gap)
gaplabelContent.place(width=100, height=30,x=250,y=140)

s1label = Label(root,text="序列1")
s1label .place(width=512, height=20,x=100,y=180)
s1labelcontent=Entry(root,textvariable=s1)
s1labelcontent.place(width=1024, height=50,x=100,y=210)

s2label = Label(root,text="序列2")
s2label .place(width=512, height=20,x=100,y=280)
s2labelcontent=Entry(root,textvariable=s2)
s2labelcontent.place(width=1024, height=50,x=100,y=310)


resultlabel = Label(root,text="比对结果")
resultlabel .place(width=512, height=20,x=100,y=430)
resultlabelcontent=Text(root)
resultlabelcontent.place(width=1024, height=300,x=100,y=460)

timelabel = Label(root,text="比对用时/s")
timelabel .place(width=512, height=20,x=100,y=780)
timelabelcontent=Entry(root, textvariable=runtime)
timelabelcontent.place(width=80, height=30,x=330,y=810)


def aligment():
    matrixIndex = 0
    if useMatrix.get() == "DNA等价矩阵":
        matrixIndex = 0
    elif useMatrix.get() == "转换-颠换矩阵":
        matrixIndex = 1
    else:
        matrixIndex = 2

    senquence1=s1.get().upper().replace("\n","")
    senquence2=s2.get().upper().replace("\n","")
    gapNum = int(gap.get())
    if alignment.get()=="全局比对":
        result, scoreMatrix, pathMatrix, allPaths, time = GlobalAlignment(replace, substitutionMatrixs[matrixIndex], senquence1, senquence2, gapNum)
        resultlabelcontent.delete('0.0',END)
        resultlabelcontent.insert('1.0', "全部比对所有结果如下:\n"+result)
        runtime.set(str(time))
        #ShowPaths(senquence1, senquence2, scoreMatrix, pathMatrix, allPaths)
        # result.set(GlobalAlignment(replace, substitutionMatrixs[matrixIndex], senquence1, senquence2, gapNum))
    else:
        result, scoreMatrix, pathMatrix, allPaths, time = localAlignment(replace, substitutionMatrixs[matrixIndex],
                                                                    senquence1, senquence2, gapNum)
        resultlabelcontent.delete('0.0', END)
        resultlabelcontent.insert('1.0', "局部比对所有结果如下:\n"+result)
        runtime.set(str(time))
        #ShowPaths(senquence1, senquence2, scoreMatrix, pathMatrix, allPaths)
        # result.set(localAlignment(replace, substitutionMatrixs[matrixIndex], senquence1, senquence2, gapNum))
# 插入一个图片
b = Button(root,relief='flat',activebackground='gray',bg='white',
           overrelief='raised',text='比对',command=aligment)#创建按钮
b.place(width=100,height=30,x=300,y=380)#放置按钮




root.mainloop()  # 主窗口进入消息事件循环【必要步骤】