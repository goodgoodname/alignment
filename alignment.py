import numpy as np
import matplotlib.pyplot as plt
import copy
import time


def GlobalScoreMatrix(replace,substitutionMatrix,senquence1,senquence2,gap):
    m = len(senquence2)+1
    n = len(senquence1)+1
    scoreMatrix =np.zeros((m,n))
    pathMatrix = np.zeros((m,n,3))
    scoreMatrix[0][0] = 0
    for i in range(1,m):
        scoreMatrix[i][0] = scoreMatrix[i-1][0] + gap
        pathMatrix[i][0] = [1,0,0]
    for j in range(1,n):
        scoreMatrix[0][j] = scoreMatrix[0][j-1] + gap
        pathMatrix[0][j] = [0, 0, 1]
    for i in range(1,m):
        for j in range(1,n):
            if senquence2[i-1]==senquence1[j-1]:
                scoreMatrix[i][j] = scoreMatrix[i-1][j-1] + substitutionMatrix[replace[senquence2[i - 1]]][replace[senquence1[j - 1]]]
                pathMatrix[i][j] = [0,1,0]
            else:
                tri = scoreMatrix[i - 1][j - 1] + substitutionMatrix[replace[senquence2[i - 1]]][replace[senquence1[j - 1]]]
                up = scoreMatrix[i - 1][j] + gap
                left = scoreMatrix[i][j - 1] + gap
                score = max(tri, up, left)
                scoreMatrix[i][j] = score
                if tri==score:
                    pathMatrix[i][j][1] = 1
                if up == score:
                    pathMatrix[i][j][0] = 1
                if left == score:
                    pathMatrix[i][j][2] = 1
    return scoreMatrix,pathMatrix

#创建局部比对得分矩阵
def LocalScoreMatrix(replace,substitutionMatrix,senquence1,senquence2,gap):
    m = len(senquence2)+1
    n = len(senquence1)+1
    scoreMatrix =np.zeros((m,n))#局部矩阵第一行及第一列均为0，不需要再初始化
    pathMatrix = np.zeros((m,n,3))
    scoreMatrix[0][0] = 0
    for i in range(1,m):
        for j in range(1,n):
            if senquence2[i-1]==senquence1[j-1]:
                scoreMatrix[i][j] = scoreMatrix[i-1][j-1] + substitutionMatrix[replace[senquence2[i - 1]]][replace[senquence1[j - 1]]]
                pathMatrix[i][j] = [0,1,0]
            else:
                tri = scoreMatrix[i - 1][j - 1] + substitutionMatrix[replace[senquence2[i - 1]]][replace[senquence1[j - 1]]]
                up = scoreMatrix[i - 1][j] + gap
                left = scoreMatrix[i][j - 1] + gap
                zero = 0
                score = max(tri, up, left, zero) #获取最大值,与全局比对不同之处在于选取最大值时将0加入选择
                scoreMatrix[i][j] = score
                #记录最大值来的方向，若最大值为0则不必记录
                if tri == score:
                    pathMatrix[i][j][1] = 1
                if up == score:
                    pathMatrix[i][j][0] = 1
                if left == score:
                    pathMatrix[i][j][2] = 1
    return scoreMatrix, pathMatrix

#创建局部比对得分矩阵
def oldLocalScoreMatrix(replace,substitutionMatrix,senquence1,senquence2,gap):
    m = len(senquence2)+1
    n = len(senquence1)+1
    scoreMatrix =np.zeros((m,n))#局部矩阵第一行及第一列均为0，不需要再初始化
    pathMatrix = np.zeros((m,n,3))
    scoreMatrix[0][0] = 0
    for i in range(1,m):
        for j in range(1,n):
            tri = scoreMatrix[i - 1][j - 1] + substitutionMatrix[replace[senquence2[i - 1]]][replace[senquence1[j - 1]]]
            up = scoreMatrix[i - 1][j] + gap
            left = scoreMatrix[i][j - 1] + gap
            zero = 0
            score = max(tri, up, left, zero) #获取最大值,与全局比对不同之处在于选取最大值时将0加入选择
            scoreMatrix[i][j] = score
            #记录最大值来的方向，若最大值为0则不必记录
            if tri == score:
                pathMatrix[i][j][1] = 1
            if up == score:
                pathMatrix[i][j][0] = 1
            if left == score:
                pathMatrix[i][j][2] = 1
    return scoreMatrix, pathMatrix


def oldGlobalScoreMatrix(replace,substitutionMatrix,senquence1,senquence2,gap):
    m = len(senquence2)+1
    n = len(senquence1)+1
    scoreMatrix =np.zeros((m,n))
    pathMatrix = np.zeros((m,n,3))
    scoreMatrix[0][0] = 0
    for i in range(1,m):
        scoreMatrix[i][0] = scoreMatrix[i-1][0] + gap
        pathMatrix[i][0] = [1,0,0]
    for j in range(1,n):
        scoreMatrix[0][j] = scoreMatrix[0][j-1] + gap
        pathMatrix[0][j] = [0, 0, 1]
    for i in range(1,m):
        for j in range(1,n):
            tri = scoreMatrix[i - 1][j - 1] + substitutionMatrix[replace[senquence2[i - 1]]][replace[senquence1[j - 1]]]
            up = scoreMatrix[i - 1][j] + gap
            left = scoreMatrix[i][j - 1] + gap
            score = max(tri, up, left)
            scoreMatrix[i][j] = score
            if tri==score:
                pathMatrix[i][j][1] = 1
            if up == score:
                pathMatrix[i][j][0] = 1
            if left == score:
                pathMatrix[i][j][2] = 1
    return scoreMatrix,pathMatrix

# 判断是否为最终路径中的点
def judgePath(point, LastPath):
    for aPath in LastPath:
        if point in aPath:
            return True
    return False


def ShowPaths(senquence1, senquence2, scoreMatrix,pathMatrix, allPaths):
    s1 = "0" + senquence1
    s2 = "0" + senquence2
    col = list(s1)
    row = list(s2)
    vals = np.array(scoreMatrix)
    # 设置画布大小
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[], )
    # 画出路径矩阵表格 涉及问题：表格与画布大小不一致 解决方式：bbox参数
    the_table = plt.table(cellText=vals, rowLabels=row, colLabels=col, rowLoc='right',
                          loc='center', cellLoc='center', bbox=[0, 0, 1, 1])
    # 设置表格文本字体大小
    the_table.set_fontsize(8)
    m = scoreMatrix.shape[0]
    n = scoreMatrix.shape[1]
    # 画出每个点的路径图
    for i in range(m):
        for j in range(n):
            if pathMatrix[i][j][0]==1:
                if judgePath((i, j), allPaths):  # 若某点在在最终路径中
                    # 画出红色箭头
                    plt.annotate('', xy=((2 * j + 1) / (2 * n), (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 xytext=((2 * j + 1) / (2 * n), (m - i) / (m + 1)),
                                 arrowprops=dict(arrowstyle="<-", color='r', connectionstyle="arc3"))
                else:
                    plt.annotate('', xy=((2 * j + 1) / (2 * n), (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 xytext=((2 * j + 1) / (2 * n), (m - i) / (m + 1)),
                                 arrowprops=dict(arrowstyle="<-", connectionstyle="arc3"))
            if pathMatrix[i][j][1]==1:
                if judgePath((i, j), allPaths):  # 若某点在在最终路径中
                    # 画出红色箭头
                    plt.annotate('', xy=((2 * j + 1) / (2 * n), (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 xytext=(j / n, (m - i) / (m + 1)),
                                 arrowprops=dict(arrowstyle="<-", color='r', connectionstyle="arc3"))
                else:
                    plt.annotate('', xy=((2 * j + 1) / (2 * n), (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 xytext=(j / n, (m - i) / (m + 1)),
                                 arrowprops=dict(arrowstyle="<-", connectionstyle="arc3"))
            if pathMatrix[i][j][2]==1:
                if judgePath((i, j), allPaths):  # 若某点在在最终路径中
                    # 画出红色箭头
                    plt.annotate('', xy=(j / n, (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 xytext=((2 * j + 1) / (2 * n), (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 arrowprops=dict(arrowstyle="->", color='r', connectionstyle="arc3"))
                else:
                    plt.annotate('', xy=(j / n, (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 xytext=((2 * j + 1) / (2 * n), (2 * m - 2 * i - 1) / (2 * (m + 1))),
                                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))



    plt.show()

#寻找全局路径
def FindGlobalPath(i,j,pathMatrix,OnePath,AllPaths):
    #递归终止条件：回到原点（0，0）
    if i == 0 and j == 0:
        OnePath.append((i, j))
        #将OnePath进行深拷贝再加入至最终路径列表AllGlobalPath中
        AllPaths.append(copy.deepcopy(OnePath))
        #将该点出栈
        OnePath.pop()
    else:
        OnePath.append((i, j))
        if pathMatrix[i][j][0]==1:
            # 进行递归记录下一个点
            FindGlobalPath(i-1, j, pathMatrix, OnePath, AllPaths)
        if pathMatrix[i][j][1]==1:
            # 进行递归记录下一个点
            FindGlobalPath(i-1, j-1, pathMatrix, OnePath, AllPaths)
        if pathMatrix[i][j][2]==1:
            # 进行递归记录下一个点
            FindGlobalPath(i, j-1, pathMatrix, OnePath, AllPaths)
        OnePath.pop()

#寻找局部路径
def FindLocalPath(i,j,pathMatrix,OnePath,AllPaths):
    #递归终止条件：回到原点（0，0）
    if pathMatrix[i][j][0]==0 and pathMatrix[i][j][1]==0 and pathMatrix[i][j][2]==0:
        OnePath.append((i, j))
        #将OnePath进行深拷贝再加入至最终路径列表AllGlobalPath中
        AllPaths.append(copy.deepcopy(OnePath))
        #将该点出栈
        OnePath.pop()

    else:
        OnePath.append((i, j))
        if pathMatrix[i][j][0]==1:
            # 进行递归记录下一个点
            FindLocalPath(i-1, j, pathMatrix, OnePath, AllPaths)
        if pathMatrix[i][j][1]==1:
            # 进行递归记录下一个点
            FindLocalPath(i-1, j-1, pathMatrix, OnePath, AllPaths)
        if pathMatrix[i][j][2]==1:
            # 进行递归记录下一个点
            FindLocalPath(i, j-1, pathMatrix, OnePath, AllPaths)
        OnePath.pop()


#输出比对后的序列
def ShowContrastResult(LastPath,senquence1,senquence2):
    #依次输出每条路径
    result = ""
    for k,aPath in enumerate(LastPath):
        rowS = ''
        colS = ''
        #每条路径倒序遍历
        for i in range(len(aPath) -1,0,-1):
            #方向从左边来
            if aPath[i][0] == aPath[i - 1][0]:
                rowS += senquence1[aPath[i - 1][1] - 1]
                colS += '-'
            #方向从上面来
            elif aPath[i][1] == aPath[i - 1][1]:
                colS += senquence2[aPath[i - 1][0] -1]
                rowS += '-'
            #方向从左上来
            else:
                rowS += senquence1[aPath[i - 1][1] - 1]
                colS += senquence2[aPath[i - 1][0] - 1]
        #依次输出每条路的序列比对结果
        result += "------------------比对结果"+str(k+1)+"------------------\n"
        result += "s1:"+str(rowS)+"\n"
        result += "s2:" + str(colS) + "\n"
    return result



def localAlignment(replace,substitutionMatrix,senquence1,senquence2,gap):
    # ------------------比对结果1 - -----------------
    # s1: TCCGA
    # s2: T - CGA
    # ------------------比对结果2 - -----------------
    # s1: TCCGA
    # s2: TC - GA
    start_time = time.time()
    scoreMatrix,pathMatrix = LocalScoreMatrix(replace,substitutionMatrix,senquence1,senquence2,gap)
    indexs = np.array(np.where(scoreMatrix == np.max(scoreMatrix)))
    num = indexs.shape[1]
    allPaths = []
    for i in range(num):
        AllLocalPaths = []
        OnePath = []
        FindLocalPath(indexs[0][i], indexs[1][i], pathMatrix, OnePath, AllLocalPaths)
        allPaths+=AllLocalPaths
    end_time = time.time()
    runtime = end_time - start_time
    result = ShowContrastResult(allPaths, senquence1, senquence2)
    ShowPaths(senquence1, senquence2, scoreMatrix,pathMatrix,allPaths)
    return result,scoreMatrix,pathMatrix,allPaths,runtime

def GlobalAlignment(replace,substitutionMatrix,senquence1,senquence2,gap):
    m = len(senquence2)
    n = len(senquence1)
    start_time = time.time()
    scoreMatrix, pathMatrix = GlobalScoreMatrix(replace, substitutionMatrix, senquence1, senquence2, gap)
    end_time = time.time()
    runtime = end_time - start_time
    allPaths = []
    OnePath = []
    FindGlobalPath(m, n, pathMatrix, OnePath, allPaths)
    end_time = time.time()
    runtime = end_time - start_time
    result = ShowContrastResult(allPaths, senquence1, senquence2)
    ShowPaths(senquence1, senquence2, scoreMatrix, pathMatrix, allPaths)
    return result,scoreMatrix,pathMatrix,allPaths,runtime

if __name__ == "__main__":
    substitutionMatrixs = [
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ],
        [
            [1, -5, -5, -1],
            [-5, 1, -1, -5],
            [-5, -1, 1, -5],
            [-1, -5, -5, 1]
        ],
        [
            [5, -4, -4, -4],
            [-4, 5, -4, -4],
            [-4, -4, 5, -4],
            [-4, -4, -4, 5]
        ]
    ]

    # 将碱基转换为集合下标
    replace = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    senquence1="atccga".upper()
    senquence2="tcga".upper()
    # senquence1="GTGCTTGCTTTAGCAGCACATATACTAAAACTGAAATGATACAGAGAAGATTCACATGGCCCCTGCTCAAGGATGACATGCAAATTCGTGAAGTCCATATTTTT".upper()
    # senquence2="TTGCTTCAGCAGCACATATACTAACTTTGGAATGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGATGACACGCAAATTCATGAAACATTCCATTAAAAA".upper()
    #senquence1="TAGACTATACTCTCAGCGACCAGCTCCCGAATTCATTACTAGAGAAGTTTCTCTGAACGTGTAGAGCTGGGGGAACCAGGAGGGGGAGGTGAGCGTTTCTCCTGAGCATGAGGTGGCTTTAGGTATTGCTTCATGGCAATTTGGCAACTGTCATGTGCCATTGAAGATGGTTCTTATCTTCCTCTGAAAGAGGAAGAGGACACAGTCTGAGTGGT".upper()
    # senquence2="AAGACTATACTCTCAGGAATCATTTCTATAGTTTTTACTAGAGAAATTTCTCTGAATGTGTAGAGCACTGGGAACTAAGAGGAGGAGCTGCGTTCTCCTCTGAGCATGAAGGGAGCCCTTGGTGTTACTTCTCTGCAACTGCCATTTGCCATTGATGATCGTTCTTCTCTTCCTCTGCAGAGTAAGAGGGTACAGGATGCAGTCTGAGTGAT".upper()
    # #senquence1="CAATGTTGGAAGTGGGGGGGGGTGTGCGGGATTACCCACATCTGACGACACGAGCTGAGGACAACCATGCACCACCTGTCACTTTGTCCCCCGAAGGGGAAGGCTCTATCTCTAGAGTTGTCAAAGGATGTCAAGATTTGGTAAGGTTCTTCGCGTTGCTTCGAATTAAACCACATGCTCCACCGCTTGTGCGGGTCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGCGGAGTGCTTAATGCGTTAGCTGCAGCACTAAGGGGCGGAAACCCCCTAACACTTAACACTCATCGTTTACGGCGTGGACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATCTCTGCGCATTTCACCGCTACACATGGAATTCCACTTTCCTCTTCTGCACTCAAGTTTTCCAGTTTCCAATGACCCTCCACGGTTGAGCCGTGGGCTTTCACATCAAACTTAAAAAACCGCCTACGCGCGCTTTACGCCCAATAAATCCGGATAACGCTTGCAA".upper()
    # senquence2="CGCCGAGAGGTCAATGGGGGGGTCGCTGGTCGGGATTAACCCAACATCTCACGACACGAGCTGACGACAACCATGCACCACCTGTCACTTTGTCCCCCGAAGGGGAAGGCTCTATCTCTAGAGTTGTCAAAGGAAGTCAAGATTTGGTAAGGTTCTTCGCGTTGCTTCGAATTAAACCACATGCTCCACCGCTTGTGCGGATCCCCATCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGCGGAGTGCTTAATGCGTTAGCTGCAGCACTAAGGGGCGGAAACCCCCTAACACTTAACACTCATCGTTTACGGCGTGGACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCTCACATCAGCGTCAGTTACAGACCGGAAAGTCGCCTTCACCACTGGTGTTCCTCCATATCTCTGCCCATTTCACCGCTACACATGGAATTCCACTTTCCTCTTCTGCACTCAAGTTTTCCAGTTTCCAATGACCCTCCACGGTTGAGCCGTGGCCTTTCACATCAAACTTAAAAAACCGCCTACGCGCGCTTTACGCCCAAAAAATCCGGATAACGCTTGC".upper()
     # localAlignment(replace,substitutionMatrixs[0],senquence1,senquence2,-1)
    result,scoreMatrix,pathMatrix,allPaths,runtime= GlobalAlignment(replace,substitutionMatrixs[0],senquence1,senquence2,-1)

    print(ShowContrastResult(allPaths, senquence1, senquence2))
    print(runtime)



# start_time = time.time()
# scoreMatrix,pathMatrix = GlobalScoreMatrix(replace,substitutionMatrixs[0],senquence1,senquence2,-1)
# end_time = time.time()
# print(end_time - start_time)
# start_time = time.time()
# scoreMatrix,pathMatrix = oldGlobalScoreMatrix(replace,substitutionMatrixs[0],senquence1,senquence2,-1)
# end_time = time.time()
# print(end_time - start_time)
# OnePath=[]
# AllGlobalPath=[]
# FindGlobalPath(len(senquence2),len(senquence1),pathMatrix,OnePath,AllGlobalPath)
# #ShowContrastResult(AllGlobalPath,senquence1,senquence2)
# #ShowPaths(senquence1, senquence2, scoreMatrix,pathMatrix)