import csv
from operator import itemgetter
from scipy import optimize, math
import  numpy as np
import os
from flask import render_template
from flask import Flask, request
from datetime import *

v = 10
#步骤1，读读降温数据（确保数据是从大到小）
def readCsv(filename):
    with open(filename) as f:
        table = []
        for line in f:
            col = line.split(',')
            col[0] = float(col[0])
            col[1] = float(col[1])
            table.append(col)
        table_sorted = sorted(table, key=itemgetter(0),reverse=False) #对降温数据进行排序，默认为升序排列
        x0 = []
        y0 = []
        for row in table_sorted:
            x0.append(row[0])
            y0.append(row[1])
    return x0,y0

def f(x,a,b):
    return a*x+b

# 步骤2，进行基线修正，确定a，b两点，返回修正后的坐标，存在x，y中
#降温数据从温度高的地方先修正
def correct(filename):
    x0, y0 =readCsv(filename)
    maxy = y0.index(max(y0))  # 最低温度点
    x = []
    y = []
    for i in range(len(y0)):
        y.append(y0[i] / y0[maxy] * 100)
        x.append(x0[i])

    # 用来寻找终点b
    count = 0
    for i in range(len(x)):
        count = count + 1  # count为5℃所包含的数据点
        if x[i] - x[0] >= 5:
            break

    # 剔除升温数据温度最后5℃包含的数据点
    for i in range(maxy, len(x) - count, 1):
        xx = [x[i], x[i + count]]
        yy = [y[i], y[i + count]]
        a, b = optimize.curve_fit(f, xx, yy)[0]
        angle = math.atan(a) * 180 / math.pi  # 将斜率转换为角度
        # print("%f-%f:%f" % (x[i], x[i+count], angle))
        if angle < 5 and angle > 0:  # 如果角度小于5，取该点为终点
            temp1 = i
            break

    # 用来寻找起点a
    # 去掉降温数据温度开始5℃的数据点
    for i in range(maxy, count + 1, -1):
        xx = [x[i - count], x[i]]
        yy = [y[i - count], y[i]]
        a, b = optimize.curve_fit(f, xx, yy)[0]
        angle = math.atan(a) * 180 / math.pi
        # print("%f-%f:%f" % (x[i - count], x[i], angle))
        if angle < 4:  # 如果角度小于5，取该点为终点
            temp = i
            break

    # # 在ab之间的点，存入x,y中
    # print('correct')
    # print(x[temp], y[temp])
    # print(x[temp1], y[temp1])

    x1 = []
    y1 = []
    for i in range(temp, temp1 + 1):
        x1.append(x[i])
        y1.append(y[i])
    return x1, y1

# 步骤3，求相对结晶度Xt，存在y4中，降序排列(因为求的基线不是近似平行的，所以求直线积分，然后减去）
def caculateXt(filename):
    x, y = correct(filename)
    # 拟合后的直线方程
    xx = [x[0], x[len(x) - 1]]
    yy = [y[0], y[len(y) - 1]]
    a,b = optimize.curve_fit(f, xx, yy)[0]
    y1=[]
    for i in range(len(x)):
        y1.append(a*x[i]+b)
    area=np.trapz(y,x)-np.trapz(y1,x)
    x11=[]
    y11=[]
    y12=[]
    Xt=[]
    for i in range(len(x)):
        x11.append(x[i])
        y11.append(y[i])
        y12.append(y1[i])
        area1=np.trapz(y11,x11)-np.trapz(y12,x11)
        area2=abs(area)-abs(area1)
        Xt.append(area2/area)   #相对结晶度从1到0存储
    return Xt

#步骤4，将温度变化从b-a转化成从b-a所用的时间，降序排列的，存在x2中,
# !!!!注意，因为下面用的是部分时间与总时间的比例，所以降温速度为多少，并不影响最后结果
def temperatureToTime(filename):
    x=correct(filename)[0]
    # x2中存储温度变化时间，降序（b-a到b-b）
    x2 = []
    for k in range(len(x)):
        x2.append(abs(x[len(x) - 1] - x[k]) / v)
    return x2

#步骤5，计算经转化后的t，Xt，存在x5，y5中
    # 横坐标为In(t/t总)，存在x5中
    # 从大到小，纵坐标为ln(-ln(1-Xt))，存在y5中
def changeTAndXt(filename):
    x2=temperatureToTime(filename)
    y4=caculateXt(filename)
    x5 = []
    y5 = []
    totalTime = x2[0]
    for i in range(1, len(x2) - 1):  # x2中最后一个元素为0，不能做分母，排除
        tr = x2[i] / totalTime
        x5.append(math.log(tr))
    for i in range(1, len(y4) - 1):  # y4中第一个元素为1，最后一个元素为0，不能放在log里面，排除
        y5.append(math.log(-(math.log(1 - y4[i]))))
    return x5,y5

#步骤6，进行曲线拟合，求交点
#6.1待拟合的两条线为直线
def f(x,a,b):
    return a*x+b

# 6.2进行-5到-3之间的直线拟合
#6.2.1找到符合条件的第一条直线的离散点
def getFirstProfitPoint(filename):
    x5 = changeTAndXt(filename)[0]
    y5 = changeTAndXt(filename)[1]
    x6 = []
    y6 = []
    for i in range(len(x5)):
        if x5[i] >= -5 and x5[i] <= -3:
            x6.append(x5[i])
            y6.append(y5[i])
    return x6,y6
#6.2.2根据离散点进行线性拟合，求出第一条直线方程
def getFirstLine(filename):
    x6,y6=getFirstProfitPoint(filename)
    a1, b1 = optimize.curve_fit(f, x6, y6)[0]
    return a1,b1
#6.2.3根据拟合求出来的第一条直线方程，取出若干对x，y
def getFirstLinePoint(filename):
    a1,b1=getFirstLine(filename)
    x = np.arange(-5,0, 0.01)
    y = a1 * x + b1
    return x,y

# 6.3进行-1.6到-1.1之间的直线拟合
#6.3.1找到符合条件的第二条直线的离散点
def getSecondProfitPoint(filename):
    x5 = changeTAndXt(filename)[0]
    y5 = changeTAndXt(filename)[1]
    x7 = []
    y7 = []
    for i in range(len(x5)):
        if x5[i] >= -1.6 and x5[i] <= -1.1:  # 因为x5中是从小到大排序的，都是小于0的，因此第一个找到的点即近似分离点
            x7.append(x5[i])
            y7.append(y5[i])
    return x7,y7
#6.3.2根据离散点进行线性拟合，求出第二条直线方程
def getSecondLine(filename):
    x7,y7=getSecondProfitPoint(filename)
    a2, b2 = optimize.curve_fit(f, x7, y7)[0]
    return a2,b2
#6.3.3根据拟合求出来的第二条直线方程，取出若干对x，y
def getSecondLinePoint(filename):
    a2,b2=getSecondLine(filename)
    x = np.arange(-4, 0, 0.01)
    y = a2 * x + b2
    return x,y

# 6.4求[-5,-3]和[-1.6,-1.1]两直线的交点
def getPointOfIntersection(filename):
    a1,b1=getFirstLine(filename)
    a2,b2=getSecondLine(filename)
    x = (b1 - b2) / (a2 - a1)
    y = a1 * x + b1
    return x,y

#步骤7，根据交点求出Xt1，Xt2
#步骤7.1，过交点做垂直于x轴的线，交于原曲线，求出此交点的y值（x0从大到小，y0从小到大）
def getY2(filename):
    x0,y0=changeTAndXt(filename)
    x,y=getPointOfIntersection(filename)
    for i in range(len(x0)):
        if x==x0[i]:  #如果存在x的值和x0[i]相等，则直接令x对的y值为x0[i]这一点的纵坐标值
            y2=y0[i]
        elif x>x0[i]:  #否则，去x对应的y值为介于x两边的纵坐标值的平均值
            y2=(y0[i]+y0[i-1])/2
            break
    return y2
#步骤7.2，由y2求出Xt2
def getXt2(filename):
    y=getY2(filename)
    Xt2 = 1 - math.exp(-math.exp(y))
    return Xt2
#步骤7.3，由步骤6中的交点中的y1求出Xt1
def getXt1(filename):
    y=getPointOfIntersection(filename)[1]
    Xt1=1 - math.exp(-math.exp(y))
    return Xt1


##步骤二，求基线修正后的升温曲线以及Tm

#步骤1，对升温数据进行基线修正
def correctHot(filename):
    x0,y0 =readCsv(filename)
    miny = y0.index(min(y0)) #最低温度点
    x=[]
    y=[]
    for i in range(len(y0)):
        y.append(-y0[i]/y0[miny]*100)
        x.append(x0[i])

    # 用来寻找终点b
    count=0
    for i in range(len(x)):
        count = count + 1  #count为5℃所包含的数据点
        if x[i] - x[0] >=5:
            break

    #剔除升温数据温度最后5℃包含的数据点
    for i in range(miny, len(x) - count, 1):
        xx = [x[i], x[i + count]]
        yy = [y[i], y[i + count]]
        a, b = optimize.curve_fit(f, xx, yy)[0]
        angle = math.atan(a) * 180 / math.pi  #将斜率转换为角度
        # print("%f-%f:%f" % (x[i], x[i+count], angle))
        if angle < 5:      #如果角度小于5，取该点为终点
            temp1 = i
            break

    #用来寻找起点a
    ang=[]
    # 去掉升温数据温度开始5℃的数据点
    for i in range(miny,count+1,-1):
        xx = [x[i-count],x[i]]
        yy = [y[i-count],y[i]]
        a, b = optimize.curve_fit(f, xx, yy)[0]
        angle=math.atan(a)*180/math.pi
        ang.append(angle)
        # print("%f-%f:%f"%(x[i-count],x[i],angle))
    maxAng=max(ang)
    for i in range(len(ang)):
        if ang[i]==maxAng:
            temp=miny-i   #两个数据之间的位置换算
            break
    # 在ab之间的点，存入x,y中
    # print('correct')
    # print(x[temp], y[temp])
    # print(x[temp1], y[temp1])

    x1 = []
    y1 = []
    for i in range(temp, temp1 + 1):
        x1.append(x[i])
        y1.append(y[i])
    return x1, y1

#步骤2，求熔融焓Hm(由前台自己输入)

#步骤3，求Tm（计算Mn和Mw时使用的不同）
#步骤3.1，计算基线修正后的温度的平均温度（Mn）
def getAverageTm(filename):
    x=correctHot(filename)[0]
    sum=0
    for i in range(len(x)):
        sum=sum+x[i]
    Tm=sum/len(x)
    return Tm

#3.2，求出纵坐标乘以横坐标除以纵坐标之和作为平均温度（Mw）
def getTmOfMw(filename):
    x,y=correctHot(filename)
    sum=0
    sum1=0
    for i in range(len(x)):
        sum=x[i]*abs(y[i])+sum
        sum1=abs(y[i])+sum1
    tm=sum/sum1
    return tm

##步骤三，计算Mn和Mw
#步骤1，计算Mn，其中deltaHm为熔融焓（不同分子量，deltaHm不同），Tm为平均温度（计算Mw和Mn使用温度不同）
#filename3为升温数据的,deltaHm为熔融焓，用户自己输入
def getMn(Xt2,filename3,Lu, Mo, pc, deltaHm0, sigmae, Tm0,deltaHm):
    deltaHm2=deltaHm*deltaHm
    Tm=getAverageTm(filename3)
    fz=2*sigmae*deltaHm0*Mo  #分子
    fm = deltaHm2 * Xt2 * Lu * pc * (1 - Tm / Tm0)
    Mn=fz/fm
    return Mn

#步骤2，计算Mw，Tm为平均温度
def getMw(Xt1,Xt2,filename3,Lu, Mo, pc, deltaHm0, sigmae, Tm0,deltaHm):
    Tm =getTmOfMw(filename3)
    deltaHm2 =deltaHm*deltaHm
    fz = 2 * sigmae * deltaHm0 * Mo  # 分子
    fm = deltaHm2 * Xt1 * Lu * pc * (1 - Tm / Tm0)
    fm1= deltaHm2 * Xt2 * Lu * pc * (1 - Tm / Tm0)
    Mi1=fz/fm
    Mi2=fz/fm1
    Mw=(Mi1+Mi2)/2
    return Mw

def getPDI(Mn,Mw):
    PDI=Mw/Mn
    return PDI

#步骤3，计算Mn，Mw和PDI
def getAlldata(filename2,filename3,Lu, Mo, pc, deltaHm0, sigmae, Tm0,deltaHm):
    xt1=getXt1(filename2)
    xt2=getXt2(filename2)
    Mn=getMn(xt2,filename3, Lu, Mo, pc, deltaHm0, sigmae, Tm0,deltaHm)
    Mw=getMw(xt1,xt2, filename3, Lu, Mo, pc, deltaHm0, sigmae, Tm0,deltaHm)
    PDI=getPDI(Mn,Mw)
    return Mn,Mw,PDI

##步骤四，后台和前端的交互

ALLOWED_EXTENSIONS = set(['txt','csv'])
#python3自动生成文件名
tmp = 'static/file'
curdir = os.path.abspath('.')  # 获得当前工作目录,如果在加一个点，是获得当前目录的父目录
UPLOAD_FOLDER = curdir + os.path.sep + tmp + os.path.sep  # 该路径为当前文件夹拼接windows下文件分隔符再拼接'tmp'文件夹，再拼接文件分隔符
if os.path.exists(UPLOAD_FOLDER) == False:
    os.makedirs(UPLOAD_FOLDER)

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

#高分子数据(Lu,Mo,pc,deltaHm0,sigmae,Tm0)
#(结构单元长度nm，结构单元分子量，晶区密度，完美熔融焓，表面能，平衡熔点k)
#JBX:聚丙烯  PEO:聚氧化乙烯
JBXData=[0.2167,42,0.936,209.2,30,459-273.15]
PEOData=[0.2783,44,1.228,220,31,69]

def getData(arr):
    Lu=float(arr[0])
    Mo=float(arr[1])
    pc=float(arr[2])
    deltaHm0=float(arr[3])
    sigmae=float(arr[4])
    Tm0=float(arr[5])
    return Lu,Mo,pc,deltaHm0,sigmae,Tm0


#用来判断上传的文件是否是csv文件
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

#生成唯一的文件名
def randomfile(filename):
    nowTime =datetime.now().strftime("%Y%m%d%H%M%S")  # 生成当前的时间
    uniqueNum=str(nowTime)
    uniqueFile=uniqueNum+filename
    return uniqueFile

#下载单个文件,返回改文件的路径
def upload_single_file(label):
    file = request.files[label]
    if file and allowed_file(file.filename):
        filename = randomfile(file.filename)
        uploadpath=os.path.join(app.config['UPLOAD_FOLDER'], filename).replace('\\','/')
        file.save(uploadpath)
    return uploadpath

#返回所有文件的路径
def upload_all_file():
    uploadpath2 = upload_single_file('file2')
    uploadpath3 = upload_single_file('file3')
    return uploadpath2,uploadpath3

#计算所得到的值，并返回
def caculate():
    path2,path3=upload_all_file()
    # v=int(request.form.get('velocity'))  #获取前台的v值
    deltaHm=float(request.form.get('Hm'))
    if request.form.get("select") == 'JuBingXi':  #获取前台select中选中的value值
        Lu, Mo, pc, deltaHm0, sigmae, Tm0 = getData(JBXData)
    if request.form.get("select") == 'JuYangHuaYiXi':
        Lu, Mo, pc, deltaHm0, sigmae, Tm0 = getData(JBXData)
    Mn,Mw,PDI=getAlldata(path2, path3, Lu, Mo, pc, deltaHm0, sigmae, Tm0,deltaHm)
    return Mn,Mw,PDI

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        Mn,Mw,PDI=caculate()
        data={
            'Mn':Mn,
            'Mw':Mw,
            'PDI':PDI
        }
        return render_template('result.html',data=data)
    return render_template('upload.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0')