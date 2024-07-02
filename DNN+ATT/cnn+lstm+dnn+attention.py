# pip install openpyxl -i https://pypi.tuna.tsinghua.edu.cn/simple/
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch
from torch import nn
import torch.utils.data as data
import torch.nn.functional as F
from torch import tensor
import torch.utils.data as Data
import math
from matplotlib import pyplot
from datetime import datetime, timedelta
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import math
import warnings
from collections import Counter

# 设置随机参数：保证实验结果可以重复
SEED = 1
import random

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)  # 适用于显卡训练
torch.cuda.manual_seed_all(SEED)  # 适用于多显卡训练
from torch.backends import cudnn

cudnn.benchmark = False
cudnn.deterministic = True

warnings.filterwarnings("ignore")
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
data=pd.read_csv("datanew.csv",encoding='gbk')
data = data.sample(frac=1)
data.dropna(axis=0, how='any')
# data = data.sample(frac=1)
print(data.columns)
# ['加速踏板位置1', '驾驶员要求的发动机百分转矩', '实际发动机百分转矩', '发动机转速', 'require_per_torque',
#        '实际换档比率', '当前档位', '油温传感器1', '输出转速', '选择的档位'],
# ['v.id', 'on road old', 'on road now', 'years', 'km', 'rating','condition', 'economy', 'top speed', 'hp', 'torque', 'current price']
data_x = data[['加速踏板位置1', '驾驶员要求的发动机百分转矩', '实际发动机百分转矩', '发动机转速', 'require_per_torque','实际换档比率', '当前档位', '油温传感器1', '输出转速']].values
data_y = data['选择的档位'].values
temp_date_y=[]
for i in data_y:
    temp_date_y.append(i-1)
data_y=temp_date_y
print(Counter(data_y))
# 数据异常值归一化处理
# d = []
# for i in data_y:
#     d.append(i / 351318)
# data_y = d
# print(len(data_y))
# 四个数据划分为一组 用前三个预测后一个

len_sequence = 3
# 只改这一个参数就可以 序列的长度
data_4_x = []
data_4_y = []


for i in range(0,len(data_y), 1):
    temp_=[]
    temp_.append(data_x[i])
    temp_.append(data_x[i])
    temp_.append(data_x[i])
    data_4_x.append(temp_)
    data_4_y.append(data_y[i])

print(len(data_4_x), len(data_4_y))

x_train, x_test, y_train, y_test = train_test_split(np.array(data_4_x), np.array(data_4_y), test_size=0.2)

print(x_train.shape)
class DataSet(Data.Dataset):
    def __init__(self, data_inputs, data_targets):
        self.inputs = torch.FloatTensor(data_inputs)
        self.label = torch.FloatTensor(data_targets)

    def __getitem__(self, index):
        return self.inputs[index], self.label[index]

    def __len__(self):
        return len(self.inputs)

Batch_Size = 16 #
DataSet = DataSet(np.array(x_train), list(y_train))
train_size = int(len(x_train) * 0.8)
test_size = len(y_train) - train_size
train_dataset, test_dataset = torch.utils.data.random_split(DataSet, [train_size, test_size])
TrainDataLoader = Data.DataLoader(train_dataset, batch_size=Batch_Size, shuffle=True, drop_last=True)
TestDataLoader = Data.DataLoader(test_dataset, batch_size=Batch_Size, shuffle=True, drop_last=True)


class CNN_LSTM_ATT_DNN_Net(nn.Module):
    def __init__(self):
        # 模型是cnn + lstm + lstm + Dense
        super(CNN_LSTM_ATT_DNN_Net, self).__init__()
        # 初始参数-------
        self.input_size = 8
        # LSTM
        self.cell_LSTM = nn.LSTM(input_size=self.input_size, hidden_size=self.input_size, num_layers=2,
                                 batch_first=True)
        # lstm输入：input: shape = [seq_length, batch_size, input_size]的张量
        # lstm输出：output.shape = [seq_length, batch_size, num_directions * hidden_size]

        # cnn conv输入数据格式:[batch, channel(=1), sequence_length, embedding_size]
        self.conv = nn.Sequential(nn.Conv2d(1, 3, (2, 2)),  # 卷积核大小为2*2
                                  # nn.Conv2d(in_channels,#输入通道数 out_channels,#输出通道数 kernel_size#卷积核大小 )
                                  nn.ReLU(), nn.MaxPool2d((2, 1)), )
        #                sequence_length - 1

        # attention 机制
        self.Linear_k = nn.Linear(self.input_size, 32)
        self.Linear_q = nn.Linear(self.input_size, 32)
        self.Linear_v = nn.Linear(self.input_size, 32)
        self.multihead_attn = nn.MultiheadAttention(embed_dim=32, num_heads=2)

        # Dense层 DNN
        self.linear_1 = nn.Linear(32, 1)
        self.linear_2 = nn.Linear(3, 4)

        # 激活函数
        self.relu = F.relu

    def forward(self, x,Batch_Size):  # x shape: (batch_size, seq_len, input_size)
        # 输如数据的格式 :torch.Size([1, 1, 128])
        # 模型是 cnn+lstm+lstm+dense
        # 处理数据：
        x=x.unsqueeze(1)
        # 开始卷积层-----------------------
        # conv输入数据格式:[batch, channel(=1), sequence_length, embedding_size]
        x = self.conv(x)
        x=x.squeeze(2)
        x = x.transpose(0, 1)

        # 开始 LSTM层---------------------
        # lstm输入：input: shape = [seq_length, batch_size,embedding_size]的张量
        # print("x.shape",x.shape)
        out, _ = self.cell_LSTM(x)
        # print("200LSTM输出,out.shape--",out.shape)# shape = [seq_length, batch_size,embedding_size]的张量 torch.Size([4, 1, 31])
        # 开始 多头attention机制------------
        att_k = self.Linear_k(out)
        att_q = self.Linear_q(out)
        att_v = self.Linear_v(out)
        # print("att_k.shape",att_k.shape)# att_k.shape torch.Size([4, 1, 32])
        attn_output, attn_output_weights = self.multihead_attn(att_k, att_q, att_v)
        # print("attn_output.shape",attn_output.shape)#att_k.shape torch.Size([4, 1, 32])
        # out = attn_output.view(Batch_Size,3,9)
        out = attn_output.transpose(0,1)

        # Dense层 DNN---------------------
        out = self.relu(self.linear_1(out))
        out = out.squeeze(2)
        out = self.linear_2(out)
        return out
model = CNN_LSTM_ATT_DNN_Net( ).to(device)  # 3 表示Sequence_length  transformer 输入数据 序列的长度
def _test():
    with torch.no_grad():
        val_epoch_loss = []
        # for i in range(0, len(x_test),batch):# batch是 1 测试用1测试就行
        for index, (inputs, targets) in enumerate(TrainDataLoader):
            # inputs = x_test[i:i+batch]
            # targets = y_test[i:i+batch]
            # if len(inputs) == batch:  # 最后一个batch可能不足长度 舍弃
            inputs = torch.tensor(inputs).to(device)
            targets = torch.tensor(targets).to(device)
            inputs = inputs.float()
            targets = targets.float()
            outputs = model(inputs,Batch_Size)
            outputs = torch.tensor(outputs, dtype=torch.float)
            targets = torch.tensor(targets, dtype=torch.long)
            loss = criterion(outputs, targets)
            val_epoch_loss.append(loss.item())
    return np.mean(val_epoch_loss)

epochs = 100  # 100 200 500 1000
optimizer = torch.optim.Adagrad(model.parameters(), lr=0.0001)  # 0.03  0.05   0.1 # 优化器 SGD和Adagrad
criterion = torch.nn.CrossEntropyLoss().to(device)

val_loss = []
train_loss = []
best_test_loss = 1000000000000
for epoch in tqdm(range(epochs)):
    train_epoch_loss = []
    # for i in range(0, len(x_train),batch):# batch是 1
    for index, (inputs, targets) in enumerate(TrainDataLoader):
        inputs = torch.tensor(inputs).to(device)
        targets = torch.tensor(targets).to(device)
        inputs = inputs.float()
        targets = targets.float()
        # print(inputs.shape)
        outputs = model(inputs,Batch_Size)
        outputs = torch.tensor(outputs, dtype=torch.float)
        targets = torch.tensor(targets, dtype=torch.long)
        # print(outputs)
        # print(targets)
        loss = criterion(outputs, targets)
        print("loss:", loss)
        loss.requires_grad_(True)
        loss.backward()
        optimizer.step()
        train_epoch_loss.append(loss.item())
    train_loss.append(np.mean(train_epoch_loss))
    val_epoch_loss = _test()
    val_loss.append(val_epoch_loss)
    print("epoch:", epoch, "train_epoch_loss:", train_epoch_loss, "val_epoch_loss:", val_epoch_loss)
    # 保存下来最好的模型：
    if val_epoch_loss < best_test_loss:
        best_test_loss = val_epoch_loss
        best_model = model
        print("best_test_loss -------------------------------------------------", best_test_loss)
        torch.save(best_model.state_dict(), 'best_Transformer_trainModel.pth')

# 画一下loss图
fig = plt.figure(facecolor='white', figsize=(10, 7))
plt.xlabel('X')
plt.ylabel('Y')
plt.xlim(xmax=len(val_loss), xmin=0)
plt.ylim(ymax=max(max(train_loss), max(val_loss)), ymin=0)
# 画两条（0-9）的坐标轴并设置轴标签x，y
x1 = [i for i in range(0, len(train_loss), 1)]  # 随机产生300个平均值为2，方差为1.2的浮点数，即第一簇点的x轴坐标
y1 = val_loss  # 随机产生300个平均值为2，方差为1.2的浮点数，即第一簇点的y轴坐标
x2 = [i for i in range(0, len(train_loss), 1)]
y2 = train_loss
colors1 = '#00CED4'  # 点的颜色
colors2 = '#DC143C'
area = np.pi * 4 ** 1  # 点面积
# 画散点图
plt.scatter(x1, y1, s=area, c=colors1, alpha=0.4, label='val_loss')
plt.scatter(x2, y2, s=area, c=colors2, alpha=0.4, label='train_loss')
plt.legend()
plt.savefig("loss图.png")
plt.show()

# 加载模型预测------
model = CNN_LSTM_ATT_DNN_Net().to(device)
model.load_state_dict(torch.load('best_Transformer_trainModel.pth'))
model.to(device)
model.eval()
# 在对模型进行评估时，应该配合使用with torch.no_grad() 与 model.eval()：
y_pred = []
y_true = []
with torch.no_grad():
    with torch.no_grad():
        val_epoch_loss = []
        for index, (inputs, targets) in enumerate(TrainDataLoader):
            inputs = torch.tensor(inputs).to(device)
            targets = torch.tensor(targets).to(device)
            inputs = inputs.float()
            targets = targets.float()
            tgt_in = torch.rand((Batch_Size, len_sequence,9))
            outputs = model(inputs, tgt_in)
            for t in np.array(outputs):
                t = np.argmax(t)
                y_pred.append(t)
            y_true.extend(targets)

y_true = np.array(y_true)
y_pred = np.array(y_pred)
d_y_true = []
d_y_pred = []

from sklearn.metrics import mean_squared_error, mean_absolute_error  # 评价指标
from sklearn.metrics import f1_score, roc_auc_score, accuracy_score

print(accuracy_score(y_pred, y_true))
from sklearn.metrics import precision_score
print(precision_score(y_pred, y_true,average='weighted'))
