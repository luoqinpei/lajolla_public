import pickle
from plyfile import PlyData

# 读取 .ply 文件
ply_data = PlyData.read('bunny.ply')

# 序列化为文件
with open('my_disney.serialized', 'wb') as f:
    pickle.dump(ply_data, f)